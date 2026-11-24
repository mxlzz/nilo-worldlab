"""Verbatim B.2b reduced-owner functions; no native solver imported."""
import hashlib,json
import numpy as np
import xarray as xr
import probe_e5_h31_owner_contract as c

EPS=.01
CP=800.
NAMES={"soil_mass_kg","grain_mass_kg","E_rest","E_roi","E_grain"}
UNITS={"soil_mass_kg":"kg","grain_mass_kg":"kg","E_rest":"J","E_roi":"J","E_grain":"J"}

def req(x,msg):
    if not bool(x):raise ValueError(msg)

def digest(ds):
    h=hashlib.sha256()
    for name in sorted(ds.variables):
        v=ds[name];h.update(name.encode()+str(v.dims).encode()+str(v.dtype).encode()+v.values.tobytes()+json.dumps(dict(v.attrs),sort_keys=True).encode())
    h.update(json.dumps({k:c.normal(v) for k,v in ds.attrs.items()},sort_keys=True).encode())
    return h.hexdigest()

def total(ds):
    return float(ds.E_rest.sum()+ds.E_roi.sum()+ds.E_grain.sum())+float(ds.attrs["E_export_J"])

def fraction(ds):
    return float(ds.attrs["roi_fraction"])

def capacities(ds):
    m=ds.soil_mass_kg.values
    n=int(ds.attrs["roi_node"]);f=fraction(ds)
    rest=m.copy();rest[n]*=1-f
    roi=np.full(ds.E_roi.shape,m[n]*f/ds.E_roi.size)
    return rest,roi

def aggregate(ds):
    e=ds.E_rest.values.copy();e[int(ds.attrs["roi_node"])]+=float(ds.E_roi.sum())
    return e

def validate(ds,expected_m,expected_ids,expected_grain,expected_step,expected_n,expected_f,expected_shape,expected_initial,expected_initial_mass,expected_incident):
    req(set(ds.data_vars)==NAMES and set(ds.coords)=={"node","grain","ty","tx"},"duplicate or missing stock")
    req(set(ds.attrs)=={"schema","step","roi_node","roi_fraction","E_export_J","E_incident_J","initial_energy_J","initial_mass_kg","mass_export_kg"},"attrs")
    req(ds.attrs["schema"]=="P1B2b-reduced-v1" and ds.attrs["step"]==expected_step,"schema/step")
    req(ds.attrs["roi_node"]==expected_n and ds.attrs["roi_fraction"]==expected_f and ds.E_roi.shape==expected_shape,"partition geometry")
    req(ds.attrs["initial_energy_J"]==expected_initial and ds.attrs["initial_mass_kg"]==expected_initial_mass and ds.attrs["E_incident_J"]==expected_incident,"ledger metadata")
    req(np.array_equal(ds.node,np.arange(len(expected_m))) and np.array_equal(ds.grain,expected_ids) and np.array_equal(ds.ty,np.arange(expected_shape[0])) and np.array_equal(ds.tx,np.arange(expected_shape[1])),"identity")
    req(np.array_equal(ds.soil_mass_kg,expected_m) and np.array_equal(ds.grain_mass_kg,expected_grain),"material stock changed")
    for name in NAMES:
        req(ds[name].attrs=={"units":UNITS[name]},"unit "+name)
        req(np.isfinite(ds[name]).all(),"nonfinite "+name)
    req(ds.soil_mass_kg.dims==("node",) and ds.E_rest.dims==("node",) and ds.E_roi.dims==("ty","tx") and ds.E_grain.dims==("grain",) and ds.grain_mass_kg.dims==("grain",),"dimensions")
    mr,mi=capacities(ds)
    req(np.all(mr>=0) and np.all(mi>=0) and np.all(expected_grain>0),"negative mass")
    reconstructed=mr.copy();reconstructed[expected_n]+=float(mi.sum())
    req(np.all(np.abs(reconstructed-expected_m)<=np.maximum(1e-6,1e-9*expected_m)),"partition mass")
    req(np.isfinite([ds.attrs["E_export_J"],ds.attrs["mass_export_kg"]]).all() and ds.attrs["E_export_J"]>=-1e-7 and ds.attrs["mass_export_kg"]>=0,"export")
    req(abs(float(ds.soil_mass_kg.sum()+ds.grain_mass_kg.sum())+ds.attrs["mass_export_kg"]-expected_initial_mass)<=.09216,"material balance")
    for e,m in [(ds.E_rest.values,mr),(ds.E_roi.values,mi),(ds.E_grain.values,expected_grain)]:
        req(np.all(e>=-1e-7) and np.all(e[m==0]==0),"energy positivity/zero mass")
        T=300+np.divide(e,CP*m,out=np.zeros_like(e),where=m>0)
        req(np.isfinite(T).all(),"temperature")
    req(abs(total(ds)-expected_initial-expected_incident)<=EPS,"energy balance")

def make(m,n,f,shape,specific,pattern=False,ids=None,gm=None):
    m=np.asarray(m,dtype=float);specific=np.broadcast_to(np.asarray(specific,dtype=float),m.shape)
    ids=np.arange(125,dtype=np.int32) if ids is None else np.asarray(ids,dtype=np.int32)
    gm=np.full(len(ids),.025) if gm is None else np.asarray(gm,dtype=float)
    rest=m*specific;roi_total=rest[n]*f;rest[n]*=1-f
    weights=np.ones(shape)
    if pattern:weights[:shape[0]//2,:shape[1]//2]=4.;weights[shape[0]//2:,shape[1]//2:]=.25
    roi=roi_total*weights/weights.sum()
    ds=xr.Dataset({"soil_mass_kg":("node",m),"grain_mass_kg":("grain",gm),"E_rest":("node",rest),"E_roi":(("ty","tx"),roi),"E_grain":("grain",np.zeros(len(ids)))},coords={"node":np.arange(len(m),dtype=np.int32),"grain":ids,"ty":np.arange(shape[0],dtype=np.int32),"tx":np.arange(shape[1],dtype=np.int32)},attrs={"schema":"P1B2b-reduced-v1","step":0,"roi_node":n,"roi_fraction":f,"E_export_J":0.,"E_incident_J":0.,"initial_energy_J":float(rest.sum()+roi.sum()),"initial_mass_kg":float(m.sum()+gm.sum()),"mass_export_kg":0.})
    for name,u in UNITS.items():ds[name].attrs={"units":u}
    return ds

def contract(ds,m,step=None):
    return (m,ds.grain.values.copy(),ds.grain_mass_kg.values.copy(),ds.attrs["step"] if step is None else step,int(ds.attrs["roi_node"]),fraction(ds),ds.E_roi.shape,float(ds.attrs["initial_energy_J"]),float(ds.attrs["initial_mass_kg"]),float(ds.attrs["E_incident_J"]))

def advance(ds,flow,step):
    before,after,e,d,qi,qo,rec,order,core=[np.asarray(flow[k]) for k in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")]
    original=contract(ds,before)
    validate(ds,*original)
    req(step==ds.attrs["step"]+1,"duplicate or out-of-order step")
    req(np.all(e>=0) and np.all(d>=0) and np.all(qi>=0) and np.all(qo>=0) and np.all(e<=before),"invalid flow")
    req(np.all((rec>=0)&(rec<len(rec))) and np.array_equal(np.sort(order),np.arange(len(rec))),"receiver/order")
    pos=np.empty(len(rec),int);pos[order]=np.arange(len(rec))
    req(all(rec[i]!=i and pos[i]<pos[rec[i]] for i in np.flatnonzero(qo>0)),"cyclic route")
    tol=np.maximum(1e-6,1e-9*np.maximum.reduce([before,after,e,d,qi,qo,np.ones_like(e)]))
    req(np.all(np.abs(after-before+e-d)<=tol) and np.all(np.abs((qi+e-qo-d)[core])<=tol[core]),"material flow balance")
    incoming_mass=np.zeros_like(qi);np.add.at(incoming_mass,rec,qo)
    req(np.all(np.abs(incoming_mass-qi)<=tol),"route balance")
    req(np.all(qo[~core]==0) and np.all(e[~core]==0) and np.all(d[~core]==0),"boundary stock")
    n=int(ds.attrs["roi_node"]);f=fraction(ds)
    rate=np.divide(e,before,out=np.zeros_like(e),where=before>0)
    removed_rest=ds.E_rest.values*rate
    removed_roi=ds.E_roi.values*rate[n]
    removed=removed_rest.copy();removed[n]+=float(removed_roi.sum())
    incoming=np.zeros_like(before);deposit=np.zeros_like(before);outgoing=np.zeros_like(before)
    for i in order:
        if not core[i]:continue
        available=qi[i]+e[i];enthalpy=incoming[i]+removed[i]
        if available==0:
            req(enthalpy==0 and d[i]==0 and qo[i]==0,"energy without flow")
            continue
        deposit[i]=enthalpy*(d[i]/available)
        outgoing[i]=enthalpy*(qo[i]/available)
        incoming[rec[i]]+=outgoing[i]
    exp=float(incoming[~core].sum())
    result=ds.copy(deep=True)
    result.soil_mass_kg.values[:]=after
    result.E_rest.values[:]=ds.E_rest.values-removed_rest+deposit
    result.E_rest.values[n]-=deposit[n]*f
    result.E_roi.values[:]=ds.E_roi.values-removed_roi+deposit[n]*f/ds.E_roi.size
    result.attrs.update(step=step,E_export_J=float(ds.attrs["E_export_J"])+exp,mass_export_kg=float(ds.attrs["mass_export_kg"])+float(qi[~core].sum()))
    validate(result,*contract(ds,after,step))
    req(abs(total(result)-total(ds))<=EPS,"per-operation energy balance")
    return result,{"export_J":exp,"removed_J":float(removed.sum()),"deposited_J":float(deposit.sum()),"outgoing_J":float(outgoing.sum()),"residual_J":total(result)-total(ds)}

def flows(m,e,d,qi,qo,rec,order,core):
    m,e,d,qi,qo=[np.asarray(x,dtype=float) for x in (m,e,d,qi,qo)]
    return dict(before=m,after=m-e+d,erosion=e,deposition=d,influx=qi,outflux=qo,receiver=np.array(rec),order=np.array(order),core=np.array(core,dtype=bool))

def remap(ds,size):
    source=ds.E_roi.values;old=len(source)
    req(old==source.shape[1] and size in (30,60,120),"unsupported grid")
    if size>=old:
        r=size//old;req(size==r*old,"not integer")
        dest=np.repeat(np.repeat(source,r,axis=0),r,axis=1)/(r*r)
    else:
        r=old//size;req(old==r*size,"not integer")
        dest=source.reshape(size,r,size,r).sum(axis=(1,3))
    other=ds.drop_vars("E_roi").drop_dims(["ty","tx"]).assign_coords(ty=np.arange(size,dtype=np.int32),tx=np.arange(size,dtype=np.int32))
    other["E_roi"]=(("ty","tx"),dest,{"units":"J"})
    validate(other,*contract(other,ds.soil_mass_kg.values))
    req(abs(total(ds)-total(other))<=EPS,"remap energy")
    return other

