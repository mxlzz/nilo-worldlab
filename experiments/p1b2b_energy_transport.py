"""P1-B.2b: test a disjoint material energy partition and passive sediment heat.
Reduced laboratory owner, not the full evolving P0-D owner. No thermal solver or mechanics.
"""
import hashlib, io, json, pathlib, platform, time, traceback, zipfile
from importlib.metadata import version
import numpy as np
import xarray as xr
from landlab.components import ConcentrationTrackerForSpace
import probe_e6_h33_resume_landlab as h33
import probe_e5_h31_owner_contract as c
import p1b2a_native_flux as a

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"p1b2b-output"
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

def analytical():
    report={}
    # Upstream 10 kg to downstream: 6 deposit, 4 export; upstream starts at 100 J/kg.
    fl=flows([20,10,0],[10,0,0],[0,6,0],[0,10,4],[10,4,0],[1,2,2],[0,1,2],[1,1,0])
    ds=make(fl["before"],1,.1,(60,60),[100,0,0]);key=digest(ds)
    res,m=advance(ds,fl,1)
    req(np.max(np.abs(aggregate(res)-[1000,600,0]))<=1e-7 and abs(m["export_J"]-400)<=1e-7,"10 kg analytical")
    req(digest(ds)==key,"input mutation");report["chain"]=m
    # Two disjoint D8 chains each delivering 3 kg to a material receiver.
    fl2=flows([10,10,10,10,0],[3,3,0,0,0],[0,0,3,3,0],[0,0,3,3,0],[3,3,0,0,0],[2,3,4,4,4],[0,1,2,3,4],[1,1,1,1,0])
    s=make(fl2["before"],2,.1,(60,60),[100,200,0,0,0]);z,_=advance(s,fl2,1)
    req(np.max(np.abs(aggregate(z)-[700,1400,300,600,0]))<=1e-7,"two chains");report["two_chains"]=True
    # Converging tributaries and simultaneous erosion/deposition with zero net stock change at node 2.
    f3=flows([10,10,10,0],[2,4,3,0],[0,0,3,0],[0,0,6,6],[2,4,6,0],[2,2,3,3],[0,1,2,3],[1,1,1,0])
    s=make(f3["before"],2,.1,(60,60),[100,200,50,0]);z,m=advance(s,f3,1)
    req(np.max(np.abs(aggregate(z)-[800,1200,350+1150/3,0]))<=1e-7 and abs(m["export_J"]-2300/3)<=1e-7,"convergence/net-zero mass");report["convergence"]=m
    zeros=flows([10,10,0],[0]*3,[0]*3,[0]*3,[0]*3,[1,2,2],[0,1,2],[1,1,0])
    s=make(zeros["before"],1,.1,(60,60),[1,2,0],True);z,m=advance(s,zeros,1)
    req(np.array_equal(s.E_roi,z.E_roi) and np.array_equal(s.E_rest,z.E_rest) and m["export_J"]==0,"K0");report["K0"]=True
    s=make(fl["before"],1,.1,(60,60),0);z,m=advance(s,fl,1);req(total(z)==0,"zero energy");report["zero_energy"]=True
    for pattern in (False,True):
        s=make([10,10,0],1,.1,(60,60),[0,2,0],pattern)
        for size in (30,120):
            z=remap(s,size);back=remap(z,60)
            req(np.max(np.abs(back.E_roi-s.E_roi))<=1e-7,"remap roundtrip")
            req(abs(float(back.E_roi.sum()-s.E_roi.sum()))<=1e-7,"remap integral")
        if pattern:
            mean=float(s.E_roi.mean());req(np.max(np.abs(s.E_roi-mean))>0,"coarse projection contrast")
    report["remap"]=True
    # A zero-mass ROI can exist only with zero energy; thermal activation is not tested here.
    s=make([10,0,0],1,.1,(60,60),[2,0,0]);validate(s,*contract(s,s.soil_mass_kg.values));report["zero_mass"]=True
    ds=make(fl["before"],1,.1,(60,60),[100,0,0])
    spec=contract(ds,fl["before"])
    bads={
        "energy_plus_1":lambda x:x.E_rest.values.__setitem__(0,x.E_rest.values[0]+1),
        "mass_plus_1":lambda x:x.soil_mass_kg.values.__setitem__(0,x.soil_mass_kg.values[0]+1),
        "units":lambda x:x.E_rest.attrs.__setitem__("units","kJ"),
        "identity":lambda x:x.assign_coords(grain=x.grain.values+1),
        "false_export":lambda x:x.attrs.__setitem__("E_export_J",1.),
        "duplicate_stock":lambda x:x.__setitem__("extra_mass",x.soil_mass_kg.copy()),
    }
    for label,mutate in bads.items():
        bad=ds.copy(deep=True);v=mutate(bad)
        if isinstance(v,xr.Dataset):bad=v
        try:validate(bad,*spec)
        except ValueError:report["reject_"+label]=True
        else:raise ValueError("corruption accepted "+label)
    for label,fn in [("repeat",lambda:advance(res,fl,1)),("receiver",lambda:advance(ds,{**fl,"receiver":np.array([99,2,2])},1))]:
        try:fn()
        except ValueError:report["reject_"+label]=True
        else:raise ValueError("invalid operation accepted "+label)
    p=OUT/"analytical_owner.nc";res.to_netcdf(p,engine="scipy")
    with xr.open_dataset(p,engine="scipy") as back:
        b=back.load();validate(b,*contract(res,res.soil_mass_kg.values));req(digest(b)==digest(res),"NetCDF roundtrip")
    report["roundtrip"]=True
    return report

def native_case(raw,seed,k,z,mode):
    key=f"{seed}_{k}";owner=h33.initial_owner(raw,seed,k)
    m=owner.soil_mass_kg.values.ravel().copy()
    y,x={7321:(10,10),8128:(10,16)}[seed];roi=y*c.SHAPE[1]+x
    req(38.5>=-50 and 41.5<=50 and 28.5>=-50 and 31.5<=50 and m[roi]>0,"ROI geometry")
    specific=np.ones_like(m) if mode=="uniform" else .1+np.arange(len(m))/len(m)
    ds=make(m,roi,9/10000,(60,60),specific,mode!="uniform")
    reference_E=aggregate(ds);history=[];arrays={}
    max_diff=0.;first_fail=None
    for block in range(1,19,3):
        grid,flow,solver,outlet=h33.model_from_owner(owner,k)
        concentration=np.divide(reference_E*c.RHO,m,out=np.zeros_like(m),where=m>0)
        tracker=ConcentrationTrackerForSpace(grid,phi=0.,fraction_fines=0.,settling_velocity=3.,concentration_initial=concentration,concentration_in_bedrock=0.)
        for step in range(block,block+3):
            fs={name:z[f"s{step}_{name}"] for name in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
            req(np.array_equal(fs["before"],m),"archived predecessor")
            np.copyto(grid.at_node["soil__depth"],m/(c.RHO*c.AREA))
            np.copyto(grid.at_node["topographic__elevation"],c.topography(owner).ravel())
            flow.run_one_step();tracker.start_tracking();suspended,bedload=solver.run_one_step(c.LANDLAB_DT);tracker.stop_tracking(c.LANDLAB_DT)
            # Independent native tracker receives identical material flows.
            node=grid.at_node;dt=c.LANDLAB_DT
            after=node["soil__depth"].copy()*c.RHO*c.AREA
            req(np.array_equal(after,fs["after"]),"native material mismatch")
            for name,field,factor in [("erosion","sediment__erosion_flux",grid.cell_area_at_node*dt*c.RHO),("deposition","sediment__deposition_flux",grid.cell_area_at_node*dt*c.RHO),("influx","sediment__influx",dt*c.RHO),("outflux","sediment__outflux",dt*c.RHO)]:
                req(np.allclose(node[field]*factor,fs[name],rtol=1e-13,atol=1e-7),"native flux mismatch "+name)
            reference_E=node["sediment_property__concentration"].copy()*after/c.RHO
            ds,metrics=advance(ds,fs,step)
            actual=aggregate(ds);diff=np.abs(reference_E-actual)
            value=float(diff.max());idx=int(np.argmax(diff));max_diff=max(max_diff,value)
            if (not np.isfinite(diff).all() or value>EPS) and first_fail is None:first_fail={"step":step,"node":idx,"difference_J":value,"native_J":float(reference_E[idx]),"ours_J":float(actual[idx]),"erosion_kg":float(fs["erosion"][idx]),"deposition_kg":float(fs["deposition"][idx]),"influx_kg":float(fs["influx"][idx]),"outflux_kg":float(fs["outflux"][idx])}
            history.append({"step":step,"max_tracker_difference_J":value,**metrics})
            arrays[f"s{step}_native_J"]=reference_E.copy();arrays[f"s{step}_ours_J"]=actual;arrays[f"s{step}_ROI_J"]=ds.E_roi.values.copy()
            candidate=owner.copy(deep=True);candidate.soil_mass_kg.values[:]=after.reshape(c.SHAPE)
            candidate.attrs.update(exported_soil_kg=float(owner.attrs["exported_soil_kg"])+float(suspended+bedload)*c.RHO,landlab_step=step,landlab_t_bench=step*c.LANDLAB_DT)
            c.validate(candidate,raw,seed,k,max_step=18);owner=candidate;m=after
        p=OUT/"native_soil.nc";owner.to_netcdf(p,engine="scipy")
        with xr.open_dataset(p,engine="scipy") as v:owner=v.load()
        p.unlink()
    p=OUT/f"{key}_{mode}.nc";ds.to_netcdf(p,engine="scipy")
    with xr.open_dataset(p,engine="scipy") as v:
        back=v.load();validate(back,*contract(ds,m));req(digest(back)==digest(ds),"native owner roundtrip")
    np.savez_compressed(OUT/f"{key}_{mode}_comparison.npz",**arrays)
    return {"max_tracker_difference_J":max_diff,"first_failure":first_fail,"history":history,"final_energy_J":total(ds),"export_J":float(ds.attrs["E_export_J"]),"final_hash":digest(ds)}

def main():
    global OUT
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    report={"experiment":"P1-B.2b","analytical":{},"cases":{},"failed_preregistered_gates":[],"runtime":{"python":platform.python_version(),"versions":{p:version(p) for p in ("landlab","numpy","scipy","xarray")}}}
    fail=report["failed_preregistered_gates"]
    try:
        report["analytical"]=analytical()
        archive=ROOT/"evidence/P1-B/B2a/14969549342.zip"
        req(hashlib.sha256(archive.read_bytes()).hexdigest()=="60028665102c3eb233b3c9d7ad15048be157e39d20189b69f831448bc76e635c","B2a digest")
        refs=json.loads((archive.parent/"result.json").read_text())
        req(not refs["failed_preregistered_gates"],"B2a gates")
        rawzip=ROOT/"evidence/archives/45389348480/14917807888.zip"
        req(hashlib.sha256(rawzip.read_bytes()).hexdigest()==a.ARCHIVES["raw"][1],"raw digest")
        with zipfile.ZipFile(rawzip) as f:
            rawp=OUT/"raw.npz";rawp.write_bytes(f.read("worldengine_two_seeds.npz"))
        raws=h33.read_raw(rawp);rawp.unlink()
        with zipfile.ZipFile(archive) as f:
            req(f.testzip() is None,"archive CRC")
            for mode in ("uniform","heterogeneous"):
                if fail:break
                for seed in (7321,8128):
                    for k in (0.,.1):
                        key=f"{seed}_{k}";name=f"{key}_fluxes.npz";data=f.read(name)
                        req(hashlib.sha256(data).hexdigest()==refs["outputs"][name],"flux digest")
                        with np.load(io.BytesIO(data),allow_pickle=False) as z:rec=native_case(raws[seed],seed,k,z,mode)
                        report["cases"][key+"_"+mode]=rec
                        if rec["first_failure"] is not None:fail.append(key+"_"+mode+":tracker_difference")
    except Exception as exc:
        traceback.print_exc();fail.append("instrument_or_analytical:"+repr(exc))
    report["elapsed_s"]=time.monotonic()-start
    report["outputs"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file()}
    (OUT/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_EVAL_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)

if __name__=="__main__":main()
