"""P1-B.2c: static radiation + archived sediment transport; no mechanical replay."""
import hashlib, io, json, pathlib, platform, sys, time, traceback, zipfile
from importlib.metadata import version
import numpy as np
import xarray as xr
from scipy.ndimage import label
from fipy import Grid2D, CellVariable, TransientTerm, DiffusionTerm
from fipy.solvers import LinearLUSolver
import p1a_thermal_response as a
import p1ar_exact_sphere_rays as rays
import p0d_granular_world as pd
import probe_h40_live_body_soil_epochs as h40
import probe_h37_sequential_bullet_contact as h37
import probe_e5_h31_owner_contract as c
import p0b_profile_restart as profile
import p1b2b_transport_kernel as kernel

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"p1b2c-output"
ARCHIVES={
"raw":("evidence/archives/45389348480/14917807888.zip","f2b204f42f7d055600587230093ad36b1a89f0f85ae758bae114a4b6b94a5170"),
"continuous":("evidence/archives/45447643891/14946621830.zip","34849933ca9b1c4fd1b2b1b5c9538f4093caef62950a71859f8583a58c7d9c44"),
"results":("evidence/archives/45447643891/14948035431.zip","0cb082e656d34aef4d2b5387d1e4d2cba81b0dd9337f100dde4a715956db5fad"),
"thermal":("evidence/archives/45488606315/14965401716.zip","7f919ae46e68c8da7d775f297f00a662c98f37752f7b78658b419eddcf79d512"),
"flux":("evidence/P1-B/B2a/14969549342.zip","60028665102c3eb233b3c9d7ad15048be157e39d20189b69f831448bc76e635c"),
"transport":("evidence/P1-B/B2b/14972836210.zip","c6fe9dd45be37cef2d420aa90825b53756105cb2572178959ec40a35c5000d6e")}
FIELDS={"E_roi","E_coarse_rest","E_grain"}
ATTRS={"p1b_schema","p1b_roi_node","p1b_roi_fraction","p1b_thermal_step","p1b_thermal_seconds","p1b_E_export_J","p1b_E_incident_J","p1b_flux_W_m2","p1b_k_W_mK","p1b_cp_J_kgK"}
REQ=kernel.req
FRACTION=9/10000
CP=800.
DT=10.
FLUX=600.
STEPS=6
LIMIT=.01

def digest_file(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def checked_archive(name):
    path,sha=ARCHIVES[name];p=ROOT/path
    REQ(digest_file(p)==sha,name+" archive digest")
    z=zipfile.ZipFile(p);REQ(z.testzip() is None,name+" archive CRC");return z
def read_nc(z,name):
    with xr.open_dataset(io.BytesIO(z.read(name)),engine="scipy") as ds:return ds.load()
def find(z,tail):
    names=[n for n in z.namelist() if n==tail or n.endswith("/"+tail)]
    REQ(len(names)==1,tail+" ambiguous/missing")
    return names[0]
def reference(z,report,raw,seed,k,epoch):
    key=f"{seed}_{k}";ds=read_nc(z,find(z,key+f"/epoch_{epoch}.nc"))
    soil=pd.extract_soil(ds);pd.validate(ds,soil,raw,seed,k,epoch)
    ref=report["cases"][key]["epochs"][epoch]
    REQ(pd.hash_owner(ds,soil)==ref["owner_sha"] and c.hash_owner(soil)==ref["soil_sha"],"P0-D reference")
    return ds
def roi(owner,seed):
    y,x=h40.PATCH[seed]
    REQ(38.5>=-50 and 41.5<=50 and 28.5>=-50 and 31.5<=50,"ROI control geometry")
    REQ(owner.is_core_node.values[y,x]!=0 and owner.soil_mass_kg.values[y,x]>0,"no soil receiver")
    return y*c.SHAPE[1]+x
def reduced(owner,seed,n):
    m=owner.soil_mass_kg.values.ravel().copy()
    s=kernel.make(m,roi(owner,seed),FRACTION,(n,n),0.,ids=owner.source_id.values.copy(),gm=owner.grain_mass_kg.values.copy())
    return s
def cellxy(n):
    dx=3./n
    x=38.5+(np.arange(n)+.5)*dx;y=28.5+(np.arange(n)+.5)*dx
    return np.stack(np.meshgrid(x,y),axis=-1).reshape(-1,2)
def hitmap(owner,seed,n):
    xy=cellxy(n);patch=h40.terrain(pd.extract_soil(owner),seed)
    previous=a.XY
    try:
        a.XY=xy
        return rays.checked(patch,owner.grain_state.values[:,:3].copy())
    finally:a.XY=previous
def irradiate(owner,seed,n,hits,flux=FLUX,conductivity=1.):
    ds=reduced(owner,seed,n)
    m=ds.soil_mass_kg.values[roi(owner,seed)]
    thick=m/(c.RHO*c.AREA)
    cap=c.RHO*thick*CP
    area=9./(n*n);dx=3./n
    REQ(np.isfinite(cap) and cap>0 and hits.shape==(n*n,),"capacity/shape")
    REQ(np.all((hits>=-1)&(hits<pd.N)),"receiver IDs")
    grid=Grid2D(nx=n,ny=n,dx=dx,dy=dx)
    REQ(np.allclose(np.asarray(grid.cellCenters).T+[38.5,28.5],cellxy(n),rtol=0,atol=1e-12),"ray/grid alignment")
    temp=CellVariable(mesh=grid,value=0.)
    load=CellVariable(mesh=grid,value=flux*(hits<0))
    eq=TransientTerm(coeff=cap)==DiffusionTerm(coeff=conductivity*thick)+load
    solver=LinearLUSolver(tolerance=1e-10)
    counts=np.bincount(hits[hits>=0],minlength=pd.N).astype(float)
    residuals=[]
    for j in range(1,STEPS+1):
        eq.solve(var=temp,dt=DT,solver=solver)
        ds.E_roi.values[:]=(np.asarray(temp.value)*cap*area).reshape(n,n)
        ds.E_grain.values[:]=counts*flux*area*DT*j
        ds.attrs["E_incident_J"]=flux*9.*DT*j
        kernel.validate(ds,*kernel.contract(ds,ds.soil_mass_kg.values))
        res=kernel.total(ds)-ds.attrs["E_incident_J"];REQ(abs(res)<=LIMIT,"radiative balance")
        residuals.append(res)
    return ds,{"max_balance_residual_J":max(map(abs,residuals)),"soil_J":float(ds.E_roi.sum()),"grains_J":float(ds.E_grain.sum()),"incident_J":float(ds.attrs["E_incident_J"]),"thickness_m_nominal":float(thick),"cap_J_m2K":float(cap)}
def deltaT(s):
    _,m=kernel.capacities(s)
    REQ(np.all(m>0),"ROI depleted")
    return s.E_roi.values/(m*CP)
def cold(s,b):
    delta=deltaT(b)-deltaT(s);mask=delta>=.025
    cc,n=label(mask);sizes=np.bincount(cc.ravel())[1:];area=9./mask.size
    return {"cold_cells":int(mask.sum()),"largest":int(sizes.max()) if n else 0,"cold_area_m2":float(mask.sum()*area),"largest_area_m2":float((sizes.max() if n else 0)*area),"max_deficit_K":float(delta.max()),"mean_deficit_K":float(delta.mean())}
def transport(ds,flows,owner1,k):
    initial=kernel.digest(ds);cur=ds
    metrics=[]
    for j in (1,2,3):
        fs={name:flows[f"s{j}_{name}"] for name in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
        cur,report=kernel.advance(cur,fs,j);metrics.append(report)
    REQ(kernel.digest(ds)==initial,"input transport mutation")
    REQ(np.array_equal(cur.soil_mass_kg.values,owner1.soil_mass_kg.values.ravel()),"archived terminal soil")
    REQ(abs(float(cur.attrs["mass_export_kg"])-float(owner1.attrs["exported_soil_kg"]))<=.09216,"archived mass export")
    if k==0:
        REQ(np.array_equal(cur.E_roi,ds.E_roi) and np.array_equal(cur.E_rest,ds.E_rest) and np.array_equal(cur.E_grain,ds.E_grain) and cur.attrs["E_export_J"]==0,"K0 energy/map unchanged")
    return cur,metrics
def extend(owner,ds,seed,flux,conductivity):
    n=len(ds.ty);result=owner.copy(deep=True).assign_coords(ty=np.arange(n,dtype=np.int32),tx=np.arange(n,dtype=np.int32))
    REQ(np.array_equal(ds.soil_mass_kg.values,owner.soil_mass_kg.values.ravel()),"material mismatch")
    REQ(set(ds.grain.values)==set(owner.source_id.values),"source IDs changed")
    lookup={int(s):j for j,s in enumerate(ds.grain.values)}
    order=[lookup[int(s)] for s in owner.source_id.values]
    REQ(np.array_equal(ds.grain_mass_kg.values[order],owner.grain_mass_kg.values),"grain mass changed")
    result["E_roi"]=(("ty","tx"),ds.E_roi.values.copy(),{"units":"J"})
    result["E_coarse_rest"]=(("y","x"),ds.E_rest.values.reshape(c.SHAPE).copy(),{"units":"J"})
    result["E_grain"]=(("grain",),ds.E_grain.values[order].copy(),{"units":"J"})
    result.attrs.update(p1b_schema="P1B-full-snapshot-v1",p1b_roi_node=roi(owner,seed),p1b_roi_fraction=FRACTION,p1b_thermal_step=STEPS,p1b_thermal_seconds=60.,p1b_E_export_J=float(ds.attrs["E_export_J"]),p1b_E_incident_J=float(ds.attrs["E_incident_J"]),p1b_flux_W_m2=float(flux),p1b_k_W_mK=float(conductivity),p1b_cp_J_kgK=CP)
    return result
def validate_full(ds,original,raw,seed,k,epoch,flux,conductivity):
    REQ(set(ds.data_vars)==set(original.data_vars)|FIELDS and set(ds.attrs)==set(original.attrs)|ATTRS and set(ds.coords)==set(original.coords)|{"ty","tx"},"full owner schema")
    projected=ds.drop_vars(list(FIELDS)).drop_dims(["ty","tx"])
    projected.attrs={key:v for key,v in ds.attrs.items() if key not in ATTRS}
    soil=pd.extract_soil(projected);pd.validate(projected,soil,raw,seed,k,epoch)
    REQ(pd.hash_owner(projected,soil)==pd.hash_owner(original,pd.extract_soil(original)),"original owner changed")
    for key in original.variables:
        REQ(np.array_equal(projected[key].values,original[key].values) and dict(projected[key].attrs)==dict(original[key].attrs),"original field "+key)
    for key in original.attrs:REQ(projected.attrs[key]==original.attrs[key],"original metadata "+key)
    n=len(ds.ty)
    REQ(n in (30,60,120) and np.array_equal(ds.ty,np.arange(n)) and np.array_equal(ds.tx,np.arange(n)),"grid coordinates")
    REQ(ds.E_roi.dims==("ty","tx") and ds.E_coarse_rest.dims==("y","x") and ds.E_grain.dims==("grain",),"energy dims")
    expected={"p1b_schema":"P1B-full-snapshot-v1","p1b_roi_node":roi(original,seed),"p1b_roi_fraction":FRACTION,"p1b_thermal_step":STEPS,"p1b_thermal_seconds":60.,"p1b_E_incident_J":flux*9.*60.,"p1b_flux_W_m2":flux,"p1b_k_W_mK":conductivity,"p1b_cp_J_kgK":CP}
    REQ(all(ds.attrs[q]==v for q,v in expected.items()),"thermal metadata")
    for q in FIELDS:REQ(ds[q].attrs=={"units":"J"} and np.isfinite(ds[q]).all() and np.all(ds[q].values>=-1e-7),"energy "+q)
    rr=reduced(original,seed,n)
    rr.E_roi.values[:]=ds.E_roi.values;rr.E_rest.values[:]=ds.E_coarse_rest.values.ravel();rr.E_grain.values[:]=ds.E_grain.values
    rr.attrs.update(step=3*epoch,E_incident_J=float(ds.attrs["p1b_E_incident_J"]),E_export_J=float(ds.attrs["p1b_E_export_J"]))
    REQ(rr.attrs["E_export_J"]>=0,"negative energy export")
    kernel.validate(rr,*kernel.contract(rr,rr.soil_mass_kg.values))
    h=hashlib.sha256(pd.hash_owner(projected,soil).encode())
    for q in sorted(FIELDS):h.update(q.encode()+ds[q].values.tobytes())
    h.update(json.dumps({q:c.normal(ds.attrs[q]) for q in sorted(ATTRS)},sort_keys=True).encode())
    return h.hexdigest()
def persist(owner0,owner1,s,t,seed,k,raw,flux,conductivity,folder,name):
    hashes=[]
    for epoch,base,red in [(0,owner0,s),(1,owner1,t)]:
        ds=extend(base,red,seed,flux,conductivity)
        sha=validate_full(ds,base,raw,seed,k,epoch,flux,conductivity)
        path=folder/f"{name}_epoch{epoch}.nc";ds.to_netcdf(path,engine="scipy")
        with xr.open_dataset(path,engine="scipy") as f:back=f.load()
        REQ(validate_full(back,base,raw,seed,k,epoch,flux,conductivity)==sha,"full NetCDF/hash")
        if epoch==0 and name=="actual":
            bad=back.copy(deep=True);bad.E_roi.values[0,0]+=1.
            try:validate_full(bad,base,raw,seed,k,epoch,flux,conductivity)
            except ValueError:pass
            else:raise ValueError("corrupted full owner accepted")
        hashes.append(sha)
    return hashes
def evaluate(owner0,owner1,seed,k,raw,flows,folder,n,full=True):
    original0=pd.hash_owner(owner0,pd.extract_soil(owner0));original1=pd.hash_owner(owner1,pd.extract_soil(owner1))
    hits=hitmap(owner0,seed,n);barehits=np.full_like(hits,-1)
    REQ(np.any(hits>=0),"no intercepted rays")
    s,sm=irradiate(owner0,seed,n,hits);b,bm=irradiate(owner0,seed,n,barehits)
    REQ(np.max(np.abs(deltaT(b)-.25))<=1e-7,"bare initial analytical")
    REQ(abs(sm["incident_J"]-324000.)<=1e-7 and abs(bm["incident_J"]-324000.)<=1e-7,"incident dose")
    t,tm=transport(s,flows,owner1,k);u,um=transport(b,flows,owner1,k)
    pre=cold(s,b);post=cold(t,u)
    info={"resolution":n,"intercepted_cells":int((hits>=0).sum()),"before":pre,"after":post,"retention_largest_area_ratio":post["largest_area_m2"]/pre["largest_area_m2"] if pre["largest_area_m2"] else None,"radiation":sm,"bare_radiation":bm,"transport":tm,"bare_transport":um,"material_energy_export_J":float(t.attrs["E_export_J"])}
    others=kernel.aggregate(t);others[roi(owner1,seed)]=0
    info["neighbor_soil_energy_J"]=float(others.sum())
    if k:REQ(info["neighbor_soil_energy_J"]>0 or info["material_energy_export_J"]>0,"no energy transported out of ROI node")
    arrays={"hits":hits,"actual_before_J":s.E_roi.values,"bare_before_J":b.E_roi.values,"actual_after_J":t.E_roi.values,"bare_after_J":u.E_roi.values}
    if full:
        REQ(pre["largest"]>=10,"cold connected component")
        info["owner_hashes_actual"]=persist(owner0,owner1,s,t,seed,k,raw,FLUX,1.,folder,"actual")
        info["owner_hashes_bare"]=persist(owner0,owner1,b,u,seed,k,raw,FLUX,1.,folder,"bare")
        zero,zm=irradiate(owner0,seed,n,hits,0.)
        REQ(kernel.total(zero)==0.,"zero flux")
        zt,_=transport(zero,flows,owner1,k)
        persist(owner0,owner1,zero,zt,seed,k,raw,0.,1.,folder,"zero_flux")
        nk,km=irradiate(owner0,seed,n,hits,FLUX,0.)
        REQ(np.max(np.abs(deltaT(nk).ravel()-.25*(hits<0)))<=1e-7,"zero conductivity analytical")
        nt,_=transport(nk,flows,owner1,k);persist(owner0,owner1,nk,nt,seed,k,raw,FLUX,0.,folder,"zero_k")
        rng=np.random.default_rng(2026112400+seed);null_pre=[];null_post=[]
        for i in range(32):
            shuffled=rng.permutation(hits)
            REQ(np.array_equal(np.sort(shuffled),np.sort(hits)),"null dose")
            ns,nm=irradiate(owner0,seed,n,shuffled)
            REQ(abs(float(ns.E_roi.sum()-s.E_roi.sum()))<=LIMIT,"null soil dose")
            nst,_=transport(ns,flows,owner1,k)
            null_pre.append(cold(ns,b)["largest"]);null_post.append(cold(nst,u)["largest"])
            arrays[f"null{i}_before_J"]=ns.E_roi.values;arrays[f"null{i}_after_J"]=nst.E_roi.values;arrays[f"null{i}_hits"]=shuffled
        p95=int(np.quantile(null_pre,.95,method="higher"));equal=sum(x>=pre["largest"] for x in null_pre)
        info["null"]={"before":null_pre,"after":null_post,"p95":p95,"equal_or_exceed":equal,"rank":(1+equal)/33}
        if k:REQ(pre["largest"]>=p95+4 and equal==0,"dose-matched coherence")
        # Conservative remap of actual transported maps, independently from native-resolution sensitivity.
        info["remap"]={}
        for n2 in (30,120):
            rt=kernel.remap(t,n2);ru=kernel.remap(u,n2)
            info["remap"][str(n2)]={**cold(rt,ru),"energy_error_J":kernel.total(rt)-kernel.total(t)}
        if seed==7321 and k:
            rs,rm=irradiate(owner0,seed,n,hitmap(owner0,seed,n));rr,rtm=transport(rs,flows,owner1,k)
            REQ(kernel.digest(rs)==kernel.digest(s) and kernel.digest(rr)==kernel.digest(t),"bit exact replay")
            info["replay_exact"]=True
        np.savez_compressed(folder/"maps_60.npz",**arrays)
    else:
        np.savez_compressed(folder/f"maps_{n}.npz",**arrays)
    REQ(pd.hash_owner(owner0,pd.extract_soil(owner0))==original0 and pd.hash_owner(owner1,pd.extract_soil(owner1))==original1,"source mutated")
    return info
def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False);report={"experiment":"P1-B.2c","cases":{},"failed_preregistered_gates":[],"runtime":{"python":platform.python_version(),"versions":{p:version(p) for p in ("numpy","scipy","xarray","fipy","pybullet")}}}
    try:
        profile.runtime(sys.executable)
        old=(ROOT/"experiments/p1b2b_energy_transport.py").read_text();new=(ROOT/"experiments/p1b2b_transport_kernel.py").read_text()
        REQ(new[new.index("EPS=.01"):]==old[old.index("EPS=.01"):old.index("def analytical():")],"tested kernel changed")
        for name in ("results","thermal","transport"):
            with checked_archive(name):pass
        with checked_archive("raw") as z:
            p=OUT/"raw.npz";p.write_bytes(z.read(find(z,"worldengine_two_seeds.npz")))
        raws=h37.read_raw(p);p.unlink()
        refs=json.loads((ROOT/"evidence/P1-B/B2a/result.json").read_text())
        REQ(not refs["failed_preregistered_gates"],"B2a gate")
        b2b=json.loads((ROOT/"evidence/P1-B/B2b/result.json").read_text())
        REQ(not b2b["failed_preregistered_gates"],"B2b gate")
        with checked_archive("continuous") as z,checked_archive("flux") as f:
            p0=json.loads(z.read(find(z,"result.json")))
            for seed in (7321,8128):
                for k in (0.,.1):
                    key=f"{seed}_{k}";folder=OUT/key;folder.mkdir()
                    try:
                        o0=reference(z,p0,raws[seed],seed,k,0);o1=reference(z,p0,raws[seed],seed,k,1)
                        name=f"{key}_fluxes.npz";data=f.read(name)
                        REQ(hashlib.sha256(data).hexdigest()==refs["outputs"][name],"flux digest")
                        with np.load(io.BytesIO(data),allow_pickle=False) as flow:
                            main=evaluate(o0,o1,seed,k,raws[seed],flow,folder,60)
                            sensitivity={}
                            for n in (30,120):
                                other=evaluate(o0,o1,seed,k,raws[seed],flow,folder,n,False)
                                other["relative_area_change"]=(other["before"]["cold_area_m2"]-main["before"]["cold_area_m2"])/main["before"]["cold_area_m2"]
                                other["relative_soil_energy_change"]=(other["radiation"]["soil_J"]-main["radiation"]["soil_J"])/main["radiation"]["soil_J"]
                                sensitivity[str(n)]=other
                            main["sensitivity"]=sensitivity;report["cases"][key]=main
                    except Exception as exc:
                        traceback.print_exc();report["failed_preregistered_gates"].append(key+":"+repr(exc))
    except Exception as exc:
        traceback.print_exc();report["failed_preregistered_gates"].append("instrument:"+repr(exc))
    report["elapsed_s"]=time.monotonic()-start
    report["outputs"]={str(p.relative_to(OUT)):digest_file(p) for p in OUT.rglob("*") if p.is_file()}
    (OUT/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_EVAL_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if report["failed_preregistered_gates"]:raise SystemExit(1)
if __name__=="__main__":main()
