"""P1-C.2a continuous producer. Existing P0-D physics, accumulated passive energy."""
import hashlib,io,json,pathlib,platform,resource,shutil,subprocess,sys,tempfile,time,traceback,zipfile
from importlib.metadata import version
import numpy as np
import pybullet as p
import xarray as xr
from fipy import Grid2D,CellVariable,TransientTerm,DiffusionTerm
from fipy.solvers import LinearLUSolver
import p0d_granular_world as pd
import probe_e5_h31_owner_contract as c
import probe_e6_h35_granular_midfeed_restart as h35
import probe_e3_h26_gravity_feed as h26
import probe_h40_live_body_soil_epochs as h40
import probe_h41_aggregate_transit as h41
import probe_h37_sequential_bullet_contact as h37
import p0b_joint_restart as b
import p0b_profile_restart as profile
import p1b2c_thermal_snapshots as thermal
import p1b2b_transport_kernel as ker
import p1a_thermal_response as a
import p1ar_exact_sphere_rays as rays

ROOT=pathlib.Path(__file__).resolve().parents[1]
WORKER=pathlib.Path(__file__).with_name("p1c_live_worker.py")
FIELDS=thermal.FIELDS
ATTRS={"p1c_schema","p1c_mode","p1c_roi_node","p1c_roi_fraction","p1c_thermal_step","p1c_thermal_seconds","p1c_E_export_J","p1c_E_incident_J","p1c_flux_W_m2","p1c_k_W_mK","p1c_cp_J_kgK"}
B2C_PATH=ROOT/"evidence/P1-B/B2c/14973541094.zip"
B2C_SHA="4a763034d7b43416e37db5d4bc56e382cb67a8534425bc6c4f1884cb680a1789"
REQ=ker.req
CP=800.;N=60;AREA=9/(N*N);FLUX=600.;DT=10.;STEPS=6;LIMIT=.01

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def normal(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    return x
def digest(ds):
    h=hashlib.sha256()
    for key in sorted(ds.variables):
        v=ds[key];arr=np.asarray(v.values)
        typ="text" if arr.dtype.kind in "OUS" else arr.dtype.newbyteorder("<").str
        h.update(json.dumps([key,list(v.dims),list(arr.shape),typ,{q:normal(z) for q,z in v.attrs.items()}],sort_keys=True).encode())
        if typ=="text":h.update(json.dumps(arr.tolist(),ensure_ascii=False).encode())
        else:h.update(arr.astype(arr.dtype.newbyteorder("<"),copy=False).tobytes())
    h.update(json.dumps({q:normal(v) for q,v in ds.attrs.items()},sort_keys=True).encode())
    return h.hexdigest()
def flux_digest(arrays):
    h=hashlib.sha256()
    for key in sorted(arrays):
        arr=np.asarray(arrays[key]);h.update(json.dumps([key,arr.dtype.str,list(arr.shape)]).encode());h.update(arr.tobytes())
    return h.hexdigest()
def project(ds):
    out=ds.drop_vars(list(FIELDS)).drop_dims(["ty","tx"])
    out.attrs={k:v for k,v in ds.attrs.items() if k not in ATTRS}
    return out
def intensity(mode):return 0. if mode=="zero_flux" else FLUX
def unpack(owner,seed,epoch,mode):
    base=project(owner);soil=pd.extract_soil(base)
    red=thermal.reduced(base,seed,N)
    red.E_roi.values[:]=owner.E_roi.values
    red.E_rest.values[:]=owner.E_coarse_rest.values.ravel()
    red.E_grain.values[:]=owner.E_grain.values
    red.attrs.update(step=epoch*3,E_export_J=float(owner.attrs["p1c_E_export_J"]),E_incident_J=float(owner.attrs["p1c_E_incident_J"]),mass_export_kg=float(soil.attrs["exported_soil_kg"]),initial_mass_kg=c.MASS0+pd.N*pd.M)
    ker.validate(red,*ker.contract(red,soil.soil_mass_kg.values.ravel()))
    return red
def publish(base,red,seed,epoch,mode):
    ids=base.source_id.values
    REQ(set(ids)==set(red.grain.values),"source IDs")
    order=[int(np.flatnonzero(red.grain.values==sid)[0]) for sid in ids]
    REQ(np.array_equal(red.soil_mass_kg.values,base.soil_mass_kg.values.ravel()),"soil stock")
    REQ(np.array_equal(red.grain_mass_kg.values[order],base.grain_mass_kg.values),"grain stock")
    out=base.copy(deep=True).assign_coords(ty=np.arange(N,dtype=np.int32),tx=np.arange(N,dtype=np.int32))
    out["E_roi"]=(("ty","tx"),red.E_roi.values.copy(),{"units":"J"})
    out["E_coarse_rest"]=(("y","x"),red.E_rest.values.reshape(c.SHAPE).copy(),{"units":"J"})
    out["E_grain"]=(("grain",),red.E_grain.values[order].copy(),{"units":"J"})
    out.attrs.update(p1c_schema="P1C-v1",p1c_mode=mode,p1c_roi_node=thermal.roi(base,seed),p1c_roi_fraction=thermal.FRACTION,p1c_thermal_step=6*(epoch+1),p1c_thermal_seconds=60.*(epoch+1),p1c_E_export_J=float(red.attrs["E_export_J"]),p1c_E_incident_J=float(red.attrs["E_incident_J"]),p1c_flux_W_m2=intensity(mode),p1c_k_W_mK=1.,p1c_cp_J_kgK=CP)
    return out
def validate(owner,raw,seed,k,epoch,mode):
    base=project(owner);soil=pd.extract_soil(base)
    pd.validate(base,soil,raw,seed,k,epoch)
    REQ(set(owner.data_vars)==set(base.data_vars)|FIELDS and set(owner.attrs)==set(base.attrs)|ATTRS and set(owner.coords)==set(base.coords)|{"ty","tx"},"owner schema")
    expected={"p1c_schema":"P1C-v1","p1c_mode":mode,"p1c_roi_node":thermal.roi(base,seed),"p1c_roi_fraction":thermal.FRACTION,"p1c_thermal_step":6*(epoch+1),"p1c_thermal_seconds":60.*(epoch+1),"p1c_E_incident_J":intensity(mode)*9.*60.*(epoch+1),"p1c_flux_W_m2":intensity(mode),"p1c_k_W_mK":1.,"p1c_cp_J_kgK":CP}
    REQ(all(owner.attrs[q]==v for q,v in expected.items()),"energy configuration/clocks")
    REQ(np.array_equal(owner.ty,np.arange(N)) and np.array_equal(owner.tx,np.arange(N)),"ROI coordinates")
    for key,dims in [("E_roi",("ty","tx")),("E_coarse_rest",("y","x")),("E_grain",("grain",))]:
        REQ(owner[key].dims==dims and dict(owner[key].attrs)=={"units":"J"},"energy field "+key)
    red=unpack(owner,seed,epoch,mode)
    REQ(abs(ker.total(red)-expected["p1c_E_incident_J"])<=LIMIT,"full energy balance")
    if mode=="zero_flux":REQ(ker.total(red)==0.,"zero irradiation")
    return digest(owner)
def irradiate(red,base,seed,mode,frozen=None):
    points=base.grain_state.values[:,:3].copy() if frozen is None else frozen
    previous=a.XY
    try:
        a.XY=thermal.cellxy(N)
        hits=rays.checked(h40.terrain(pd.extract_soil(base),seed),points)
    finally:a.XY=previous
    REQ(np.all((hits>=-1)&(hits<pd.N)),"ray identity")
    # Map positional ray winner to permanent source ID, then to red view index.
    lookup={int(s):i for i,s in enumerate(red.grain.values)}
    mapped=np.full(hits.shape,-1,dtype=np.int32)
    idx=np.flatnonzero(hits>=0)
    mapped[idx]=[lookup[int(base.source_id.values[j])] for j in hits[idx]]
    m=float(red.soil_mass_kg.values[thermal.roi(base,seed)])
    thick=m/(c.RHO*c.AREA);cap=c.RHO*thick*CP
    REQ(cap>0 and np.isfinite(cap),"ROI mass")
    grid=Grid2D(nx=N,ny=N,dx=.05,dy=.05)
    REQ(np.allclose(np.asarray(grid.cellCenters).T+[38.5,28.5],thermal.cellxy(N),rtol=0,atol=1e-12),"ray grid")
    temp=CellVariable(mesh=grid,value=red.E_roi.values.ravel()/(cap*AREA))
    flux=intensity(mode)
    eq=TransientTerm(coeff=cap)==DiffusionTerm(coeff=thick)+CellVariable(mesh=grid,value=flux*(hits<0))
    solver=LinearLUSolver(tolerance=1e-10)
    before=red.copy(deep=True);initial_inc=float(red.attrs["E_incident_J"])
    counts=np.bincount(mapped[mapped>=0],minlength=pd.N).astype(float)
    residual=[]
    for j in range(1,STEPS+1):
        eq.solve(var=temp,dt=DT,solver=solver)
        red.E_roi.values[:]=(np.asarray(temp.value)*cap*AREA).reshape(N,N)
        red.E_grain.values[:]=before.E_grain.values+counts*flux*AREA*DT*j
        red.attrs["E_incident_J"]=initial_inc+flux*9.*DT*j
        ker.validate(red,*ker.contract(red,red.soil_mass_kg.values))
        err=ker.total(red)-red.attrs["E_incident_J"];REQ(abs(err)<=LIMIT,"thermal balance");residual.append(err)
    return red,mapped,{"residual_max_J":max(map(abs,residual)),"grain_added_J":float((red.E_grain-before.E_grain).sum()),"intercepted_cells":int((hits>=0).sum()),"max_grain_K":float((300+red.E_grain.values/(CP*red.grain_mass_kg.values)).max()),"max_roi_K":float(300+thermal.deltaT(red).max()),"roi_J":float(red.E_roi.sum()),"rest_J":float(red.E_rest.sum()),"grain_J":float(red.E_grain.sum()),"export_J":float(red.attrs["E_export_J"]),"incident_J":float(red.attrs["E_incident_J"])}
def request(python,rawpath,raw,seed,k,folder,epoch,old):
    target=3*epoch;dest=folder/f"block_{target}"
    REQ(not dest.exists(),"duplicate soil request")
    prev="-" if not epoch else str(folder/f"block_{target-3}"/"soil.nc")
    REQ((old is None)==(epoch==0),"predecessor")
    proc=subprocess.run([python,str(WORKER),str(rawpath),str(seed),str(k),prev,str(target),str(dest)],capture_output=True,text=True)
    if proc.returncode:raise RuntimeError(proc.stdout+"\n"+proc.stderr)
    event=json.loads((dest/"event.json").read_text())
    line=next(l for l in proc.stdout.splitlines() if l.startswith("P1C_SOIL_JSON "))
    REQ(event==json.loads(line.split(" ",1)[1]),"worker event")
    for key in ("soil.nc","flux.npz"):REQ(sha(dest/key)==event["files"][key],"worker digest")
    soil=b.read_owner(dest/"soil.nc");c.validate(soil,raw,seed,k,max_step=18)
    REQ(event["seed"]==seed and event["K"]==k and event["to"]==target and event["from"]==(None if not epoch else target-2) and event["input_sha"]==(None if old is None else c.hash_owner(old)) and event["output_sha"]==c.hash_owner(soil),"worker lineage")
    with np.load(dest/"flux.npz",allow_pickle=False) as z:arrays={key:z[key].copy() for key in z.files}
    return soil,event,arrays
def compare_flux(arrays,ref,epoch):
    if not epoch:REQ(not arrays,"epoch0 flux");return {}
    result={}
    for key,v in arrays.items():
        expect=ref[key]
        if v.dtype.kind in "biu":REQ(np.array_equal(v,expect),"flow topology "+key)
        else:
            REQ(np.allclose(v,expect,rtol=1e-13,atol=1e-7),"flow numeric "+key)
            result[key]=float(np.max(np.abs(v-expect)))
    names=("before","after","erosion","deposition","influx","outflux","bedrock","receiver","order","active","core")
    expected={f"s{j}_{name}" for j in range(3*epoch-2,3*epoch+1) for name in names}
    REQ(set(arrays)==expected,"three native steps: exact field set")
    return {"max_abs":max(result.values(),default=0.),"fields":len(arrays)}
def save(folder,owner,state,recipe,history,raw,seed,k,mode):
    stage=folder.with_name(folder.name+".staging");stage.mkdir()
    try:
        owner.to_netcdf(stage/"owner.nc",engine="scipy")
        p.saveBullet(str(stage/"native.bullet"),physicsClientId=state["cid"])
        key=validate(owner,raw,seed,k,pd.MID,mode)
        back=b.read_owner(stage/"owner.nc");REQ(validate(back,raw,seed,k,pd.MID,mode)==key,"checkpoint roundtrip")
        meta={"profile":profile.PROFILE,"seed":seed,"K":k,"epoch":pd.MID,"mode":mode,"owner_sha":key,"recipe":recipe,"history":history,"floor":state["floor"],"completed_soil_step":9,"files":{f.name:sha(f) for f in stage.iterdir()}}
        b.dump(stage/"manifest.json",meta);stage.replace(folder)
    finally:
        if stage.exists():shutil.rmtree(stage)
def run(rawpath,python,raw,seed,k,folder,out,reference,flux_ref,epoch0_ref,mode="actual",observer=True):
    start=time.monotonic();state=None;old=None;owner=None;recipe=[];history=[];epochs=[];requests=[];frozen=None
    out.mkdir(parents=True,exist_ok=False)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            private=pathlib.Path(tmp)
            for epoch in range(7):
                soil,event,arrays=request(python,rawpath,raw,seed,k,private,epoch,old)
                try:request(python,rawpath,raw,seed,k,private,epoch,old)
                except ValueError:pass
                else:raise ValueError("duplicate accepted")
                # References are used for post-operation comparison only.
                ref=reference["epochs"][epoch];REQ(c.hash_owner(soil)==ref["soil_sha"],"soil reference")
                flux_comparison=compare_flux(arrays,flux_ref,epoch)
                if epoch:
                    red=unpack(owner,seed,epoch-1,mode);tm=[]
                    for j in range(3*epoch-2,3*epoch+1):
                        fs={key:arrays[f"s{j}_{key}"] for key in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
                        red,m=ker.advance(red,fs,j);tm.append(m)
                    REQ(np.array_equal(red.soil_mass_kg.values,soil.soil_mass_kg.values.ravel()),"worker final material")
                    REQ(abs(float(red.attrs["mass_export_kg"])-float(soil.attrs["exported_soil_kg"]))<=.09216,"worker export")
                else:red=None;tm=[]
                newpatch=h40.terrain(soil,seed);recipe.append(newpatch.tolist());trace=hashlib.sha256()
                if epoch==0:
                    state,patch=pd.setup(soil,seed)
                    plan_sha=h41.digest(np.asarray([xyz for _,xyz in state["plan"]]))
                    for step in range(1,h26.MAX_STEPS+1):
                        done=h35.step_scene(state,step);pd.trace_step(trace,state,step)
                        if done:break
                    else:raise ValueError("feed timeout")
                    form=step;REQ(not state["plan"] and state["max_ledger_error_kg"]<=1e-12,"feed ledger")
                else:
                    pd.integrity(state);pd.replace(state,newpatch);pd.integrity(state)
                    for step in range(1,pd.BLOCK+1):
                        p.stepSimulation(physicsClientId=state["cid"]);pd.trace_step(trace,state,form+(epoch-1)*pd.BLOCK+step)
                    if k:REQ(float(soil.attrs["exported_soil_kg"])>float(old.attrs["exported_soil_kg"]) and not np.array_equal(soil.soil_mass_kg.values,old.soil_mass_kg.values),"erosion stopped")
                    else:REQ(float(soil.attrs["exported_soil_kg"])==0. and np.array_equal(soil.soil_mass_kg.values,old.soil_mass_kg.values),"K0 changed")
                patch=newpatch;pd.integrity(state);cid=state["cid"]
                z=h40.height(patch,70.,70.);hit=p.rayTest([70.,70.,z+10],[70.,70.,z-10],physicsClientId=cid)[0]
                REQ(hit[0]==state["floor"] and abs(hit[3][2]-z)<.001,"terrain ray")
                base=pd.snapshot(soil,state,raw,seed,k,epoch,form,plan_sha)
                base_sha=pd.hash_owner(base,soil);trace_sha=trace.hexdigest()
                REQ(base_sha==ref["owner_sha"] and trace_sha==ref["trace_sha"],"P0-D owner/trace mismatch")
                if not epoch:
                    red=thermal.reduced(base,seed,N)
                    if mode=="frozen":frozen=base.grain_state.values[:,:3].copy()
                red,hits,rad=irradiate(red,base,seed,mode,frozen)
                owner=publish(base,red,seed,epoch,mode)
                owner_sha=validate(owner,raw,seed,k,epoch,mode)
                REQ(pd.hash_owner(project(owner),soil)==base_sha,"P0-D projection")
                if not epoch and mode=="actual":
                    for field in FIELDS:REQ(np.array_equal(owner[field].values,epoch0_ref[field].values),"B.2c epoch0 "+field)
                if observer:
                    obs=pd.observe(base,soil,state,patch)
                    pose=h41.snap(cid,state["ids"]).copy();prior=digest(owner)
                    for _ in range(2):
                        check=unpack(owner,seed,epoch,mode);_temp=thermal.deltaT(check);REQ(np.isfinite(_temp).all(),"energy query")
                    REQ(digest(owner)==prior and np.array_equal(pose,h41.snap(cid,state["ids"])),"energy query mutated")
                else:obs=None
                path=out/f"epoch_{epoch}.nc";owner.to_netcdf(path,engine="scipy")
                back=b.read_owner(path);REQ(validate(back,raw,seed,k,epoch,mode)==owner_sha,"owner roundtrip")
                np.savez_compressed(out/f"flux_{epoch}.npz",**arrays)
                np.savez_compressed(out/f"ray_{epoch}.npz",source_id=np.where(hits>=0,red.grain.values[np.maximum(hits,0)],-1))
                history.append(trace_sha);requests.append(event)
                rec={"epoch":epoch,"steps":form+epoch*pd.BLOCK,"owner_sha":owner_sha,"projection_sha":base_sha,"soil_sha":c.hash_owner(soil),"trace_sha":trace_sha,"flux_sha":flux_digest(arrays),"ray_sha":hashlib.sha256(hits.tobytes()).hexdigest(),"radiation":rad,"transport":tm,"flux_comparison":flux_comparison,"observation":obs}
                epochs.append(rec)
                if epoch==pd.MID and mode=="actual" and observer:save(folder,owner,state,recipe,history,raw,seed,k,mode)
                old=soil
            REQ([e["to"] for e in requests]==list(range(0,19,3)),"request order")
            return {"epochs":epochs,"history":history,"formation_steps":form,"new_bullet_steps":form+6*pd.BLOCK,"requests":requests,"wall_seconds":time.monotonic()-start,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"child_peak_rss_kib":resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    finally:
        if state is not None:p.disconnect(physicsClientId=state["cid"])
def main(python):
    root=ROOT/"p1c-producer";checkpoints=ROOT/"p1c-checkpoints";root.mkdir();checkpoints.mkdir()
    result={"experiment":"P1-C.2a","cases":{},"controls":{},"failed_preregistered_gates":[],"runtime":{}}
    failures=result["failed_preregistered_gates"]
    try:
        result["runtime"]={"bullet_profile":profile.runtime(sys.executable),"landlab_profile":profile.runtime(python),"versions":{n:version(n) for n in ("numpy","scipy","fipy","pybullet","xarray")},"python":platform.python_version()}
        REQ(sha(B2C_PATH)==B2C_SHA,"B.2c archive")
        with thermal.checked_archive("raw") as z:
            rawpath=root/"raw.npz";rawpath.write_bytes(z.read(thermal.find(z,"worldengine_two_seeds.npz")))
        raws=h37.read_raw(rawpath)
        with thermal.checked_archive("continuous") as z,thermal.checked_archive("flux") as f,zipfile.ZipFile(B2C_PATH) as th:
            refs=json.loads(z.read(thermal.find(z,"result.json")));REQ(not refs["failed_preregistered_gates"],"P0-D reference")
            fref=json.loads((ROOT/"evidence/P1-B/B2a/result.json").read_text())
            bref=json.loads((ROOT/"evidence/P1-B/B2c/result.json").read_text());REQ(not bref["failed_preregistered_gates"],"B.2c gates")
            for seed in (7321,8128):
                for k in (0.,.1):
                    key=f"{seed}_{k}"
                    try:
                        data=f.read(key+"_fluxes.npz");REQ(hashlib.sha256(data).hexdigest()==fref["outputs"][key+"_fluxes.npz"],"flux reference")
                        refpath=thermal.find(th,key+"/actual_epoch0.nc")
                        refdata=th.read(refpath);REQ(hashlib.sha256(refdata).hexdigest()==bref["outputs"][key+"/actual_epoch0.nc"],"epoch0 reference")
                        ref0=thermal.read_nc(th,refpath)
                        with np.load(io.BytesIO(data),allow_pickle=False) as flux_ref:
                            result["cases"][key]=run(rawpath,python,raws[seed],seed,k,checkpoints/key,root/key,refs["cases"][key],flux_ref,ref0)
                    except Exception as e:traceback.print_exc();failures.append(key+":"+repr(e))
            if not failures:
                key="7321_0.1";expected=result["cases"][key];ref0=thermal.read_nc(th,thermal.find(th,key+"/actual_epoch0.nc"))
                with np.load(io.BytesIO(f.read(key+"_fluxes.npz")),allow_pickle=False) as flux_ref:
                    for mode,obs in [("actual",False),("zero_flux",True),("frozen",True)]:
                        label="observerless" if mode=="actual" else mode
                        rr=run(rawpath,python,raws[7321],7321,.1,root/"unused",root/label,refs["cases"][key],flux_ref,ref0,mode,obs)
                        result["controls"][label]=rr
                        REQ(all((a["projection_sha"],a["trace_sha"])==(b["projection_sha"],b["trace_sha"]) for a,b in zip(rr["epochs"],expected["epochs"])),"control physical trajectories")
                        if label=="observerless":
                            REQ(all({k:v for k,v in a.items() if k!="observation"}=={k:v for k,v in b.items() if k!="observation"} for a,b in zip(rr["epochs"],expected["epochs"])),"observer altered dynamics")
                        if mode=="frozen":
                            differences=[]
                            for epoch in range(7):
                                x=b.read_owner(root/key/f"epoch_{epoch}.nc");y=b.read_owner(root/label/f"epoch_{epoch}.nc")
                                differences.append(max(float(np.abs(x[field].values-y[field].values).max()) for field in FIELDS))
                            REQ(differences[0]<=1e-7 and max(differences[1:])>1e-7,"frozen geometry contrast")
                            rr["max_field_differences_J"]=differences
        rawpath.unlink()
    except Exception as e:traceback.print_exc();failures.append("instrument:"+repr(e))
    result["outputs"]={str(q.relative_to(root)):sha(q) for q in root.rglob("*") if q.is_file()}
    result["checkpoints"]={str(q.relative_to(checkpoints)):sha(q) for q in checkpoints.rglob("*") if q.is_file()}
    (root/"result.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("WORLD_EVAL_JSON "+json.dumps(result,sort_keys=True),flush=True)
    if failures:raise SystemExit(1)
if __name__=="__main__":main(*sys.argv[1:])
