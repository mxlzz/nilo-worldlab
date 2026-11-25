"""P1-C.2b consumer: no producer futures, restore native checkpoint, advance suffix."""
import hashlib,json,pathlib,platform,resource,shutil,sys,tempfile,time,traceback,zipfile
import numpy as np
import pybullet as p
import p1c_live_world as live
from p1c_live_world import pd,c,h35,h40,h41,h37,b,profile,thermal,ker,REQ,sha

ROOT=live.ROOT
OUT=ROOT/"p1c-consumer"
CHECK=ROOT/"evidence/P1-C/C2a/14977883611.zip"
CHECK_SHA="a81ecf3bb6091b2e4a0578baaa3b84702a1e40f2ee8db93e7a91f1577e1f5be6"
RAW=ROOT/"evidence/archives/45389348480/14917807888.zip"
RAW_SHA="f2b204f42f7d055600587230093ad36b1a89f0f85ae758bae114a4b6b94a5170"
KEYS=("7321_0.0","7321_0.1","8128_0.0","8128_0.1")

def verify(folder,raw,seed,k):
    meta=json.loads((folder/"manifest.json").read_text())
    REQ(set(meta)=={"profile","seed","K","epoch","mode","owner_sha","recipe","history","floor","completed_soil_step","files"},"manifest schema")
    profile.verify_profile(meta)
    REQ((meta["seed"],meta["K"],meta["epoch"],meta["mode"],meta["completed_soil_step"])==(seed,k,3,"actual",9),"checkpoint identity")
    REQ(set(meta["files"])=={"owner.nc","native.bullet"} and {q.name for q in folder.iterdir()}=={"owner.nc","native.bullet","manifest.json"},"checkpoint files")
    for name,h in meta["files"].items():REQ(sha(folder/name)==h,"checkpoint digest "+name)
    owner=b.read_owner(folder/"owner.nc")
    REQ(live.validate(owner,raw,seed,k,3,"actual")==meta["owner_sha"],"checkpoint logical digest")
    ds=live.project(owner);soil=pd.extract_soil(ds)
    recipe=np.asarray(meta["recipe"])
    REQ(recipe.shape==(4,2,2) and np.isfinite(recipe).all() and np.array_equal(recipe[-1],h40.terrain(soil,seed)),"recipe")
    REQ(len(meta["history"])==4 and all(isinstance(x,str) and len(x)==64 for x in meta["history"]),"prefix")
    return meta,owner,ds,soil

def rejected(fn,label):
    try:fn()
    except ValueError:return True
    raise ValueError("accepted corruption "+label)

def controls(folder,owner,raw,seed,k):
    out={}
    def mutate(label,fn):
        bad=owner.copy(deep=True)
        fn(bad)
        out[label]=rejected(lambda:live.validate(bad,raw,seed,k,3,"actual"),label)
    mutate("energy_plus_1",lambda x:x.E_roi.values.__setitem__((0,0),x.E_roi.values[0,0]+1.))
    mutate("duplicate_stock",lambda x:x.__setitem__("extra_mass",x.soil_mass_kg.copy()))
    mutate("source_id",lambda x:x.source_id.values.__setitem__(0,x.source_id.values[1]))
    mutate("clock",lambda x:x.attrs.__setitem__("p1c_thermal_step",25))
    with tempfile.TemporaryDirectory() as tmp:
        dest=pathlib.Path(tmp)/"bad";shutil.copytree(folder,dest)
        q=dest/"native.bullet";arr=bytearray(q.read_bytes());arr[-1]^=1;q.write_bytes(arr)
        out["cache"]=rejected(lambda:verify(dest,raw,seed,k),"cache")
        shutil.rmtree(dest);shutil.copytree(folder,dest)
        bad=owner.copy(deep=True);e=bad.E_roi.values.ravel();i=int(np.argmax(e));j=(i+1)%len(e)
        REQ(e[i]>1.,"corruption control stock")
        e[i]-=1.;e[j]+=1.
        # Valid total but different spatial distribution: digest must reject.
        REQ(live.validate(bad,raw,seed,k,3,"actual")!=live.digest(owner),"compensated corruption hash")
        bad.to_netcdf(dest/"owner.nc",engine="scipy")
        out["compensated_file"]=rejected(lambda:verify(dest,raw,seed,k),"compensated file")
        meta=json.loads((dest/"manifest.json").read_text());meta["files"]["owner.nc"]=sha(dest/"owner.nc")
        (dest/"manifest.json").write_text(json.dumps(meta,sort_keys=True))
        out["compensated_logical"]=rejected(lambda:verify(dest,raw,seed,k),"compensated logical")
    return out

def run(rawpath,python,raw,seed,k,folder,out):
    start=time.monotonic();state=None;out.mkdir()
    before={q.name:sha(q) for q in folder.iterdir()}
    meta,owner,ds,old=verify(folder,raw,seed,k)
    tests=controls(folder,owner,raw,seed,k)
    recipe=meta["recipe"];history=list(meta["history"]);epochs=[];requests=[]
    try:
        with tempfile.TemporaryDirectory() as tmp:
            private=pathlib.Path(tmp)
            # Recipe reconstructs object topology, not mechanical time.
            initial=old.copy(deep=True);mask=initial.is_core_node.values!=0
            initial.soil_mass_kg.values[:]=np.where(mask,c.MASS0/int(mask.sum()),0.)
            state,patch=pd.setup(initial,seed)
            REQ(np.array_equal(patch,np.asarray(recipe[0])),"initial recipe")
            plan={sid:xyz for sid,xyz in state["plan"]}
            plan_sha=h41.digest(np.asarray([xyz for _,xyz in state["plan"]]))
            REQ(plan_sha==ds.attrs["p0d_plan_sha"],"source recipe")
            state["release"]=list(zip(ds.release_step.values.tolist(),ds.source_id.values.tolist()))
            for _,sid in state["release"]:h35.add_grain(state,sid,plan[sid])
            state["plan"]=[];state["last_release"]=state["release"][-1][0]
            for vertices in recipe[1:]:pd.replace(state,np.asarray(vertices))
            REQ(state["ids"]==ds.native_id.values.tolist() and state["floor"]==meta["floor"],"native recipe IDs")
            p.restoreState(fileName=str(folder/"native.bullet"),physicsClientId=state["cid"])
            REQ(np.array_equal(h41.snap(state["cid"],state["ids"]),ds.grain_state.values),"restore exact")
            pd.integrity(state)
            form=int(ds.attrs["p0d_formation_steps"]);patch=np.asarray(recipe[-1])
            prev=private/"block_9";prev.mkdir();old.to_netcdf(prev/"soil.nc",engine="scipy")
            REQ(c.hash_owner(b.read_owner(prev/"soil.nc"))==c.hash_owner(old),"checkpoint projection")
            for epoch in (4,5,6):
                soil,event,arrays=live.request(python,rawpath,raw,seed,k,private,epoch,old)
                rejected(lambda:live.request(python,rawpath,raw,seed,k,private,epoch,old),"duplicate request")
                names=("before","after","erosion","deposition","influx","outflux","bedrock","receiver","order","active","core")
                REQ(set(arrays)=={f"s{j}_{name}" for j in range(3*epoch-2,3*epoch+1) for name in names},"native fields")
                red=live.unpack(owner,seed,epoch-1,"actual");tm=[]
                for j in range(3*epoch-2,3*epoch+1):
                    fs={key:arrays[f"s{j}_{key}"] for key in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
                    red,m=ker.advance(red,fs,j);tm.append(m)
                REQ(np.array_equal(red.soil_mass_kg.values,soil.soil_mass_kg.values.ravel()),"material")
                REQ(abs(float(red.attrs["mass_export_kg"])-float(soil.attrs["exported_soil_kg"]))<=.09216,"export")
                newpatch=h40.terrain(soil,seed);trace=hashlib.sha256()
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
                red,hits,rad=live.irradiate(red,base,seed,"actual")
                owner=live.publish(base,red,seed,epoch,"actual")
                owner_sha=live.validate(owner,raw,seed,k,epoch,"actual")
                REQ(pd.hash_owner(live.project(owner),soil)==base_sha,"projection")
                obs=pd.observe(base,soil,state,patch)
                pose=h41.snap(cid,state["ids"]).copy();prior=live.digest(owner)
                for _ in range(2):
                    check=live.unpack(owner,seed,epoch,"actual");REQ(np.isfinite(thermal.deltaT(check)).all(),"energy query")
                REQ(live.digest(owner)==prior and np.array_equal(pose,h41.snap(cid,state["ids"])),"query mutation")
                path=out/f"epoch_{epoch}.nc";owner.to_netcdf(path,engine="scipy")
                REQ(live.validate(b.read_owner(path),raw,seed,k,epoch,"actual")==owner_sha,"roundtrip")
                np.savez_compressed(out/f"flux_{epoch}.npz",**arrays)
                np.savez_compressed(out/f"ray_{epoch}.npz",source_id=np.where(hits>=0,red.grain.values[np.maximum(hits,0)],-1))
                history.append(trace_sha);requests.append(event)
                epochs.append({"epoch":epoch,"steps":form+epoch*pd.BLOCK,"owner_sha":owner_sha,"projection_sha":base_sha,"soil_sha":c.hash_owner(soil),"trace_sha":trace_sha,"flux_sha":live.flux_digest(arrays),"ray_sha":hashlib.sha256(hits.tobytes()).hexdigest(),"radiation":rad,"transport":tm,"observation":obs})
                old=soil
            REQ([e["to"] for e in requests]==[12,15,18] and [e["from"] for e in requests]==[10,13,16],"suffix only")
            REQ(before=={q.name:sha(q) for q in folder.iterdir()},"checkpoint mutated")
            return {"epochs":epochs,"history":history,"formation_steps":form,"new_bullet_steps":3*pd.BLOCK,"new_thermal_steps":18,"new_soil_steps":9,"requests":requests,"controls":tests,"checkpoint_owner_sha":meta["owner_sha"],"wall_seconds":time.monotonic()-start,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    finally:
        if state is not None:p.disconnect(physicsClientId=state["cid"])

def main(python):
    OUT.mkdir();report={"experiment":"P1-C.2b-consumer","cases":{},"failed_preregistered_gates":[]}
    fail=report["failed_preregistered_gates"]
    try:
        REQ(not (ROOT/".git").exists() and not (ROOT/"evidence/P1-C/C2a/14977708755.zip").exists() and not (ROOT/"evidence/P1-C/C2a/result.json").exists(),"future data visible")
        REQ({str(q.relative_to(ROOT)) for q in (ROOT/"evidence").rglob("*") if q.is_file()}=={str(CHECK.relative_to(ROOT)),str(RAW.relative_to(ROOT))},"unexpected evidence input")
        report["runtime"]={"bullet":profile.runtime(sys.executable),"landlab":profile.runtime(python),"python":platform.python_version()}
        REQ(sha(CHECK)==CHECK_SHA and sha(RAW)==RAW_SHA,"input archive digest")
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp);rawpath=root/"raw.npz"
            with zipfile.ZipFile(RAW) as z:rawpath.write_bytes(z.read(thermal.find(z,"worldengine_two_seeds.npz")))
            raws=h37.read_raw(rawpath)
            with zipfile.ZipFile(CHECK) as z:
                REQ(z.testzip() is None,"checkpoint CRC")
                names={n for n in z.namelist() if not n.endswith("/")}
                REQ(names=={key+"/"+f for key in KEYS for f in ("owner.nc","native.bullet","manifest.json")},"checkpoint archive contents")
                z.extractall(root)
            for key in KEYS:
                seed,k=key.split("_");seed=int(seed);k=float(k)
                try:report["cases"][key]=run(rawpath,python,raws[seed],seed,k,root/key,OUT/key)
                except Exception as e:traceback.print_exc();fail.append(key+":"+repr(e))
        REQ(sha(CHECK)==CHECK_SHA and sha(RAW)==RAW_SHA,"input archive mutated")
    except Exception as e:traceback.print_exc();fail.append("instrument:"+repr(e))
    report["outputs"]={str(q.relative_to(OUT)):sha(q) for q in OUT.rglob("*") if q.is_file()}
    (OUT/"result.json").write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("WORLD_EVAL_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(*sys.argv[1:])
