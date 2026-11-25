"""E7.1b: prospective feed-order/force redistribution on retained exogenous soil."""
import hashlib,json,pathlib,resource,sys,tempfile,time,traceback,zipfile
from importlib.metadata import version
import numpy as np
import pybullet as p
import xarray as xr
import p1c_live_world as live
from p1c_live_world import pd,h35,h26,h40,h41,h37,c,thermal,profile
from e7_force_instrument import contact_rows,cdf,l1,GRID,G
ROOT=live.ROOT
P=ROOT/"evidence/P1-C/C2a/14977708755.zip"
RAW=ROOT/"evidence/archives/45389348480/14917807888.zip"
P_SHA="fd0ba07f97a2e0f9121bc5cdf6b19ba37d2dcd06a1c751a193fa60d40b271d7d"
RAW_SHA="f2b204f42f7d055600587230093ad36b1a89f0f85ae758bae114a4b6b94a5170"
REQ=live.REQ
N,R,M=pd.N,pd.R,pd.M
CENTER=np.array([40.,30.])
SEEDS=(7321,8128)
FEEDS=(2026112501,2026112502)
MODES=("alternate","grouped")
TEMPLATES=("monolayer","pyramid")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def array_sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def key(t,f,k,m):return f"{t}_{f}_{k}_{m}"
def load_inputs():
    REQ(sha(P)==P_SHA and sha(RAW)==RAW_SHA,"archive digest")
    with tempfile.TemporaryDirectory() as tmp:
        path=pathlib.Path(tmp)/"raw.npz"
        with zipfile.ZipFile(RAW) as z:
            REQ(z.testzip() is None,"raw CRC")
            path.write_bytes(z.read(thermal.find(z,"worldengine_two_seeds.npz")))
        raws=h37.read_raw(path)
    soils={}
    with zipfile.ZipFile(P) as z:
        REQ(z.testzip() is None,"producer CRC")
        ref=json.loads(z.read("result.json"));REQ(not ref["failed_preregistered_gates"],"producer gates")
        for t in SEEDS:
            for k in (0.,.1):
                history=[]
                for epoch,rec in enumerate(ref["cases"][f"{t}_{k}"]["epochs"]):
                    n=f"{t}_{k}/epoch_{epoch}.nc"
                    REQ(hashlib.sha256(z.read(n)).hexdigest()==ref["outputs"][n],"reference file")
                    owner=thermal.read_nc(z,n)
                    REQ(live.digest(owner)==rec["owner_sha"],"reference owner")
                    base=live.project(owner);soil=pd.extract_soil(base)
                    pd.validate(base,soil,raws[t],t,k,epoch)
                    REQ(pd.hash_owner(base,soil)==rec["projection_sha"] and c.hash_owner(soil)==rec["soil_sha"],"reference projection")
                    history.append(soil)
                REQ(len(history)==7,"soil epochs")
                if k==0:
                    REQ(all(np.array_equal(s.soil_mass_kg,history[0].soil_mass_kg) and np.array_equal(c.topography(s),c.topography(history[0])) and s.attrs["exported_soil_kg"]==0 for s in history),"K0 soil")
                soils[t,k]=history
    return raws,soils
def feed_points(seed,patch,mode):
    u=np.random.default_rng(seed).random((N,2))
    rad=.17*np.sqrt(u[:,0]);angle=2*np.pi*u[:,1]
    pts=np.column_stack([np.where(np.arange(N)<63,39.55,40.45)+rad*np.cos(angle),30+rad*np.sin(angle),np.full(N,h40.height(patch,*CENTER)+1.5)])
    order=list(range(N)) if mode=="grouped" else [j for i in range(62) for j in (i,i+63)]+[62]
    REQ(sorted(order)==list(range(N)) and pts.shape==(N,3),"feed identity")
    return pts,order
def template_points(patch,mode):
    d=.152;loc=[]
    if mode=="monolayer":
        pnts=[(i-5.,j-5.5,i,j) for i in range(11) for j in range(12)]
        pnts=sorted(pnts,key=lambda x:(x[0]**2+x[1]**2,x[2],x[3]))[:N]
        loc=[[x*d,y*d,.077] for x,y,i,j in pnts]
    else:
        for lev,n in enumerate((7,6,5,4)):
            for i in range(n):
                for j in range(n):
                    if lev==3 and i==j==3:continue
                    loc.append([(i-(n-1)/2)*d,(j-(n-1)/2)*d,.077+lev*d/np.sqrt(2)])
    loc=np.asarray(loc,dtype=float);REQ(loc.shape==(N,3),"template inventory")
    loc[:,:2]-=loc[:,:2].mean(axis=0)
    height=h40.height(patch,*CENTER)
    a=h40.height(patch,CENTER[0]+1,CENTER[1])-height
    b=h40.height(patch,CENTER[0],CENTER[1]+1)-height
    n=np.array([-a,-b,1.]);n/=np.linalg.norm(n)
    u=np.array([1.,0.,a]);u/=np.linalg.norm(u)
    v=np.cross(n,u)
    origin=np.array([*CENTER,height]);pts=origin+loc@np.stack([u,v,n])
    ds=np.linalg.norm(pts[:,None,:]-pts[None,:,:],axis=2);np.fill_diagonal(ds,np.inf)
    REQ(np.min(ds)>=2*R and np.all((pts-origin)@n>=R) and np.all(pts[:,0]>pts[:,1]) and np.all((pts[:,:2]>2)&(pts[:,:2]<98)),"template placement")
    return pts
class Window:
    def __init__(self):
        self.force=[];self.speeds=[];self.react=np.zeros((2,len(GRID)));self.mass=np.zeros_like(self.react);self.h=hashlib.sha256()
    def sample(self,state,patch,step):
        cid=state["cid"];ids=state["ids"];q=h41.snap(cid,ids)
        rows=contact_rows(cid,state["floor"],ids,patch)
        if len(self.force) in (0,239):REQ(np.array_equal(q,h41.snap(cid,ids)),"observation mutated")
        REQ(len(rows)>0 and np.all(rows[:,12]>=0),"vertical reaction CDF weights")
        f=rows[:,10:13].sum(axis=0);REQ(f[2]>0,"positive support")
        self.force.append(f);self.speeds.append(np.linalg.norm(q[:,3:6],axis=1))
        com=q[:,:2].mean(axis=0)
        for j in range(2):
            self.react[j]+=cdf(rows[:,1+j]-com[j],rows[:,12])
            self.mass[j]+=cdf(q[:,j]-com[j],np.full(N,M*G))
        self.h.update(np.array([step,len(rows),len(q)],dtype=np.int64).tobytes());self.h.update(rows.tobytes());self.h.update(q.tobytes())
    def finish(self):
        REQ(len(self.force)==len(self.speeds)==240,"window length")
        f=np.asarray(self.force);v=np.asarray(self.speeds);mean=f.mean(axis=0)
        res=float(np.linalg.norm(mean-[0,0,N*M*G])/(N*M*G))
        p90=float(np.quantile(v,.9));last=float(np.quantile(v[-1],.9))
        a=self.react/240;b=self.mass/240;delta=a-b
        info={"valid":bool(res<=.01 and p90<=.08 and last<=.08),"reaction_N":mean.tolist(),"relative_residual":res,"p90_speed":p90,"last_p90_speed":last,"redistribution_x":l1(delta[0]),"redistribution_y":l1(delta[1]),"window_sha":self.h.hexdigest()}
        return info,{"forces":f,"speeds":v,"reaction_cdf":a,"mass_cdf":b,"R":delta}
def structure(state,source_ids,patch,initial_edges):
    cid=state["cid"];ids=state["ids"];q=h41.snap(cid,ids)
    lookup=dict(zip(ids,source_ids));edges=set()
    for x in p.getContactPoints(physicsClientId=cid):
        if x[1] in lookup and x[2] in lookup and x[1]!=x[2] and x[9]>1e-5:
            edges.add(tuple(sorted((lookup[x[1]],lookup[x[2]]))))
    parent=list(range(N))
    def root(i):
        while parent[i]!=i:i=parent[i]
        return i
    for i,j in sorted(edges):parent[root(i)]=root(j)
    counts=np.bincount([root(i) for i in range(N)],minlength=N)
    eig=np.linalg.eigvalsh(np.cov(q[:,:2],rowvar=False,bias=True))
    REQ(eig.sum()>0,"footprint covariance")
    high=sum(float(row[2])-h40.height(patch,*row[:2])>.26 for row in q)
    if initial_edges is None:initial_edges=edges.copy()
    union=initial_edges|edges
    return {"high":int(high),"anisotropy":float((eig[1]-eig[0])/eig.sum()),"largest_component":int(counts.max()),"edge_count":len(edges),"edge_jaccard":len(initial_edges&edges)/len(union) if union else 1.},q,np.asarray(sorted(edges),dtype=np.int32).reshape(-1,2),initial_edges
def template_owner(soil,state,t,k,mode,epoch,initial):
    ds=soil.copy(deep=True).assign_coords(grain=np.arange(N,dtype=np.int32),motion=np.arange(13,dtype=np.int32))
    ds["grain_state"]=(("grain","motion"),h41.snap(state["cid"],state["ids"]))
    ds["grain_mass_kg"]=(("grain",),np.full(N,M),{"units":"kg"})
    ds["grain_radius_m_nominal"]=(("grain",),np.full(N,R),{"units":"m_nominal"})
    ds["source_id"]=(("grain",),np.arange(N,dtype=np.int32))
    ds["native_id"]=(("grain",),np.asarray(state["ids"],dtype=np.int32))
    ds["grain_reservoir_kg"]=xr.DataArray(0.,attrs={"units":"kg"})
    ds.attrs.update(e7_template=mode,e7_epoch=epoch,e7_steps=1800+epoch*720,e7_initial_sha=initial,e7_geometry="nominal")
    REQ(all(ds[k].identical(soil[k]) for k in soil.variables) and abs(float(ds.grain_mass_kg.sum())-N*M)<1e-12,"template stock")
    return ds
def scene(folder,raw,soils,t,k,feed,mode):
    start=time.monotonic();folder.mkdir();state=None;old=h41.GRAIN_SEEDS[t]
    is_template=mode in TEMPLATES
    try:
        if not is_template:h41.GRAIN_SEEDS[t]=feed
        state=h35.setup();cid=state["cid"];patch=h40.terrain(soils[0],t)
        p.removeBody(state["floor"],physicsClientId=cid);state["floor"]=h41.make_ground(cid,patch)
        if is_template:
            pts=template_points(patch,mode);order=list(range(N));state["plan"]=[]
            for i in order:h35.add_grain(state,i,pts[i])
        else:
            pts,order=feed_points(feed,patch,mode)
            state["plan"]=[(i,pts[i].copy()) for i in order]
        plan_sha=array_sha(pts[order]);initial_sha=array_sha(pts)
        win=Window()
        if is_template:
            for step in range(1,1801):
                p.stepSimulation(physicsClientId=cid)
                if step>1560:win.sample(state,patch,step)
            form=1800;source_ids=list(range(N))
        else:
            for step in range(1,h26.MAX_STEPS+1):
                done=h35.step_scene(state,step)
                if state["last_release"] is not None and step-state["last_release"]>1560:win.sample(state,patch,step)
                if done:break
            else:raise ValueError("feed timeout")
            form=step;source_ids=[sid for _,sid in state["release"]]
            REQ(not state["plan"] and sorted(source_ids)==list(range(N)) and state["max_ledger_error_kg"]<=1e-12,"feed stock")
        REQ(len(state["ids"])==N and sorted(source_ids)==list(range(N)),"grain IDs")
        pd.integrity(state)
        rank=float(np.mean([i for i,s in enumerate(source_ids) if s<63])) if not is_template else None
        side={"terrain":t,"K":k,"feed":feed,"mode":mode,"point_sha":initial_sha,"plan_sha":plan_sha,"points":pts.tolist(),"order":order,"release":state["release"],"formation_steps":form,"mean_left_rank":rank,"soil_sha":[c.hash_owner(s) for s in soils],"recipe":[h40.terrain(s,t).tolist() for s in soils]}
        (folder/"recipe.json").write_text(json.dumps(side,indent=2,sort_keys=True)+"\n")
        epochs=[];arrays={};initial_edges=None
        for epoch,soil in enumerate(soils):
            if epoch:
                patch=h40.terrain(soil,t);pd.replace(state,patch);pd.integrity(state)
                win=Window()
                for j in range(1,721):
                    p.stepSimulation(physicsClientId=cid)
                    if j>480:win.sample(state,patch,form+(epoch-1)*720+j)
            pd.integrity(state)
            info,wa=win.finish();struct,q,edges,initial_edges=structure(state,source_ids,patch,initial_edges)
            owner=template_owner(soil,state,t,k,mode,epoch,initial_sha) if is_template else pd.snapshot(soil,state,raw,t,k,epoch,form,plan_sha)
            dig=live.digest(owner);path=folder/f"owner_{epoch}.nc";owner.to_netcdf(path,engine="scipy")
            with xr.open_dataset(path,engine="scipy") as f:back=f.load()
            REQ(live.digest(back)==dig,"owner roundtrip")
            if not is_template:pd.validate(back,soil,raw,t,k,epoch)
            info.update(struct,epoch=epoch,owner_sha=dig,soil_sha=c.hash_owner(soil),steps=form+epoch*720)
            epochs.append(info)
            for n,v in {**wa,"states":q,"edges":edges}.items():arrays[f"{n}_{epoch}"]=v
            # Partial results retained after each epoch.
            np.savez_compressed(folder/"profiles.npz",**arrays)
            (folder/"epochs.json").write_text(json.dumps(epochs,indent=2,sort_keys=True)+"\n")
        return {"epochs":epochs,"recipe":side,"wall_seconds":time.monotonic()-start,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},arrays
    finally:
        h41.GRAIN_SEEDS[t]=old
        if state is not None:p.disconnect(physicsClientId=state["cid"])
def evaluate_pair(t,f,cases,profiles):
    a,b=[key(t,f,.1,m) for m in MODES]
    ts=[key(t,0,.1,m) for m in TEMPLATES]
    if not all(n in cases for n in [a,b,*ts]):return {"eligible":False,"reason":"missing case"}
    x,y=cases[a],cases[b]
    REQ(x["recipe"]["point_sha"]==y["recipe"]["point_sha"] and x["recipe"]["plan_sha"]!=y["recipe"]["plan_sha"] and x["recipe"]["soil_sha"]==y["recipe"]["soil_sha"],"paired source/soil identity")
    diff=abs(x["recipe"]["mean_left_rank"]-y["recipe"]["mean_left_rank"])
    valid=all(cases[n]["epochs"][e]["valid"] for n in [a,b,*ts] for e in (0,6))
    eligible=diff>=15 and valid
    function=[l1(profiles[a][f"R_{e}"][0]-profiles[b][f"R_{e}"][0]) for e in (0,6)]
    high=[abs(x["epochs"][e]["high"]-y["epochs"][e]["high"]) for e in (0,6)]
    anis=[abs(x["epochs"][e]["anisotropy"]-y["epochs"][e]["anisotropy"]) for e in (0,6)]
    structural=all(v>=10 for v in high) or all(v>=.15 for v in anis)
    control={}
    for n in (a,b):
        red=[cases[n]["epochs"][e]["redistribution_x"] for e in (0,6)]
        ds={q:[l1(profiles[n][f"R_{e}"][0]-profiles[q][f"R_{e}"][0]) for e in (0,6)] for q in ts}
        control[n]={"redistribution":red,"template_distances":ds,"nontrivial":all(v>.15 for v in red) and all(v>.075 for r in ds.values() for v in r)}
    pass_gate=bool(eligible and function[0]>=.15 and function[1]>=.075 and structural and any(v["nontrivial"] for v in control.values()))
    return {"eligible":eligible,"effective_rank_difference":diff,"valid_endpoints":valid,"function_distances":function,"high_differences":high,"anisotropy_differences":anis,"structural":structural,"controls":control,"passed":pass_gate}
def main(dest):
    out=pathlib.Path(dest);out.mkdir();start=time.monotonic()
    result={"experiment":"E7.1b","cases":{},"pairs":{},"failed_preregistered_gates":[],"scientific_outcome":"inconclusive"}
    profiles={};fail=result["failed_preregistered_gates"]
    try:
        result["runtime"]={"numpy":profile.runtime(sys.executable),"versions":{n:version(n) for n in ("numpy","scipy","xarray","pybullet","fipy")}}
        raw,soils=load_inputs()
        jobs=[(t,f,.1,m) for t in SEEDS for f in FEEDS for m in MODES]+[(t,0,.1,m) for t in SEEDS for m in TEMPLATES]+[(7321,FEEDS[0],0.,m) for m in MODES]
        for t,f,k,m in jobs+[(7321,FEEDS[0],.1,"alternate_replay")]:
            name=key(t,f,k,m);mode=m.removesuffix("_replay")
            try:
                v,ar=scene(out/name,raw[t],soils[t,k],t,k,f,mode)
                result["cases"][name]=v;profiles[name]=ar
                print("E7_SCENE "+json.dumps({"case":name,"wall_seconds":v["wall_seconds"],"valid_windows":sum(e["valid"] for e in v["epochs"])}),flush=True)
            except Exception as e:traceback.print_exc();fail.append(name+":"+repr(e))
        base=key(7321,FEEDS[0],.1,"alternate");rep=base+"_replay"
        if base in profiles and rep in profiles:
            a,b=result["cases"][base],result["cases"][rep]
            exact=a["recipe"]==b["recipe"] and a["epochs"]==b["epochs"] and set(profiles[base])==set(profiles[rep]) and all(x.dtype==profiles[rep][n].dtype and x.shape==profiles[rep][n].shape and x.tobytes()==profiles[rep][n].tobytes() for n,x in profiles[base].items())
            result["replay_exact"]=bool(exact)
            if not exact:fail.append("replay mismatch")
        else:fail.append("replay unavailable")
        for t in SEEDS:
            for f in FEEDS:result["pairs"][f"{t}_{f}"]=evaluate_pair(t,f,result["cases"],profiles)
        eligible=sum(bool(v.get("eligible")) for v in result["pairs"].values());passed=sum(bool(v.get("passed")) for v in result["pairs"].values())
        result["eligible_pairs"]=eligible;result["passed_pairs"]=passed
        result["scientific_outcome"]="inconclusive" if fail or eligible<3 else "positive_limited" if passed>=3 else "negative"
        REQ(sha(P)==P_SHA and sha(RAW)==RAW_SHA,"input archives mutated")
    except Exception as e:traceback.print_exc();fail.append("instrument:"+repr(e));result["scientific_outcome"]="inconclusive"
    result["wall_seconds"]=time.monotonic()-start
    result["outputs"]={str(q.relative_to(out)):sha(q) for q in out.rglob("*") if q.is_file()}
    (out/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("WORLD_E7_JSON "+json.dumps({"experiment":result["experiment"],"scientific_outcome":result["scientific_outcome"],"failed_preregistered_gates":fail,"pairs":result["pairs"],"eligible_pairs":result.get("eligible_pairs"),"passed_pairs":result.get("passed_pairs"),"replay_exact":result.get("replay_exact"),"wall_seconds":result["wall_seconds"],"case_summary":{k:{"valid_windows":sum(e["valid"] for e in v["epochs"]),"wall_seconds":v["wall_seconds"],"epochs":[{"epoch":e["epoch"],"valid":e["valid"],"redistribution_x":e["redistribution_x"],"high":e["high"],"anisotropy":e["anisotropy"],"relative_residual":e["relative_residual"],"p90_speed":e["p90_speed"]} for e in v["epochs"]]} for k,v in result["cases"].items()}},sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
