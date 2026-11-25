"""E7.1a instrument calibration only, no prospective feed scenes."""
import hashlib,json,pathlib,sys,time,traceback
import numpy as np
import pybullet as p
import probe_e6_h35_granular_midfeed_restart as h35
import probe_h41_aggregate_transit as h41
import probe_h40_live_body_soil_epochs as h40
import probe_e3_h26_gravity_feed as h26

GRID=np.linspace(-3.,3.,241)
G=abs(h26.GRAVITY)
# H26 GRAVITY is scalar, negative.
def req(ok,msg):
    if not bool(ok):raise ValueError(msg)
def cdf(x,w):
    x,w=np.asarray(x,dtype=float),np.asarray(w,dtype=float)
    req(x.ndim==1 and x.shape==w.shape and len(x)>0 and np.isfinite(x).all() and np.isfinite(w).all(),"CDF inputs")
    req(np.all((x>=GRID[0])&(x<=GRID[-1])) and np.all(w>=0) and float(w.sum())>0,"CDF bounds/weights")
    order=np.argsort(x,kind="stable")
    v=np.r_[0.,np.cumsum(w[order])]
    return v[np.searchsorted(x[order],GRID,side="right")]/v[-1]
def l1(x):return float(np.sum(np.abs(x[:-1]))*.025)
def contact_rows(cid,floor,ids,patch):
    out=[]
    for c in p.getContactPoints(bodyB=floor,physicsClientId=cid):
        req(c[2]==floor and c[1] in ids and len(c)>=14,"contact identity")
        point=np.asarray(c[6],dtype=float)
        normal=np.asarray(c[7],dtype=float)*float(c[9])
        tangent=np.asarray(c[11],dtype=float)*float(c[10])+np.asarray(c[13],dtype=float)*float(c[12])
        f=normal+tangent
        req(np.isfinite(np.r_[point,normal,tangent,f]).all() and c[9]>=0 and f[2]>=-1e-8,"invalid force")
        req(abs(point[2]-h40.height(patch,*point[:2]))<=.02,"contact point outside ground")
        out.append(np.r_[ids.index(c[1]),point,normal,tangent,f])
    return np.asarray(out,dtype=float).reshape(-1,13)
def fixture(name,patch,positions,out):
    state=h35.setup()
    cid=state["cid"]
    try:
        p.removeBody(state["floor"],physicsClientId=cid)
        floor=h41.make_ground(cid,patch)
        state["floor"]=floor
        for j,pos in enumerate(positions):h35.add_grain(state,j,np.array(pos,dtype=float))
        ids=state["ids"];n=len(ids);weight=n*h26.M*G
        req(len(set(ids))==n and p.getNumBodies(physicsClientId=cid)==n+1 and all(p.getDynamicsInfo(i,-1,physicsClientId=cid)[0]==h26.M for i in ids),"body inventory")
        req(len(contact_rows(cid,floor,ids,patch))==0,"free-contact control")
        forces=[];speeds=[];allrows=[];states=[];byid=[];tangent_norm=[]
        for step in range(1,1801):
            p.stepSimulation(physicsClientId=cid)
            if step>1560:
                before=h41.snap(cid,ids)
                rows=contact_rows(cid,floor,ids,patch)
                # Explicit per-body query crosscheck on first/last window sample.
                if step in (1561,1800):
                    flat=[]
                    for i in ids:
                        flat.extend(p.getContactPoints(bodyA=i,bodyB=floor,physicsClientId=cid))
                    req(len(flat)==len(rows) and all(c[1] in ids and c[2]==floor for c in flat),"per-body contact query")
                req(np.array_equal(before,h41.snap(cid,ids)),"force observation mutated pose")
                f=rows[:,10:13].sum(axis=0)
                req(f[2]>0 and len(rows)>0,"no positive support")
                forces.append(f)
                speeds.append(np.linalg.norm(before[:,3:6],axis=1))
                states.append(before)
                allrows.extend(np.c_[np.full(len(rows),step),rows].tolist())
                byid.append([float(rows[rows[:,0]==j,12].sum()) for j in range(n)])
                tangent_norm.append(float(np.linalg.norm(rows[:,7:10].sum(axis=0))))
        req(h41.valid_bodies(cid,ids,floor,np.full(n,h26.R)),"body integrity")
        force=np.asarray(forces);speed=np.asarray(speeds);avg=force.mean(axis=0);byid=np.asarray(byid).mean(axis=0)
        residual=float(np.linalg.norm(avg-[0,0,weight])/weight)
        req(residual<=.01,"vector reaction balance "+str(residual))
        req(float(np.quantile(speed,.9))<=.08 and float(np.quantile(speed[-1],.9))<=.08,"not stable")
        if name=="separated":req(np.all(np.abs(byid-h26.M*G)<=.01*h26.M*G),"separated grain reaction")
        if name=="stacked":req(sum(byid>1e-6)==1 and abs(max(byid)-weight)<=.01*weight,"stack support")
        if name=="inclined":req(float(np.mean(tangent_norm))>1e-6,"inclined tangential contribution")
        path=out/(name+".npz")
        np.savez_compressed(path,force=force,speed=speed,contacts=np.asarray(allrows),states=np.asarray(states),patch=patch)
        return {"grains":n,"weight_N":weight,"average_reaction_N":avg.tolist(),"relative_residual":residual,"by_grain_vertical_N":byid.tolist(),"mean_tangent_N":float(np.mean(tangent_norm)),"p90_speed":float(np.quantile(speed,.9)),"window_steps":len(force),"contact_records":len(allrows),"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
    finally:p.disconnect(physicsClientId=cid)
def main(dest):
    out=pathlib.Path(dest);out.mkdir();start=time.monotonic()
    report={"experiment":"E7.1a","fixtures":{},"failed_preregistered_gates":[]}
    try:
        m=cdf([-.5,.5],[.5,.5])
        req(l1(m-m)==0 and abs(l1(cdf([.5],[1])-m)-.5)<=1e-12,"synthetic CDF")
        try:cdf([3.01],[1])
        except ValueError:pass
        else:raise ValueError("CDF accepted outside domain")
        report["cdf_controls"]=True
        flat=np.zeros((2,2));incl=np.array([[0.,5.],[0.,5.]])
        for name,patch,points in [
            ("single",flat,[(40.,30.,.25)]),
            ("separated",flat,[(39.5,30.,.25),(40.5,30.,.25)]),
            ("stacked",flat,[(40.,30.,.25),(40.,30.,.402)]),
            ("inclined",incl,[(40.,30.,h40.height(incl,40.,30.)+.25)])]:
            try:report["fixtures"][name]=fixture(name,patch,points,out)
            except Exception as e:traceback.print_exc();report["failed_preregistered_gates"].append(name+":"+repr(e))
    except Exception as e:traceback.print_exc();report["failed_preregistered_gates"].append("instrument:"+repr(e))
    report["wall_seconds"]=time.monotonic()-start
    report["outputs"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob("*.npz")}
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_FORCE_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if report["failed_preregistered_gates"]:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
