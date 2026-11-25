"""P3.1: measured link fluxes, finite infiltration storage, same-host restart.
No erosion or independent water inventory. Native Landlab 2.11.
"""
import hashlib,json,pathlib,resource,sys,time,traceback
from importlib.metadata import version
import numpy as np
import xarray as xr
from landlab import RasterModelGrid,NodeStatus
from landlab.components import OverlandFlow,SoilInfiltrationGreenAmpt

DX=2.;DT=.25;SHAPE=(12,12);N=144;STEPS=960;MID=480
FILM=1e-8;THICK=.2;RHO=1590.;BASE=.03;CAP=.08
CASES=("active","closed","noinf","norain","saturated")
KEYS=("step","rain_m3","out_m3","in_m3","surface_m3","subsurface_m3","reject_m3","residual_m3","max_flow_local_m","max_inf_local_m","min_h_m","max_h_m","q_abs_max_m2s","hash")
def req(ok,msg):
    if not bool(ok):raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(v):return v.item() if isinstance(v,np.generic) else v
def digest(ds):
    h=hashlib.sha256(json.dumps({k:norm(v) for k,v in ds.attrs.items()},sort_keys=True).encode())
    for n in sorted(ds.variables):
        a=ds[n];v=a.values
        h.update(json.dumps([n,list(a.dims),dict(a.attrs)],sort_keys=True).encode())
        h.update(np.asarray(v,dtype="<f8" if v.dtype.kind=="f" else "<i8").tobytes())
    return h.hexdigest()
def rate(case,step):
    if case=="norain":return 0.
    return (72. if step<=240 else 36. if 480<step<=720 else 0.)/3.6e6
def geometry(case):
    g=RasterModelGrid(SHAPE,xy_spacing=DX)
    g.set_closed_boundaries_at_grid_edges(True,True,True,True)
    if case!="closed":g.status_at_node[5]=NodeStatus.FIXED_VALUE
    core=g.status_at_node==NodeStatus.CORE
    bed=1.+.005*g.node_y+.001*g.node_x+.02*np.exp(-((g.node_x-11.)**2+(g.node_y-11.)**2)/16.)
    mass=np.where(core,RHO*THICK*DX*DX,0.)
    return g,core,bed,mass
def initial(case):
    g,core,bed,mass=geometry(case)
    s=np.where(core,CAP if case=="saturated" else .04,0.)
    h=np.where(core,.005,FILM)
    ds=xr.Dataset({
        "bedrock_m":("node",bed,{"units":"m_nominal"}),
        "soil_mass_kg":("node",mass,{"units":"kg"}),
        "node_status":("node",g.status_at_node.astype(np.int32)),
        "core_mask":("node",core.astype(np.int32)),
        "node_x_m":("node",g.node_x,{"units":"m_nominal"}),
        "node_y_m":("node",g.node_y,{"units":"m_nominal"}),
        "tail":("link",g.node_at_link_tail.astype(np.int32)),
        "head":("link",g.node_at_link_head.astype(np.int32)),
        "surface_depth_m":("node",h,{"units":"m"}),
        "subsurface_water_depth_m":("node",s,{"units":"m"}),
        "link_q_m2_s":("link",np.zeros(g.number_of_links),{"units":"m2/s"}),
    },coords={"node":np.arange(g.number_of_nodes,dtype=np.int32),"link":np.arange(g.number_of_links,dtype=np.int32)},
    attrs={"schema":"P3.1-v1","case":case,"dx_m":DX,"dt_s":DT,"step":0,"time_s":0.,"rain_m3":0.,"out_m3":0.,"in_m3":0.,"capacity_rejected_m3":0.,"initial_water_m3":float((h[core]+s[core]).sum()*DX*DX),"landlab_version":version("landlab"),"soil_capacity_m":CAP})
    validate(ds)
    return ds
def validate(ds,expected=None):
    case=ds.attrs["case"];req(case in CASES,"case")
    g,core,bed,mass=geometry(case)
    fields={"bedrock_m","soil_mass_kg","node_status","core_mask","node_x_m","node_y_m","tail","head","surface_depth_m","subsurface_water_depth_m","link_q_m2_s"}
    attrs={"schema","case","dx_m","dt_s","step","time_s","rain_m3","out_m3","in_m3","capacity_rejected_m3","initial_water_m3","landlab_version","soil_capacity_m"}
    req(set(ds.data_vars)==fields and set(ds.coords)=={"node","link"} and set(ds.attrs)==attrs,"schema or duplicated stock")
    req(ds.attrs["schema"]=="P3.1-v1" and ds.attrs["dx_m"]==DX and ds.attrs["dt_s"]==DT and ds.attrs["soil_capacity_m"]==CAP and ds.attrs["landlab_version"]=="2.11.0","metadata")
    req(np.array_equal(ds.node.values,np.arange(N)) and np.array_equal(ds.link.values,np.arange(g.number_of_links)),"coordinates")
    for k,v in {"bedrock_m":bed,"soil_mass_kg":mass,"node_status":g.status_at_node,"core_mask":core,"node_x_m":g.node_x,"node_y_m":g.node_y,"tail":g.node_at_link_tail,"head":g.node_at_link_head}.items():
        req(np.array_equal(ds[k].values,v),"fixed geometry "+k)
    for k,u in {"bedrock_m":"m_nominal","soil_mass_kg":"kg","node_x_m":"m_nominal","node_y_m":"m_nominal","surface_depth_m":"m","subsurface_water_depth_m":"m","link_q_m2_s":"m2/s"}.items():req(ds[k].attrs=={"units":u},"units "+k)
    h=ds.surface_depth_m.values;s=ds.subsurface_water_depth_m.values;q=ds.link_q_m2_s.values
    req(h.shape==s.shape==(N,) and q.shape==(g.number_of_links,) and ds.surface_depth_m.dims==ds.subsurface_water_depth_m.dims==("node",) and ds.link_q_m2_s.dims==("link",),"shape")
    req(np.isfinite(h).all() and np.isfinite(s).all() and np.isfinite(q).all(),"finite")
    req(np.all(h[core]>=0) and np.all(h[~core]==FILM) and np.all(s[~core]==0) and np.all(s[core]>= (CAP if case=="saturated" else .04)) and np.all(s[core]<=CAP),"water bounds")
    inactive=np.ones(len(q),dtype=bool);inactive[g.active_links]=False
    req(np.all(q[inactive]==0),"inactive q")
    step=int(ds.attrs["step"])
    req(ds.attrs["step"]==step and 0<=step<=STEPS and ds.attrs["time_s"]==step*DT and (expected is None or step==expected),"clock")
    expected_rain=sum(rate(case,i)*DT*DX*DX*int(core.sum()) for i in range(1,step+1))
    req(abs(ds.attrs["rain_m3"]-expected_rain)<=1e-10,"rain schedule")
    init=(.005+(CAP if case=="saturated" else .04))*DX*DX*int(core.sum())
    req(abs(ds.attrs["initial_water_m3"]-init)<=1e-12,"initial ledger")
    req(all(np.isfinite(ds.attrs[k]) and ds.attrs[k]>=0 for k in ("rain_m3","out_m3","in_m3","capacity_rejected_m3")),"ledger")
    residual=float((h[core]+s[core]).sum()*DX*DX+ds.attrs["out_m3"]-ds.attrs["in_m3"]-ds.attrs["rain_m3"]-ds.attrs["initial_water_m3"])
    req(abs(residual)<=1e-8,"water balance "+repr(residual))
    if case=="closed":req(ds.attrs["out_m3"]==ds.attrs["in_m3"]==0,"closed boundary")
    if case in ("noinf","saturated"):req(np.all(s[core]==(CAP if case=="saturated" else .04)),"infiltration control")
    return residual
class Model:
    def __init__(self,ds):
        self.case=ds.attrs["case"];self.g,self.core,bed,mass=geometry(self.case)
        g=self.g
        g.at_node["topographic__elevation"]=bed+THICK
        g.at_node["surface_water__depth"]=np.zeros(N)
        g.at_node["soil_water_infiltration__depth"]=np.full(N,.01)
        self.water=OverlandFlow(g,h_init=FILM,alpha=.7,theta=.8,mannings_n=.03,steep_slopes=False,rainfall_intensity=0.)
        self.inf=SoilInfiltrationGreenAmpt(g,hydraulic_conductivity=np.where(self.core,0. if self.case=="noinf" else 1e-6,0.),soil_bulk_density=RHO,rock_density=2650.,initial_soil_moisture_content=.15,soil_type="sandy loam",surface_water_minimum_depth=FILM)
        self.h=g.at_node["surface_water__depth"];self.q=g.at_link["surface_water__discharge"];self.F=g.at_node["soil_water_infiltration__depth"]
        self.s=ds.subsurface_water_depth_m.values.copy()
        np.copyto(self.h,ds.surface_depth_m.values);np.copyto(self.q,ds.link_q_m2_s.values)
        self.F[self.core]=self.s[self.core]-BASE
        self.owner=ds.copy(deep=True)
        t,h=g.node_at_link_tail,g.node_at_link_head
        self.bsign=self.core[t].astype(int)-self.core[h].astype(int)
        self.bsign[~np.isin(np.arange(g.number_of_links),g.active_links)]=0
        self.area=DX*DX
        req(int(np.count_nonzero(self.bsign))==(0 if self.case=="closed" else 1),"boundary link count")
    def advance(self,step):
        g,core=self.g,self.core;old=self.h.copy();rain=rate(self.case,step)
        req(self.owner.attrs["step"]==step-1 and DT<=self.water.calc_time_step(),"step/CFL")
        self.water.rainfall_intensity=rain
        req(self.water.run_one_step(DT)==DT,"duration")
        # q per unit width, not volume. Signed native boundary crossing.
        q=self.q;post=self.h.copy();div=g.calc_flux_div_at_node(q)
        local=float(np.max(np.abs((post-old-(rain-div)*DT)[core])))
        req(local<=1e-10 and np.all(post[core]>=0),"flow local or negative "+repr(local))
        signed=q*self.bsign*DX*DT
        out=float(signed[signed>0].sum());inc=float(-signed[signed<0].sum())
        # Check separately from actual link flux, not closing by residual.
        req(abs(float(div[core].sum()*self.area*DT)-(out-inc))<=1e-10,"boundary/divergence")
        f0=self.F.copy();s0=self.s.copy()
        self.inf.run_one_step(DT)
        raw=self.F-f0;raw_h=post-self.h
        inf_local=float(np.max(np.abs(raw-raw_h)))
        req(inf_local<=1e-12 and np.all(raw>=0) and np.all(raw[~core]==0) and np.all(self.h[core]>=0),"native infiltration transfer")
        # Native proposed transfer accepted only into remaining pore capacity.
        accepted=np.minimum(raw[core],CAP-s0[core])
        req(np.all(accepted>=0),"capacity")
        self.s[core]=np.minimum(s0[core]+accepted,CAP)
        actual=self.s[core]-s0[core]
        rejected=float((raw[core]-actual).sum()*self.area)
        req(rejected>=-1e-12 and np.all(np.abs(actual-accepted)<=1e-12),"cap transfer")
        # Compute rejection as proposed minus bounded accepted, not budget correction.
        self.h[core]=post[core]-actual
        self.F[core]=self.s[core]-BASE
        req(np.all(self.F[~core]==.01) and np.all(np.abs((post[core]-self.h[core])-actual)<=1e-12),"adaptation transfer")
        a=self.owner.attrs
        a["step"]=step;a["time_s"]=step*DT
        a["rain_m3"]+=rain*DT*int(core.sum())*self.area
        a["out_m3"]+=out;a["in_m3"]+=inc
        # Store nonnegative capacity-limited proposal separately; arithmetic
        # error is not a source or sink of water.
        cap_reject=float(np.maximum(raw[core]-accepted,0).sum()*self.area)
        a["capacity_rejected_m3"]+=cap_reject
        np.copyto(self.owner.surface_depth_m.values,self.h)
        np.copyto(self.owner.subsurface_water_depth_m.values,self.s)
        np.copyto(self.owner.link_q_m2_s.values,self.q)
        # Full immutable geometry/schema validated at boundaries; per-step
        # conservation/finite/physical bounds checked directly.
        req(np.isfinite(self.h).all() and np.isfinite(self.q).all() and np.isfinite(self.s).all() and np.all(self.h[core]>=0) and np.all(self.h[~core]==FILM) and np.all(self.s[core]<=CAP) and np.all(self.s[core]>=s0[core]),"water bounds")
        req(np.all(q[g.status_at_link!=0]==0),"inactive q")
        surf=float(self.h[core].sum()*self.area);sub=float(self.s[core].sum()*self.area)
        residual=surf+sub+a["out_m3"]-a["in_m3"]-a["rain_m3"]-a["initial_water_m3"]
        req(abs(residual)<=1e-8,"water budget step "+str(step)+" "+repr(residual))
        if self.case=="closed":req(out==inc==0,"closed")
        if self.case in ("noinf","saturated"):req(np.array_equal(s0,self.s),"control storage")
        return dict(zip(KEYS,[step,a["rain_m3"],a["out_m3"],a["in_m3"],surf,sub,a["capacity_rejected_m3"],residual,local,inf_local,float(self.h[core].min()),float(self.h[core].max()),float(np.max(np.abs(q))),digest(self.owner)]))
def snapshot(ds,path):
    validate(ds);h=digest(ds)
    ds.to_netcdf(path,engine="scipy")
    with xr.open_dataset(path,engine="scipy") as f:back=f.load()
    validate(back);req(digest(back)==h,"roundtrip")
    return {"sha256":sha(path),"owner_sha":h}
def read(path,meta=None,step=None):
    if meta is not None:req(sha(path)==meta["sha256"],"file digest")
    with xr.open_dataset(path,engine="scipy") as f:ds=f.load()
    validate(ds,step)
    if meta is not None:req(digest(ds)==meta["owner_sha"],"logical digest")
    return ds
def series(rows,path):
    np.savez_compressed(path,**{k:np.asarray([r[k] for r in rows]) for k in KEYS})
def obs(ds):
    h=digest(ds);core=ds.core_mask.values!=0
    info={"surface_m3":float(ds.surface_depth_m.values[core].sum()*DX*DX),"subsurface_m3":float(ds.subsurface_water_depth_m.values[core].sum()*DX*DX),"q_max":float(np.abs(ds.link_q_m2_s.values).max())}
    req(digest(ds)==h,"observer mutated")
    return info
def producer(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir()
    start=time.monotonic();report={"experiment":"P3.1","cases":{},"failed_preregistered_gates":[]}
    fail=report["failed_preregistered_gates"];final={}
    for case in CASES:
        folder=out/case;folder.mkdir();rows=[]
        try:
            ds=initial(case);snap0=snapshot(ds,folder/"initial.nc");m=Model(ds)
            for step in range(1,STEPS+1):
                rows.append(m.advance(step))
                if step==MID:
                    mid=snapshot(m.owner,folder/"mid.nc");obs(m.owner)
                    if case=="active":
                        meta=snapshot(m.owner,cp/"owner.nc")
                        req(meta["owner_sha"]==mid["owner_sha"],"checkpoint identity")
                        (cp/"manifest.json").write_text(json.dumps({"step":MID,"case":case,**meta},sort_keys=True)+"\n")
            end=snapshot(m.owner,folder/"final.nc");obs(m.owner)
            series(rows,folder/"series.npz")
            final[case]=m.owner.copy(deep=True)
            report["cases"][case]={"initial":snap0,"mid":mid,"final":end,"last":rows[-1],"max_residual_m3":max(abs(x["residual_m3"]) for x in rows),"max_flow_local_m":max(x["max_flow_local_m"] for x in rows),"max_inf_local_m":max(x["max_inf_local_m"] for x in rows),"min_h_m":min(x["min_h_m"] for x in rows),"steps":len(rows)}
        except Exception as e:
            traceback.print_exc();fail.append(case+":"+repr(e))
            if rows:series(rows,folder/"partial.npz")
    if len(final)==5:
        a=report["cases"]["active"]["last"]
        def gate(k,x):
            if not bool(x):fail.append(k)
        gate("active_out",a["out_m3"]>.01);gate("active_infiltration",a["subsurface_m3"]-16.>.01)
        gate("closed",report["cases"]["closed"]["last"]["out_m3"]==report["cases"]["closed"]["last"]["in_m3"]==0)
        gate("noinf",report["cases"]["noinf"]["last"]["subsurface_m3"]==report["cases"]["noinf"]["last"]["subsurface_m3"] and np.array_equal(final["noinf"].subsurface_water_depth_m,initial("noinf").subsurface_water_depth_m))
        gate("norain",report["cases"]["norain"]["last"]["rain_m3"]==0)
        gate("saturated",report["cases"]["saturated"]["last"]["reject_m3"]>0 and np.array_equal(final["saturated"].subsurface_water_depth_m,initial("saturated").subsurface_water_depth_m))
        for k in ("closed","noinf"):
            d=float(np.max(np.abs(final[k].surface_depth_m.values-final["active"].surface_depth_m.values)))
            report[k+"_difference_m"]=d;gate(k+"_difference",d>1e-6)
    report["wall_seconds"]=time.monotonic()-start;report["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_WATER_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
def consumer(checkpoint,dest):
    cp=pathlib.Path(checkpoint);out=pathlib.Path(dest);out.mkdir()
    meta=json.loads((cp/"manifest.json").read_text());req(meta["step"]==MID and meta["case"]=="active","checkpoint metadata")
    ds=read(cp/"owner.nc",meta,MID)
    controls={}
    def reject(name,fn):
        try:fn()
        except ValueError:controls[name]=True
        else:raise ValueError("accepted corruption "+name)
    bad=out/"bad.nc";b=bytearray((cp/"owner.nc").read_bytes());b[-1]^=1;bad.write_bytes(b)
    reject("file",lambda:read(bad,meta,MID));bad.unlink()
    alt=ds.copy(deep=True);g,core,_,_=geometry("active")
    alt.link_q_m2_s.values[g.active_links[0]]+=1e-4
    def verify(x):
        validate(x,MID);req(digest(x)==meta["owner_sha"],"logical digest")
    reject("q",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.subsurface_water_depth_m.values[np.flatnonzero(core)[0]]=CAP+.01
    reject("capacity",lambda:verify(alt))
    alt=ds.copy(deep=True);alt["duplicate_water"]=alt.surface_depth_m.copy()
    reject("duplicate",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.attrs["step"]+=1
    reject("clock",lambda:verify(alt))
    m=Model(ds);req(digest(m.owner)==meta["owner_sha"],"restored state")
    req(np.array_equal(m.q,ds.link_q_m2_s.values) and np.array_equal(m.h,ds.surface_depth_m.values),"restored native fields")
    rows=[m.advance(step) for step in range(MID+1,STEPS+1)]
    fin=snapshot(m.owner,out/"final.nc");series(rows,out/"series.npz")
    result={"controls":controls,"input_owner_sha":meta["owner_sha"],"final":fin,"steps":len(rows),"first_step":rows[0]["step"],"last":rows[-1],"outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("WORLD_WATER_RESTART_JSON "+json.dumps(result,sort_keys=True),flush=True)
def compare(reference,consumer_path,dest):
    p=pathlib.Path(reference);c=pathlib.Path(consumer_path)
    a=json.loads((p/"result.json").read_text());b=json.loads((c/"result.json").read_text())
    req(not a["failed_preregistered_gates"] and len(b["controls"])==5 and all(b["controls"].values()),"source/control")
    req(b["steps"]==STEPS-MID and b["first_step"]==MID+1 and b["input_owner_sha"]==a["cases"]["active"]["mid"]["owner_sha"],"restart suffix")
    req(b["final"]["owner_sha"]==a["cases"]["active"]["final"]["owner_sha"],"restart owner")
    with np.load(p/"active/series.npz",allow_pickle=False) as x,np.load(c/"series.npz",allow_pickle=False) as y:
        req(set(x.files)==set(y.files)==set(KEYS),"series fields")
        req(all(np.array_equal(x[k][MID:],y[k]) for k in KEYS),"restart series")
    req(digest(read(p/"active/final.nc"))==digest(read(c/"final.nc")),"final state")
    result={"suffix_exact":True,"steps":STEPS-MID,"controls":b["controls"],"final_owner_sha":b["final"]["owner_sha"]}
    pathlib.Path(dest).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("WORLD_WATER_COMPARE_JSON "+json.dumps(result,sort_keys=True),flush=True)
if __name__=="__main__":
    {"produce":producer,"consume":consumer,"compare":compare}[sys.argv[1]](*sys.argv[2:])
