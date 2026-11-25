"""P3.3 reciprocal hydromorphodynamics; single owner, native operators."""
import hashlib,json,pathlib,sys,time,resource,traceback
import numpy as np
import xarray as xr
from landlab.components import OverlandFlow
from p3_exner import grid,req,DX,DT,PHI,FILM,RHO_S
import zipfile
from importlib.metadata import version
def load_provider(dest):
    manifest_path=pathlib.Path("research/P3/RBD_SOURCE_MANIFEST.json")
    archive=pathlib.Path("evidence/P3/RBD_SOURCE/source.zip")
    m=json.loads(manifest_path.read_text());a=json.loads((archive.parent/"archive_manifest.json").read_text())
    req(a["commit"]==m["commit"] and sha(manifest_path)==a["source_manifest_sha256"],"source manifest identity")
    req(sha(archive)==a["source_zip_sha256"] and archive.stat().st_size==a["bytes"],"source archive identity")
    root=pathlib.Path("rbd-upstream");pkg=root/"river_bed_dynamics";pkg.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        req(z.testzip() is None and set(z.namelist())==set(m["files"])==set(a["files"]),"source archive members")
        for name,meta in m["files"].items():
            b=z.read(name)
            req(len(b)==meta["bytes"] and hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()==meta["git_blob_sha1"] and hashlib.sha256(b).hexdigest()==a["files"][name],"source member "+name)
            target=root/"LICENSE.md" if name=="LICENSE.md" else pkg/pathlib.Path(name).name
            target.write_bytes(b)
    (dest/"upstream_license.md").write_bytes((root/"LICENSE.md").read_bytes())
    (dest/"source_manifest.json").write_bytes(manifest_path.read_bytes())
    (dest/"source_archive_manifest.json").write_bytes((archive.parent/"archive_manifest.json").read_bytes())
    sys.path.insert(0,str(root.resolve()))
    from river_bed_dynamics import RiverBedDynamics
    req(version("landlab")=="2.11.0","landlab version")
    return RiverBedDynamics
from p3_pore import remap,ga,geometry,domain
from p3_water import digest,sha
A=DX*DX;STEPS=240;MID=120
CASES=("active","frozen","no_sediment","no_rain")
NODE=("V","W","H");LINK=("q","hl")
LEDGER=("rain","out","inc","solid_out","erosion","deposition","release","infiltration")
def rain(case,k):
    return (0. if case=="no_rain" else 72. if k<=60 else 36. if 120<k<=180 else 0.)/3.6e6
def initial(case):
    g,c,act,t,hd,b=grid();z=g.at_node["topographic__elevation"].copy()
    data={"bedrock":("node",z-.2),"z0":("node",z),"core":("node",c.astype(np.int32)),"active":("link",act.astype(np.int32)),"tail":("link",t.astype(np.int32)),"head":("link",hd.astype(np.int32)),"status":("node",np.asarray(g.status_at_node,dtype=np.int32)),"x":("node",g.node_x),"y":("node",g.node_y),"V":("node",np.where(c,.2*(1-PHI)*A,0.)),"W":("node",np.where(c,.04*A,0.)),"H":("node",np.where(c,.04*A,0.)),"q":("link",np.zeros(len(t))),"hl":("link",np.full(len(t),FILM))}
    ds=xr.Dataset(data,coords={"node":np.arange(len(c),dtype=np.int32),"link":np.arange(len(t),dtype=np.int32)},attrs={"schema":"P3.3-v1","case":case,"step":0,"time_s":0.,"dx_m":DX,"dt_s":DT,**{k:0. for k in LEDGER}})
    return ds
def zphys(ds):
    c=ds.core.values.astype(bool)
    return np.where(c,ds.bedrock.values+geometry(ds.V.values,A)[0],ds.z0.values)
def validate(ds,expected=None):
    case=ds.attrs["case"];req(case in CASES,"case")
    ref=initial(case)
    req(set(ds.variables)==set(ref.variables) and set(ds.attrs)==set(ref.attrs),"schema/duplicate")
    req(all(np.array_equal(ds[k].values,ref[k].values) for k in ref.variables if k not in NODE+LINK),"topology")
    req(all(ds[k].dims==ref[k].dims and ds[k].shape==ref[k].shape for k in ref.variables),"dimensions")
    req(all(ds.attrs[k]==ref.attrs[k] for k in ("schema","dx_m","dt_s","case")),"metadata")
    c=ds.core.values.astype(bool);act=ds.active.values.astype(bool)
    V,W,H,q,hl=(ds[k].values for k in NODE+LINK)
    req(all(np.isfinite(x).all() for x in (V,W,H,q,hl)) and all(np.all(x>=0) and np.all(x[~c]==0) for x in (V,W,H)),"stocks")
    req(np.all(W<=geometry(V,A)[1]+1e-12) and domain(V,W,A,c)[0],"capacity/domain")
    req(np.all(q[~act]==0) and np.all(hl[act]>0) and np.all(hl[~act]==FILM),"cache")
    k=int(ds.attrs["step"]);a=ds.attrs
    req(a["step"]==k and 0<=k<=STEPS and a["time_s"]==k*DT and (expected is None or k==expected),"clock")
    req(all(np.isfinite(a[n]) and a[n]>=0 for n in LEDGER),"ledger")
    req(abs(a["rain"]-sum(rain(case,i)*DT*A*int(c.sum()) for i in range(1,k+1)))<=1e-10,"rain schedule")
    req(a["solid_out"]==0,"solid boundary")
    sr=float(V.sum()+a["solid_out"]-ref.V.values.sum());wr=float((W+H).sum()+a["out"]-a["inc"]-a["rain"]-(ref.W.values+ref.H.values).sum())
    req(abs(sr)<=1e-9 and abs(wr)<=1e-8,"global balance")
    if case=="no_sediment":req(np.array_equal(V,ref.V.values) and all(a[n]==0 for n in ("release","erosion","deposition")),"no sediment")
    return sr,wr
def snap(ds,p):
    validate(ds);d=digest(ds);ds.to_netcdf(p,engine="scipy")
    with xr.open_dataset(p,engine="scipy") as f:back=f.load()
    req(digest(back)==d,"NC roundtrip");validate(back)
    return {"sha256":sha(p),"owner_sha":d}
def read(p,meta=None,step=None):
    if meta:req(sha(p)==meta["sha256"],"file digest")
    with xr.open_dataset(p,engine="scipy") as f:ds=f.load()
    validate(ds,step)
    if meta:req(digest(ds)==meta["owner_sha"],"owner digest")
    return ds
def observe(ds):
    d=digest(ds);v=(float(ds.V.sum()),float(ds.W.sum()),float(ds.H.sum()),float(np.max(zphys(ds))))
    req(digest(ds)==d,"observation");return v
class Model:
    def __init__(self,ds,RBD):
        self.ds=ds.copy(deep=True);self.RBD=RBD
        g,c,act,t,hd,b=grid();self.g=g;self.c=c;self.act=act;self.t=t;self.hd=hd;self.b=b
        g.at_node["surface_water__depth"]=np.zeros(len(c))
        self.of=OverlandFlow(g,h_init=FILM,alpha=.7,theta=.8,mannings_n=.03,steep_slopes=False,rainfall_intensity=0.)
        self.h=g.at_node["surface_water__depth"];self.q=g.at_link["surface_water__discharge"]
        np.copyto(self.h,np.where(c,ds.H.values/A,FILM));np.copyto(self.q,ds.q.values);np.copyto(self.of._h_links,ds.hl.values)
    def flow(self,k,freeze=False):
        ds=self.ds;c=self.c;g=self.g
        req(k==ds.attrs["step"]+1,"step")
        z=ds.z0.values if freeze else zphys(ds)
        np.copyto(g.at_node["topographic__elevation"],z)
        np.copyto(self.h,np.where(c,ds.H.values/A,FILM));req(np.array_equal(self.q,ds.q.values) and np.array_equal(self.of._h_links,ds.hl.values),"hydraulic cache")
        old=self.h.copy();prev=np.zeros_like(self.q);prev[self.act]=self.q[self.act]/ds.hl.values[self.act]
        rr=rain(ds.attrs["case"],k);self.of.rainfall_intensity=rr
        req(DT<=self.of.calc_time_step(),"CFL");req(self.of.run_one_step(DT)==DT,"duration")
        q=self.q.copy();hl=np.asarray(self.of._h_links).copy();hh=self.h.copy()
        div=g.calc_flux_div_at_node(q);local=float(np.max(np.abs((hh-old-(rr-div)*DT)[c])))
        bf=q*self.b*DX*DT;out=float(bf[bf>0].sum());inc=float(-bf[bf<0].sum())
        req(local<=1e-10 and abs(float(div[c].sum()*A*DT)-(out-inc))<=1e-10,"hydraulic local")
        req(np.all(hh>=0) and np.all(hh[~c]==FILM) and np.all(q[~self.act]==0) and np.all(hl[self.act]>0),"hydraulic bounds")
        u=np.zeros_like(q);u[self.act]=q[self.act]/hl[self.act]
        req(np.isfinite(u).all(),"velocity")
        return hh,q,hl,u,prev,rr,out,inc,local
    def advance(self,k):
        ds=self.ds;c=self.c;case=ds.attrs["case"]
        h,q,hl,u,prev,rr,out,inc,flowlocal=self.flow(k,case=="frozen")
        z=zphys(ds);V,W=ds.V.values,ds.W.values;H=np.where(c,h*A,0.)
        qb=np.zeros_like(q);nativez=z.copy()
        if case!="no_sediment":
            g,*_=grid();g.at_node["topographic__elevation"]=z.copy()
            g.at_node["surface_water__depth"]=h.copy();g.at_link["surface_water__depth"]=hl.copy();g.at_link["surface_water__velocity"]=u.copy()
            r=self.RBD(g,gsd=[[3,100],[2,50],[1,0]],rho=1000,rho_s=RHO_S,bedload_equation="MPM",lambda_p=PHI,variable_critical_shear_stress=False,outlet_boundary_condition="fixedValue",dt=DT,surface_water__velocity_prev_time_link=prev.copy(),track_stratigraphy=False)
            req(np.array_equal(r._surface_water__velocity_prev_time_link,prev),"prev velocity")
            r.run_one_step();qb=np.asarray(r._sed_transp__bedload_rate_link).copy();nativez=g.at_node["topographic__elevation"].copy()
            for x,y in ((g.at_node["surface_water__depth"],h),(g.at_link["surface_water__depth"],hl),(g.at_link["surface_water__velocity"],u),(r._surface_water__velocity_prev_time_link,prev)):
                req(np.array_equal(x,y),"RBD hydraulic mutation")
        req(np.isfinite(qb).all() and np.all(qb[~self.act]==0) and np.all(qb[self.b!=0]==0) and np.array_equal(nativez[~c],z[~c]),"sediment boundary")
        rem,a=remap(V,W,H,qb*DX*DT,self.t,self.hd,c,self.act,A)
        dzv=(nativez-z)*(1-PHI)*A;dV=a["V_after"]-V
        exlocal=float(np.max(np.abs((dzv-dV)[c])))
        newz=np.where(c,ds.bedrock.values+a["thickness_after"],ds.z0.values)
        req(exlocal<=1e-10 and np.max(np.abs(newz-nativez))<=1e-12,"Exner/geometry")
        gr,ga_a=ga(a,h,1e-6)
        req(rem["GA_domain_valid"] and rem["solid_out_m3"]==0,"pore domain/boundary")
        nxt=ds.copy(deep=True)
        for key,val in (("V",a["V_after"]),("W",ga_a["GA_W_after"]),("H",ga_a["GA_H_after"]),("q",q),("hl",hl)):np.copyto(nxt[key].values,val)
        at=nxt.attrs;at["step"]=k;at["time_s"]=k*DT
        adds={"rain":rr*DT*A*int(c.sum()),"out":out,"inc":inc,"solid_out":rem["solid_out_m3"],"erosion":float(-dV[dV<0].sum()),"deposition":float(dV[dV>0].sum()),"release":rem["released_water_m3"],"infiltration":gr["infiltration_m3"]}
        for name,x in adds.items():at[name]+=x
        # Validate each committed state. No residual repair.
        sr,wr=validate(nxt);self.ds=nxt
        metrics={"step":k,**{name:at[name] for name in LEDGER},"solid_res":sr,"water_res":wr,"flow_local":flowlocal,"exner_local":exlocal,"pore_local":rem["max_local_water_residual_m3"],"ga_local":gr["max_local_water_residual_m3"],"min_F":rem["min_F_equiv_m"],"max_dz":float(np.max(np.abs(newz-ds.z0.values))),"owner_sha":digest(nxt)}
        trace={"V":nxt.V.values.copy(),"W":nxt.W.values.copy(),"H":nxt.H.values.copy(),"q":q,"hl":hl,"flow_h":h,"qb":qb,"z":newz}
        return metrics,trace
def save_series(rows,traces,p):
    np.savez_compressed(p,**{k:np.asarray([v[k] for v in rows]) for k in rows[0]},**{"trace_"+k:np.asarray([v[k] for v in traces]) for k in traces[0]})
def producer(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir()
    start=time.monotonic();RBD=load_provider(out);rep={"experiment":"P3.3","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"];finals={}
    for case in CASES:
        folder=out/case;folder.mkdir();rows=[];traces=[]
        try:
            m=Model(initial(case),RBD);s0=snap(m.ds,folder/"initial.nc")
            for k in range(1,STEPS+1):
                if case=="active" and k==MID+1:
                    d=digest(m.ds);a=Model(m.ds,RBD).flow(k,False);b=Model(m.ds,RBD).flow(k,True)
                    branch={"q_diff":float(np.max(np.abs(a[1]-b[1]))),"h_diff":float(np.max(np.abs(a[0]-b[0]))),"input_owner_sha":d}
                    req(digest(m.ds)==d,"fork mutation")
                    np.savez_compressed(folder/"fork.npz",real_q=a[1],frozen_q=b[1],real_h=a[0],frozen_h=b[0])
                row,tr=m.advance(k);rows.append(row);traces.append(tr)
                if case=="active" and k==MID+1:
                    branch["real_exact"]=bool(np.array_equal(a[1],tr["q"]) and np.array_equal(a[0],tr["flow_h"]))
                    rep["fork"]=branch
                if k==MID:
                    sm=snap(m.ds,folder/"mid.nc");observe(m.ds)
                    if case=="active":
                        meta=snap(m.ds,cp/"owner.nc");req(meta["owner_sha"]==sm["owner_sha"],"checkpoint")
                        (cp/"manifest.json").write_text(json.dumps({"step":MID,"case":case,**meta},sort_keys=True)+"\n")
            sf=snap(m.ds,folder/"final.nc");observe(m.ds);save_series(rows,traces,folder/"series.npz")
            finals[case]=m.ds.copy(deep=True)
            rep["cases"][case]={"initial":s0,"mid":sm,"final":sf,"steps":len(rows),"last":rows[-1],"max_residuals":{key:max(abs(r[key]) for r in rows) for key in ("solid_res","water_res","flow_local","exner_local","pore_local","ga_local")},"min_F":min(r["min_F"] for r in rows)}
        except Exception as e:
            traceback.print_exc();fail.append(case+":"+repr(e))
            if rows:save_series(rows,traces,folder/"partial.npz")
    def gate(name,v):
        if not bool(v):fail.append(name)
    if len(finals)==4:
        a=rep["cases"]["active"]["last"]
        gate("active_solid",a["erosion"]>.001 and a["deposition"]>.001 and a["max_dz"]>.0001)
        gate("active_water",a["release"]>1e-5 and a["infiltration"]>1e-5 and a["out"]>.01)
        gate("rain",all(abs(rep["cases"][n]["last"]["rain"]-(0 if n=="no_rain" else .018))<=1e-10 for n in CASES))
        rep["comparisons"]={}
        for case in ("frozen","no_sediment"):
            delta=float(np.max(np.abs(finals["active"].q.values-finals[case].q.values)));rep["comparisons"][case+"_q_diff"]=delta;gate(case+"_q_diff",delta>1e-7)
        delta=float(np.max(np.abs(finals["active"].H.values-finals["no_rain"].H.values))/A);rep["comparisons"]["rain_h_diff"]=delta;gate("rain_h_diff",delta>1e-6)
        f=rep["fork"];gate("geometry_isolated",f["q_diff"]>1e-7 and f["h_diff"]>1e-8 and f["real_exact"])
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_RECIPROCAL_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
def consumer(checkpoint,dest):
    cp=pathlib.Path(checkpoint);out=pathlib.Path(dest);out.mkdir();meta=json.loads((cp/"manifest.json").read_text())
    req(meta["step"]==MID and meta["case"]=="active","checkpoint metadata");ds=read(cp/"owner.nc",meta,MID)
    controls={}
    def reject(name,fn):
        try:fn()
        except ValueError:controls[name]=True
        else:raise ValueError("accepted "+name)
    bad=out/"bad.nc";b=bytearray((cp/"owner.nc").read_bytes());b[-1]^=1;bad.write_bytes(b);reject("file",lambda:read(bad,meta,MID));bad.unlink()
    def verify(x):
        validate(x,MID);req(digest(x)==meta["owner_sha"],"logical digest")
    alt=ds.copy(deep=True);alt.q.values[np.flatnonzero(alt.active.values)[0]]+=1e-4;reject("q",lambda:verify(alt))
    alt=ds.copy(deep=True);alt["duplicate"]=alt.H.copy();reject("duplicate",lambda:verify(alt))
    alt=ds.copy(deep=True);i=np.flatnonzero(alt.core.values)[0];alt.W.values[i]=geometry(alt.V.values,A)[1][i]+.01;reject("capacity",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.attrs["step"]+=1;reject("clock",lambda:verify(alt))
    RBD=load_provider(out);m=Model(ds,RBD);req(digest(m.ds)==meta["owner_sha"],"restored")
    rows=[];traces=[]
    for k in range(MID+1,STEPS+1):
        a,b=m.advance(k);rows.append(a);traces.append(b)
    fin=snap(m.ds,out/"final.nc");save_series(rows,traces,out/"series.npz")
    rep={"steps":len(rows),"first_step":rows[0]["step"],"input_owner_sha":meta["owner_sha"],"final":fin,"controls":controls,"outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_RECIPROCAL_RESTART_JSON "+json.dumps(rep,sort_keys=True),flush=True)
def compare(reference,consumer_path,dest):
    p=pathlib.Path(reference);c=pathlib.Path(consumer_path);a=json.loads((p/"result.json").read_text());b=json.loads((c/"result.json").read_text())
    req(not a["failed_preregistered_gates"] and b["steps"]==MID and b["first_step"]==MID+1 and len(b["controls"])==5 and all(b["controls"].values()),"source or controls")
    req(a["cases"]["active"]["mid"]["owner_sha"]==b["input_owner_sha"] and a["cases"]["active"]["final"]["owner_sha"]==b["final"]["owner_sha"],"owner")
    with np.load(p/"active/series.npz",allow_pickle=False) as x,np.load(c/"series.npz",allow_pickle=False) as y:
        req(set(x.files)==set(y.files) and all(np.array_equal(x[k][MID:],y[k]) for k in x.files),"suffix arrays")
    req(digest(read(p/"active/final.nc"))==digest(read(c/"final.nc")),"final state")
    rep={"suffix_exact":True,"steps":MID,"controls":b["controls"],"owner_sha":b["final"]["owner_sha"]}
    pathlib.Path(dest).write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_RECIPROCAL_COMPARE_JSON "+json.dumps(rep,sort_keys=True),flush=True)
if __name__=="__main__":{"produce":producer,"consume":consumer,"compare":compare}[sys.argv[1]](*sys.argv[2:])
