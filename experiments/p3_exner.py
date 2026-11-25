"""P3.2a: isolated upstream RBD; native Exner and measured link mass.
No modified solver, no bidirectional coupling, no general substrate model.
"""
import hashlib,json,pathlib,resource,sys,time,traceback
from importlib.metadata import version
from urllib.request import urlopen
import numpy as np
from landlab import RasterModelGrid,NodeStatus
from landlab.components import OverlandFlow

DX=2.;DT=.1;PHI=.4;RHO_S=2650.;FILM=1e-8
SHAPE=(7,7);N=49
MANIFEST=pathlib.Path("research/P3/RBD_SOURCE_MANIFEST.json")
def req(ok,message):
    if not bool(ok):raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load_provider(dest):
    manifest=json.loads(MANIFEST.read_text())
    root=pathlib.Path("rbd-upstream");pkg=root/"river_bed_dynamics";pkg.mkdir(parents=True)
    base="https://raw.githubusercontent.com/angelmons/landlab/"+manifest["commit"]+"/"
    for path,meta in manifest["files"].items():
        with urlopen(base+path,timeout=30) as response:
            b=response.read(meta["bytes"]+1)
        req(len(b)==meta["bytes"],"source size "+path)
        req(hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()==meta["git_blob_sha1"],"source digest "+path)
        target=root/"LICENSE.md" if path=="LICENSE.md" else pkg/pathlib.Path(path).name
        target.write_bytes(b)
    (dest/"upstream_license.md").write_bytes((root/"LICENSE.md").read_bytes())
    (dest/"source_manifest.json").write_bytes(MANIFEST.read_bytes())
    sys.path.insert(0,str(root.resolve()))
    from river_bed_dynamics import RiverBedDynamics
    req(version("landlab")=="2.11.0","landlab version")
    return RiverBedDynamics
def grid(flat=False):
    g=RasterModelGrid(SHAPE,xy_spacing=DX)
    g.set_closed_boundaries_at_grid_edges(True,True,True,True)
    g.status_at_node[3]=NodeStatus.FIXED_VALUE
    z=np.ones(N) if flat else 1.+.025*g.node_y+.003*g.node_x+.04*np.exp(-((g.node_x-6.)**2+(g.node_y-6.)**2)/4.)
    g.at_node["topographic__elevation"]=z.copy()
    core=np.asarray(g.status_at_node)==NodeStatus.CORE
    t=np.asarray(g.node_at_link_tail);h=np.asarray(g.node_at_link_head)
    active=np.zeros(g.number_of_links,dtype=bool);active[g.active_links]=True
    bsign=core[t].astype(int)-core[h].astype(int);bsign[~active]=0
    req(int(core.sum())==25 and np.count_nonzero(bsign)==1 and bsign[bsign!=0][0]==-1,"boundary topology")
    return g,core,active,t,h,bsign
def one(RBD,name,h=None,hlink=None,u=None,prev=None):
    g,core,active,t,head,bsign=grid(name=="flat")
    n=g.number_of_links
    if h is None:
        h=np.full(N,FILM if name=="dry" else .04)
        hlink=np.where(active,FILM if name=="dry" else .04,0.)
        u=np.zeros(n);prev=np.zeros(n)
    h=np.asarray(h).copy();hlink=np.asarray(hlink).copy()
    u=np.asarray(u).copy();prev=np.asarray(prev).copy()
    req(all(np.isfinite(a).all() for a in (h,hlink,u,prev)) and np.all(h>=0) and np.all(hlink>=0),"hydraulic input")
    z0=g.at_node["topographic__elevation"].copy()
    thick=1e-8 if name=="depleted" else .2
    bedrock=z0-thick;stock0=np.where(core,thick*(1-PHI)*DX*DX,0.)
    owner0=stock0.copy()
    g.at_node["surface_water__depth"]=h.copy()
    g.at_link["surface_water__depth"]=hlink.copy()
    g.at_link["surface_water__velocity"]=u.copy()
    fixed=np.zeros(n)
    if name=="imposed_export":fixed[bsign!=0]=-1e-4
    r=RBD(g,gsd=[[3,100],[2,50],[1,0]],rho=1000,rho_s=RHO_S,bedload_equation="MPM",lambda_p=PHI,variable_critical_shear_stress=False,outlet_boundary_condition="fixedValue",dt=DT,sed_transp__bedload_rate_fix_link=fixed,surface_water__velocity_prev_time_link=prev.copy(),track_stratigraphy=False)
    req(np.array_equal(r._surface_water__velocity_prev_time_link,prev),"previous velocity identity")
    r.run_one_step()
    z=g.at_node["topographic__elevation"].copy()
    qb=np.asarray(r._sed_transp__bedload_rate_link).copy()
    net=np.asarray(r._sed_transp__net_bedload_node).copy()
    shear=np.asarray(r._surface_water__shear_stress_link).copy()
    req(all(np.isfinite(a).all() for a in (z,qb,net,shear)),"nonfinite native sediment")
    req(np.all(qb[~active]==0) and np.array_equal(z[~core],z0[~core]),"inactive or halo")
    for x,y in ((g.at_node["surface_water__depth"],h),(g.at_link["surface_water__depth"],hlink),(g.at_link["surface_water__velocity"],u),(r._surface_water__velocity_prev_time_link,prev)):
        req(np.array_equal(x,y),"hydraulic mutation")
    # Independent, signed, finite-volume balance from actual link flux.
    vol=qb*DX*DT
    delta=np.zeros(N);np.add.at(delta,t,-vol);np.add.at(delta,head,vol)
    dz_vol=(z-z0)*(1-PHI)*DX*DX
    local=float(np.max(np.abs(dz_vol[core]-delta[core])))
    native_local=float(np.max(np.abs(net[core]*DT+delta[core])))
    b=vol*bsign;out=float(b[b>0].sum());inc=float(-b[b<0].sum())
    residual=float(dz_vol[core].sum()+out-inc)
    req(local<=1e-10 and native_local<=1e-10 and abs(residual)<=1e-9,"Exner balance")
    gross=np.zeros(N);np.add.at(gross,t,np.maximum(vol,0));np.add.at(gross,head,np.maximum(-vol,0))
    candidate=stock0+delta
    accepted=bool(np.all(gross[core]<=stock0[core]) and np.all(candidate[core]>=0) and np.all(z[core]>=bedrock[core]-1e-12))
    owner=stock0.copy()
    if accepted:owner[core]=candidate[core]
    else:req(np.array_equal(owner,owner0),"rejected candidate changed owner")
    erosion=float(-dz_vol[core][dz_vol[core]<0].sum());dep=float(dz_vol[core][dz_vol[core]>0].sum())
    result={"accepted":accepted,"erosion_solid_m3":erosion,"deposition_solid_m3":dep,"out_solid_m3":out,"in_solid_m3":inc,"residual_solid_m3":residual,"max_local_m3":local,"max_native_net_error_m3":native_local,"max_qb_m2s":float(np.max(np.abs(qb))),"max_shear_pa":float(np.max(shear)),"max_dz_m":float(np.max(np.abs(z-z0))),"min_candidate_thickness_m":float(np.min((z-bedrock)[core])),"initial_dry_mass_kg":float(stock0[core].sum()*RHO_S),"accepted_final_dry_mass_kg":float(owner[core].sum()*RHO_S),"outlet_qb_m2s":float(qb[bsign!=0][0]),"median_grain_mm":float(np.median(r._bed_surf__median_size_link[active]))}
    arrays={"z_before":z0,"z_after":z,"bedrock":bedrock,"stock_before":stock0,"stock_candidate":candidate,"stock_accepted":owner,"qb":qb,"delta_solid_volume":delta,"dz_solid_volume":dz_vol,"gross_out":gross,"native_net":net,"shear":shear,"h_node":h,"h_link":hlink,"u":u,"u_prev":prev,"tail":t,"head":head,"core":core,"active":active,"boundary_sign":bsign}
    return result,arrays
def live(RBD,out,report):
    g,core,active,t,head,bsign=grid()
    g.at_node["surface_water__depth"]=np.zeros(N)
    w=OverlandFlow(g,h_init=FILM,alpha=.7,theta=.8,mannings_n=.03,steep_slopes=False,rainfall_intensity=0.)
    h=g.at_node["surface_water__depth"];q=g.at_link["surface_water__discharge"]
    h[:]=np.where(core,.04,FILM);u=np.zeros(g.number_of_links)
    initial=float(h[core].sum()*DX*DX);outcum=incum=0.;rows=[]
    samples={}
    for step in range(1,121):
        old=h.copy();prev=u.copy()
        req(DT<=w.calc_time_step(),"live CFL")
        req(w.run_one_step(DT)==DT,"live duration")
        div=g.calc_flux_div_at_node(q)
        local=float(np.max(np.abs((h-old+div*DT)[core])))
        b=q*bsign*DX*DT;outcum+=float(b[b>0].sum());incum+=float(-b[b<0].sum())
        res=float(h[core].sum()*DX*DX+outcum-incum-initial)
        req(local<=1e-10 and abs(res)<=1e-8 and np.all(h>=0) and np.all(h[~core]==FILM) and np.all(q[~active]==0),"live water balance")
        hl=np.asarray(w._h_links).copy()
        req(np.all(hl[active]>0),"native flow depth")
        u=np.zeros_like(q);u[active]=q[active]/hl[active]
        req(np.isfinite(u).all(),"velocity")
        rows.append((step,outcum,incum,float(h[core].sum()*DX*DX),res,local,float(h[core].min()),float(np.max(np.abs(q)))))
        if step in (40,80,120):
            name="live_"+str(step);result,arrays=one(RBD,name,h,hl,u,prev)
            req(result["accepted"],"live sediment finite stock")
            arrays["water_q"]=q.copy()
            np.savez_compressed(out/(name+".npz"),**arrays)
            samples[name]=result
            report["cases"][name]=result
    np.savez_compressed(out/"live_series.npz",series=np.asarray(rows),columns=np.array(["step","out_m3","in_m3","surface_m3","residual_m3","max_local_m","min_h_m","max_q_m2s"]))
    return {"steps":120,"samples":list(samples),"out_m3":outcum,"in_m3":incum,"max_residual_m3":max(abs(r[4]) for r in rows),"max_local_m":max(r[5] for r in rows),"any_live_erosion_deposition":any(x["erosion_solid_m3"]>1e-8 and x["deposition_solid_m3"]>1e-8 for x in samples.values())}
def main(dest):
    out=pathlib.Path(dest);out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    report={"experiment":"P3.2a","cases":{},"failed_preregistered_gates":[]}
    fail=report["failed_preregistered_gates"]
    try:
        RBD=load_provider(out)
        fixtures={}
        for name in ("active","flat","dry","imposed_export","depleted","active_replay"):
            try:
                result,arrays=one(RBD,"active" if name=="active_replay" else name)
                np.savez_compressed(out/(name+".npz"),**arrays)
                report["cases"][name]=result;fixtures[name]=(result,arrays)
                req(result["accepted"]==(name!="depleted"),"stock acceptance")
                if name in ("active","active_replay"):
                    req(result["erosion_solid_m3"]>1e-6 and result["deposition_solid_m3"]>1e-6,"active erosion/deposition")
                    req(result["out_solid_m3"]==result["in_solid_m3"]==0,"native sediment closed outlet")
                if name in ("flat","dry"):req(result["max_qb_m2s"]==result["max_dz_m"]==0,"zero control")
                if name=="imposed_export":req(abs(result["out_solid_m3"]-.00002)<=1e-12 and result["in_solid_m3"]==0,"imposed export")
            except Exception as e:
                traceback.print_exc();fail.append(name+":"+repr(e))
        if "active" in fixtures and "active_replay" in fixtures:
            a,b=fixtures["active"],fixtures["active_replay"]
            report["replay_exact"]=a[0]==b[0] and all(np.array_equal(a[1][k],b[1][k]) for k in a[1])
            if not report["replay_exact"]:fail.append("replay")
        try:
            report["live"]=live(RBD,out,report)
            if not report["live"]["any_live_erosion_deposition"]:fail.append("no live sediment movement above gate")
        except Exception as e:
            traceback.print_exc();fail.append("live:"+repr(e))
    except Exception as e:
        traceback.print_exc();fail.append("provider:"+repr(e))
    report["wall_seconds"]=time.monotonic()-start;report["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report["outputs"]={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_EXNER_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
