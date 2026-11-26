"""P6.2B: larger coherent terrain using unchanged P6.1-R operators."""
import json,pathlib,resource,sys,time,traceback
import numpy as np
from scipy.ndimage import gaussian_filter
from landlab import RasterModelGrid,NodeStatus
import p6_horizon_r as b
import p3_reciprocal as hydro
import p3_pore as pore
from p3_exner import req
from p3_water import digest,sha
MODES=("control","seed26061103","seed26061104")
SEEDS={"control":0,"seed26061103":26061103,"seed26061104":26061104}
STEPS=240;MID=120;CORE=225
old_initial=b.initial
old_total=b.forcing_total
def grid(flat=False):
    g=RasterModelGrid((17,17),xy_spacing=b.DX)
    g.set_closed_boundaries_at_grid_edges(True,True,True,True)
    g.status_at_node[8]=NodeStatus.FIXED_VALUE
    z=np.ones(g.number_of_nodes) if flat else 1.+.025*g.node_y+.003*g.node_x
    g.at_node["topographic__elevation"]=z.copy()
    c=np.asarray(g.status_at_node)==NodeStatus.CORE
    t=np.asarray(g.node_at_link_tail);h=np.asarray(g.node_at_link_head)
    act=np.zeros(g.number_of_links,dtype=bool);act[g.active_links]=True
    bs=c[t].astype(int)-c[h].astype(int);bs[~act]=0
    req(c.sum()==CORE and np.count_nonzero(bs)==1 and bs[bs!=0][0]==-1,"topology")
    return g,c,act,t,h,bs
def fractions(ds):
    return np.where(ds.core.values.astype(bool),np.where(ds.x.values<=16,.25,.75),.5)
def initial(mode):
    ds=old_initial("periodic");seed=SEEDS[mode];c=ds.core.values.astype(bool)
    if seed:
        noise=gaussian_filter(np.random.Generator(np.random.PCG64(seed)).standard_normal((15,15)),2.,mode="reflect").ravel()
        noise-=noise.mean();noise*=.04/np.max(np.abs(noise))
        ds.z0.values[c]+=noise
        ds.bedrock.values[c]=ds.z0.values[c]-.2
    ds.attrs.update(thermal_schema="P6.2B-v1",thermal_mode=mode,terrain_seed=seed,terrain_amplitude_m=.04 if seed else 0.,terrain_nx=17,terrain_ny=17)
    return ds
def configure():
    hydro.grid=grid;pore.grid=grid;b.grid=grid;b.frac_initial=fractions
    b.initial=initial;b.MODES=MODES;b.STEPS=STEPS;b.MID=MID;b.RESTART={"seed26061104":MID}
    b.forcing_total=lambda mode,k:tuple(x*9 for x in old_total("periodic",k))
def top_share(a):
    a=np.asarray(a);n=max(1,int(np.ceil(len(a)*.2)))
    return float(np.sort(a)[-n:].sum()/a.sum()) if a.sum()>0 else 0.
def measures(ds):
    c=ds.core.values.astype(bool);z=b.zphys(ds);dz=z-ds.z0.values
    return {"elevation_std_m":float(np.std(z[c])),"rms_dz_m":float(np.sqrt(np.mean(dz[c]**2))),"max_dz_m":float(np.max(np.abs(dz[c]))),"fine_change":float(np.max(np.abs(b.fraction(ds)-fractions(ds))))}
def produce(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir()
    start=time.monotonic();RBD=b.load_provider(out);rep={"experiment":"P6.2B","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"]
    ref=initial("control");initials={};finals={}
    for mode in MODES:
        folder=out/mode;folder.mkdir();rows=[];snaps={};m=None;k=0
        try:
            t0=time.monotonic();ds=initial(mode);c=ds.core.values.astype(bool)
            req(digest(ds)==digest(initial(mode)),"generation exact")
            req(np.array_equal(ds.V_class.values,ref.V_class.values) and np.array_equal(ds.W.values,ref.W.values) and np.array_equal(ds.H.values,ref.H.values) and abs(ds.attrs["thermal_initial_J"]-ref.attrs["thermal_initial_J"])<=1e-5,"same initial inventory")
            req(np.array_equal(ds.z0.values[~c],ref.z0.values[~c]) and np.max(np.abs(b.zphys(ds)-ds.z0.values))<=1e-12,"halo/initial derived relief")
            initials[mode]=ds.copy(deep=True);m=b.Model(ds,RBD)
            flow=np.zeros(len(ds.link));erosion=np.zeros(len(c));dep=erosion.copy()
            inside=c[m.t]&c[m.hd]&m.act
            snaps["0"]=b.snap(m.ds,folder/"step0000.nc");b.observe(m.ds)
            for k in range(1,STEPS+1):
                old=b.zphys(m.ds)
                if mode=="seed26061104" and k==MID+1:
                    before=digest(m.ds)
                    a=b.Model(m.ds,RBD).flow(k,False);f=b.Model(m.ds,RBD).flow(k,True)
                    req(digest(m.ds)==before,"fork mutation")
                    qdiff=float(np.max(np.abs(a[1]-f[1])));hdiff=float(np.max(np.abs(a[0]-f[0])))
                    np.savez_compressed(folder/"fork.npz",q_live=a[1],q_frozen=f[1],h_live=a[0],h_frozen=f[0])
                row,tr=m.advance(k);rows.append(row)
                if mode=="seed26061104" and k==MID+1:
                    exact=bool(np.array_equal(a[1],tr["q"]) and np.array_equal(a[0],tr["flow_h"]))
                    rep["fork"]={"q_diff":qdiff,"h_diff":hdiff,"real_exact":exact,"owner_sha":before}
                    req(qdiff>1e-7 and hdiff>1e-8 and exact,"live relief causal fork")
                flow+=np.abs(tr["q"])*b.DX*b.DT
                delta=(tr["z"]-old)*(1-b.PHI)*b.A
                erosion+=np.maximum(-delta,0);dep+=np.maximum(delta,0)
                if k in (MID,STEPS):
                    snaps[str(k)]=b.snap(m.ds,folder/f"step{k:04d}.nc");b.observe(m.ds)
                    if mode=="seed26061104" and k==MID:
                        d=cp/mode;d.mkdir();meta=b.snap(m.ds,d/"owner.nc")
                        (d/"manifest.json").write_text(json.dumps({"step":k,"mode":mode,**meta},sort_keys=True)+"\n")
            mm=measures(m.ds);finals[mode]=m.ds.copy(deep=True)
            if mode!="control":req(mm["max_dz_m"]>1e-4 and mm["fine_change"]>1e-5 and m.ds.attrs["vapor_kg"]>0 and m.ds.attrs["out"]>0,"evolution")
            req(abs(erosion.sum()-m.ds.attrs["erosion"])<=1e-9 and abs(dep.sum()-m.ds.attrs["deposition"])<=1e-9,"spatial diagnostic ledger")
            np.savez_compressed(folder/"maps.npz",z0=ds.z0.values,z=b.zphys(m.ds),core=c,tail=m.t,head=m.hd,interior_links=inside,cumulative_abs_water_m3=flow,erosion_m3=erosion,deposition_m3=dep)
            rep["cases"][mode]={"initial":measures(ds),"final":mm,"interior_flow_top20_share":top_share(flow[inside]),"erosion_top20_share":top_share(erosion[c]),"wall_seconds":time.monotonic()-t0}
        except Exception as e:
            traceback.print_exc();fail.append(mode+":"+repr(e))
            if m is not None:
                snaps["last_valid"]=b.snap(m.ds,folder/"last_valid.nc")
                rep.setdefault("stops",{})[mode]={"attempted_step":k,"last_valid_step":int(m.ds.attrs["step"]),"error":repr(e)}
        b.series(rows,folder/"series.npz")
        rep["cases"].setdefault(mode,{}).update(steps=len(rows),snapshots=snaps,last=rows[-1] if rows else None,max_residuals={n:max(abs(v[n]) for v in rows) for n in ("class_res","solid_res","water_res","energy_res","max_energy_op")} if rows else {},limited_steps=sum(v["limited_nodes"]>0 for v in rows))
    if len(finals)==3:
        pairs={}
        for i,u in enumerate(MODES):
            for v in MODES[i+1:]:
                a=initials[u];d=initials[v]
                z0=float(np.max(np.abs(a.z0.values-d.z0.values)))
                du=b.zphys(finals[u])-a.z0.values;dv=b.zphys(finals[v])-d.z0.values
                resp=float(np.max(np.abs(du-dv)))
                pairs[u+":"+v]={"initial_z_max_diff_m":z0,"response_max_diff_m":resp}
                if not(z0>.001 and resp>1e-5):fail.append("diversity "+u+":"+v)
        rep["pairs"]=pairs
    else:fail.append("incomplete")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_TERRITORY_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
def compare(reference,consumer,dest):
    b.compare(reference,consumer,dest)
    a=b.read(pathlib.Path(reference)/"seed26061104/step0240.nc")
    c=b.read(pathlib.Path(consumer)/"seed26061104/step0240.nc")
    req(digest(a)==digest(c),"final snapshot readback")
if __name__=="__main__":
    configure()
    {"produce":produce,"consume":b.consumer,"compare":compare}[sys.argv[1]](*sys.argv[2:])
