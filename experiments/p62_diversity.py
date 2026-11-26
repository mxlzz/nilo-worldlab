"""P6.2A: seeded mass redistribution; original P6.1-R operators unchanged."""
import hashlib,io,json,pathlib,resource,sys,time,traceback,zipfile
import numpy as np
import xarray as xr
from scipy.ndimage import gaussian_filter
import p6_horizon_r as b
from p3_exner import req
from p3_water import digest,sha
MODES=("control","seed26061101","seed26061102")
SEEDS={"control":0,"seed26061101":26061101,"seed26061102":26061102}
STEPS=480; MID=240
old_initial=b.initial
old_rate=b.rate;old_sun=b.sun;old_total=b.forcing_total
def initial(mode):
    ds=old_initial("periodic")
    seed=SEEDS[mode];c=ds.core.values.astype(bool)
    if seed:
        g=gaussian_filter(np.random.Generator(np.random.PCG64(seed)).standard_normal((5,5)),1.0,mode="reflect").ravel()
        f=b.frac_initial(ds)[c]
        for group in (.25,.75):
            mask=f==group;g[mask]-=g[mask].mean()
        g*=.02/np.max(np.abs(g))
        factor=(.2+g)/.2
        M=ds.V_class.values.copy();M[c]*=factor[:,None]
        ds.V_class.values[:]=M
        z=ds.z0.values.copy()
        z[c]=ds.bedrock.values[c]+M[c].sum(axis=1)/((1-b.PHI)*b.A)
        ds.z0.values[:]=z
        E=b.capacity(M.sum(axis=1),ds.W.values,ds.H.values)*(300.-b.TR)
        ds.E_sensible_J.values[:]=E
        ds.attrs["thermal_initial_J"]=float(E.sum())
    ds.attrs.update(thermal_schema="P6.2A-v1",thermal_mode=mode,terrain_seed=seed,terrain_amplitude_m=.02 if seed else 0.)
    return ds
def configure():
    b.initial=initial;b.MODES=MODES;b.STEPS=STEPS;b.MID=MID
    b.RESTART={"seed26061102":MID}
    b.rate=lambda mode,k:old_rate("periodic",k)
    b.sun=lambda mode,k:old_sun("periodic",k)
    b.forcing_total=lambda mode,k:old_total("periodic",k)
def metrics(ds):
    c=ds.core.values.astype(bool);z=b.zphys(ds);dz=z-ds.z0.values
    return {"roughness_m":float(np.std(z[c])),"rms_dz_m":float(np.sqrt(np.mean(dz[c]**2))),"max_dz_m":float(np.max(np.abs(dz[c]))),"fine_change":float(np.max(np.abs(b.fraction(ds)-b.frac_initial(ds))))}
def archive_control():
    p=pathlib.Path("evidence/P6/P61R/p6-horizon-r.zip")
    req(sha(p)=="a1309b755e908093755f3f29bbc83a5ab9820f741195c1eb2f88c6efa8e3adf9","P61R ZIP")
    with zipfile.ZipFile(p) as z:
        name="reference/periodic/step0480.nc";data=z.read(name)
        req(hashlib.sha256(data).hexdigest()=="8e028fe40b29e0eaf9b350b52549bad8d33153954a19be66d46c619612fd4aee","P61R member")
    with xr.open_dataset(io.BytesIO(data),engine="scipy") as d:return d.load()
def produce(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir()
    t=time.monotonic();RBD=b.load_provider(out);rep={"experiment":"P6.2A","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"]
    start={};end={};ref=old_initial("periodic")
    for mode in MODES:
        folder=out/mode;folder.mkdir();rows=[];snaps={};m=None;k=0
        try:
            ds=initial(mode);c=ds.core.values.astype(bool)
            req(digest(ds)==digest(initial(mode)),"generation exact")
            req(np.array_equal(ds.bedrock.values,ref.bedrock.values) and np.array_equal(ds.z0.values[~c],ref.z0.values[~c]),"bedrock/halo")
            req(np.max(np.abs(ds.V_class.values.sum(axis=0)-ref.V_class.values.sum(axis=0)))<=1e-12,"initial class stock")
            req(np.array_equal(ds.W.values,ref.W.values) and np.array_equal(ds.H.values,ref.H.values),"initial water")
            req(abs(ds.attrs["thermal_initial_J"]-ref.attrs["thermal_initial_J"])<=1e-5,"initial energy")
            thick=b.total(ds)[c]/((1-b.PHI)*b.A)
            req(np.all((thick>=.18-1e-12)&(thick<=.22+1e-12)),"initial thickness")
            req(np.max(np.abs(b.zphys(ds)-ds.z0.values))<=1e-12,"initial relief derivation")
            if mode=="control":req(all(np.array_equal(ds[n].values,ref[n].values) for n in ref.variables),"control initial")
            start[mode]=ds.copy(deep=True);m=b.Model(ds,RBD)
            snaps["0"]=b.snap(m.ds,folder/"step0000.nc");b.observe(m.ds)
            for k in range(1,STEPS+1):
                row,_=m.advance(k);rows.append(row)
                if k in (MID,STEPS):
                    snaps[str(k)]=b.snap(m.ds,folder/f"step{k:04d}.nc");b.observe(m.ds)
                    if mode=="seed26061102" and k==MID:
                        d=cp/mode;d.mkdir();meta=b.snap(m.ds,d/"owner.nc")
                        (d/"manifest.json").write_text(json.dumps({"step":k,"mode":mode,**meta},sort_keys=True)+"\n")
            end[mode]=m.ds.copy(deep=True);mm=metrics(m.ds)
            if mode!="control":req(mm["max_dz_m"]>1e-4 and mm["fine_change"]>1e-5 and m.ds.attrs["vapor_kg"]>0 and m.ds.attrs["out"]>0,"active evolution")
        except Exception as e:
            traceback.print_exc();fail.append(mode+":"+repr(e))
            if m is not None:
                snaps["last_valid"]=b.snap(m.ds,folder/"last_valid.nc")
                rep.setdefault("stops",{})[mode]={"attempted_step":k,"last_valid_step":int(m.ds.attrs["step"]),"error":repr(e)}
        b.series(rows,folder/"series.npz")
        rep["cases"][mode]={"steps":len(rows),"snapshots":snaps,"initial":metrics(start[mode]) if mode in start else None,"final":metrics(end[mode]) if mode in end else None,"last":rows[-1] if rows else None,"max_residuals":{n:max(abs(v[n]) for v in rows) for n in ("class_res","solid_res","water_res","energy_res","max_energy_op")} if rows else {},"limited_steps":sum(v["limited_nodes"]>0 for v in rows)}
    if len(end)==3:
        ref480=archive_control();tol={"V_class":1e-10,"W":1e-10,"H":1e-10,"q":1e-9,"hl":1e-10,"E_sensible_J":.001}
        diffs={n:float(np.max(np.abs(end["control"][n].values-ref480[n].values))) for n in b.DYNAMIC}
        rep["control_diffs"]=diffs
        if not all(diffs[n]<=tol[n] for n in diffs):fail.append("P61R control")
        pairs={}
        for i,u in enumerate(MODES):
            for v in MODES[i+1:]:
                z0=float(np.max(np.abs(start[u].z0.values-start[v].z0.values)))
                du=b.zphys(end[u])-start[u].z0.values;dv=b.zphys(end[v])-start[v].z0.values
                d=float(np.max(np.abs(du-dv)))
                pairs[u+":"+v]={"initial_z_max_diff_m":z0,"response_max_diff_m":d}
                if not(z0>.001 and d>1e-5):fail.append("diversity "+u+":"+v)
        rep["pairs"]=pairs
    else:fail.append("incomplete")
    rep["wall_seconds"]=time.monotonic()-t;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_DIVERSITY_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":
    configure()
    {"produce":produce,"consume":b.consumer,"compare":b.compare}[sys.argv[1]](*sys.argv[2:])
