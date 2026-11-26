"""P7.1: functional concentration on erosion-created zones vs fixed hydraulic geometry."""
import json,pathlib,resource,sys,time,traceback
import numpy as np
import p62_territory as t
import p6_horizon_r as b
from p3_exner import req
from p3_water import digest,sha
TERRAINS={"control":0,"seed26061105":26061105,"seed26061106":26061106}
MODES=tuple(a+"_"+v for a in TERRAINS for v in ("live","frozen"))
STEPS=480;MID=240;RESERVED="seed26061106_live"
Base=b.Model
class Model(Base):
    def flow(self,k,freeze=False):
        return super().flow(k,freeze or self.ds.attrs["thermal_mode"].endswith("_frozen"))
def initial(mode):
    ds=t.initial(mode)
    ds.attrs["thermal_schema"]="P7.1-v1"
    return ds
def configure():
    t.MODES=MODES;t.SEEDS={m:TERRAINS[m.rsplit("_",1)[0]] for m in MODES}
    t.configure()
    b.initial=initial;b.Model=Model;b.MODES=MODES;b.STEPS=STEPS;b.MID=MID;b.RESTART={RESERVED:MID}
def region(ds):
    c=ds.core.values.astype(bool)
    return c&(ds.x.values>=6)&(ds.x.values<=26)&(ds.y.values>=6)&(ds.y.values<=26)
def select(ds,z):
    ids=np.flatnonzero(region(ds));lower=ds.z0.values-z
    req(len(ids)==121,"region")
    chosen=ids[np.lexsort((ids,-lower[ids]))[:25]]
    mask=np.zeros(len(lower),bool);mask[chosen]=True
    return mask
def share(a,mask,roi):
    den=float(a[roi].sum());req(den>0 and np.isfinite(a).all(),"nonzero functional flow")
    return float(a[mask].sum()/den)
def producer(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir()
    start=time.monotonic();RBD=b.load_provider(out);rep={"experiment":"P7.1","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"];maps={}
    for terrain in TERRAINS:
        a=initial(terrain+"_live");d=initial(terrain+"_frozen")
        req(all(np.array_equal(a[n].values,d[n].values) for n in a.variables),"paired initial arrays")
        req(all(a.attrs[k]==d.attrs[k] for k in a.attrs if k!="thermal_mode"),"paired initial metadata")
    for mode in MODES:
        folder=out/mode;folder.mkdir();rows=[];snaps={};m=None;k=0
        try:
            t0=time.monotonic();ds=initial(mode);req(digest(ds)==digest(initial(mode)),"initial deterministic")
            m=Model(ds,RBD);c=ds.core.values.astype(bool);roi=region(ds)
            inside=c[m.t]&c[m.hd]&m.act;tail=m.t[inside];head=m.hd[inside]
            ingress=np.zeros((3,len(c)));zmid=None
            snaps["0"]=b.snap(m.ds,folder/"step0000.nc");b.observe(m.ds)
            for k in range(1,STEPS+1):
                row,tr=m.advance(k);rows.append(row)
                q=tr["q"][inside];destnode=np.where(q>=0,head,tail)
                window=0 if k<=240 else 1 if k<=360 else 2
                np.add.at(ingress[window],destnode,np.abs(q)*b.DX*b.DT)
                if k in (MID,STEPS):
                    snaps[str(k)]=b.snap(m.ds,folder/f"step{k:04d}.nc");b.observe(m.ds)
                    if k==MID:
                        zmid=b.zphys(m.ds).copy()
                        if mode==RESERVED:
                            d=cp/mode;d.mkdir();meta=b.snap(m.ds,d/"owner.nc")
                            (d/"manifest.json").write_text(json.dumps({"step":k,"mode":mode,**meta},sort_keys=True)+"\n")
            zfinal=b.zphys(m.ds);mask=select(ds,zmid)
            req(roi.sum()==121 and mask.sum()==25 and np.all(ingress>=0) and np.isfinite(ingress).all(),"diagnostic bounds")
            maps[mode]=(ds.z0.values.copy(),zmid,zfinal,roi,mask,ingress)
            np.savez_compressed(folder/"maps.npz",z0=ds.z0.values,z240=zmid,z480=zfinal,roi=roi,own_mask240=mask,ingress_m3=ingress)
            rep["cases"][mode]={"wall_seconds":time.monotonic()-t0,"final":t.measures(m.ds),"own_mask_shares":[share(i,mask,roi) for i in ingress]}
        except Exception as e:
            traceback.print_exc();fail.append(mode+":"+repr(e))
            if m is not None:
                snaps["last_valid"]=b.snap(m.ds,folder/"last_valid.nc")
                rep.setdefault("stops",{})[mode]={"attempted_step":k,"last_valid_step":int(m.ds.attrs["step"]),"error":repr(e)}
        b.series(rows,folder/"series.npz")
        rep["cases"].setdefault(mode,{}).update(steps=len(rows),snapshots=snaps,last=rows[-1] if rows else None,max_residuals={n:max(abs(v[n]) for v in rows) for n in ("class_res","solid_res","water_res","energy_res","max_energy_op")} if rows else {},limited_steps=sum(v["limited_nodes"]>0 for v in rows))
    rep["pairs"]={}
    if len(maps)==len(MODES) and not fail:
        try:
            for terrain in TERRAINS:
                a=maps[terrain+"_live"];f=maps[terrain+"_frozen"];z0,zm,zf,roi,mask,ing=a
                req(np.array_equal(a[0],f[0]) and np.array_equal(roi,f[3]),"paired maps")
                live=[share(i,mask,roi) for i in ing];frozen=[share(i,mask,roi) for i in f[5]]
                diff=[x-y for x,y in zip(live,frozen)]
                nmid=int(np.count_nonzero((z0-zm>=.001)&mask));nfinal=int(np.count_nonzero((z0-zf>=.001)&mask))
                structural=nmid>=13 and nfinal>=13;functional=diff[1]>=.02 and diff[2]>=.02
                rep["pairs"][terrain]={"selected_nodes":np.flatnonzero(mask).tolist(),"n_lowered_240":nmid,"n_lowered_480":nfinal,"structural_gate":structural,"live_share":live,"frozen_share":frozen,"share_difference":diff,"functional_gate":functional,"joint_gate":structural and functional}
            count=sum(rep["pairs"][s]["joint_gate"] for s in TERRAINS if s!="control")
            rep["functional_pass_count"]=count;rep["scientific_outcome"]={0:"negative_limited",1:"inconclusive",2:"positive_limited"}[count]
        except Exception as e:
            traceback.print_exc();fail.append("diagnostics:"+repr(e))
    else:fail.append("incomplete")
    if fail:rep["scientific_outcome"]="invalid_or_incomplete"
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_FUNCTION_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":
    configure()
    {"produce":producer,"consume":b.consumer,"compare":b.compare}[sys.argv[1]](*sys.argv[2:])
