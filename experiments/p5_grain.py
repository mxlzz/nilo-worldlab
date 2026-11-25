"""P5.1 native WC flux proposals from retained P4.2, conservative class owner."""
import hashlib,io,json,pathlib,sys,time,resource,zipfile,traceback
import numpy as np
from p3_exner import grid,req,DX,DT,PHI,RHO_S
from p3_reciprocal import load_provider
from p3_pore import input_hash
from p3_water import digest,sha
import xarray as xr
ARCH=pathlib.Path("evidence/P4/P42/15030052011.zip")
ARCH_SHA="8cbf531b85232ec030ca433a0074706c1520a988c1a49fca762c7506a1dd55d8"
STEPS=(40,80,120)
MODES=("uniform","heterogeneous","mirror")
def move(M,f,p,t,head,c,act):
    args=(M,f,p,t,head,c,act);before=input_hash(*args)
    req(M.shape==(len(c),2) and p.shape==(len(f),2) and f.shape==t.shape==head.shape==act.shape,"shape")
    req(all(np.isfinite(x).all() for x in (M,f,p)) and np.all(M>=0) and np.all(M[~c]==0),"bounds")
    req(np.all(p>=0) and np.all(np.abs(p[f!=0].sum(axis=1)-1)<=1e-12),"fractions")
    req(np.all(f[~act]==0),"inactive")
    source=np.where(f>=0,t,head);target=np.where(f>=0,head,t);mag=np.abs(f)
    req(np.all(mag[~c[source]]==0) and np.all(mag[~c[target]]==0),"solid boundary")
    amount=mag[:,None]*p
    gross=np.zeros_like(M);np.add.at(gross,source,amount)
    req(np.all(gross<=M+1e-12),"class overdraw")
    delta=np.zeros_like(M);np.add.at(delta,source,-amount);np.add.at(delta,target,amount)
    new=M+delta
    req(np.all(new>=0) and np.all(new[~c]==0),"candidate bounds")
    residual=np.max(np.abs(new.sum(axis=0)-M.sum(axis=0)))
    req(residual<=1e-9 and input_hash(*args)==before,"class conservation/mutation")
    return new,delta,amount,gross,float(residual)
def fixture():
    M=np.array([[.6,.4],[.4,.6]]);f=np.array([.2]);p=np.array([[.25,.75]])
    t=np.array([0]);hd=np.array([1]);c=np.array([True,True]);a=np.array([True])
    args=(M,f,p,t,hd,c,a)
    N,d,am,g,res=move(*args)
    req(np.allclose(N,[[.55,.25],[.45,.75]],rtol=0,atol=1e-12),"analytic")
    req(np.array_equal(move(M,f*0,p,t,hd,c,a)[0],M),"zero")
    controls={}
    def reject(name,args,why):
        before=input_hash(*args)
        try:move(*args)
        except ValueError as e:req(str(e)==why,"wrong rejection "+name+":"+str(e));controls[name]=True
        else:raise ValueError("accepted "+name)
        req(input_hash(*args)==before,"rejection mutation")
    reject("absent",(np.array([[1.,0.],[.4,.6]]),f,p,t,hd,c,a),"class overdraw")
    reject("overdraw",(M,f*10,p,t,hd,c,a),"class overdraw")
    reject("fraction",(M,f,p*.5,t,hd,c,a),"fractions")
    reject("nan",(M,f*np.nan,p,t,hd,c,a),"bounds")
    reject("inactive",(M,f,p,t,hd,c,~a),"inactive")
    reject("inflow",(np.array([[0.,0.],[.4,.6]]),f,p,t,hd,np.array([False,True]),a),"solid boundary")
    return {"delta":d.tolist(),"controls":controls,"zero":True,"residual":res}
def retained():
    req(sha(ARCH)==ARCH_SHA,"archive sha")
    with zipfile.ZipFile(ARCH) as z:
        req(z.testzip() is None,"archive CRC")
        r=json.loads(z.read("reference/result.json"));req(not r["failed_preregistered_gates"],"source gates")
        data={}
        for name in ("active/initial.nc","active/series.npz"):
            b=z.read("reference/"+name);req(hashlib.sha256(b).hexdigest()==r["outputs"][name],"source output sha");data[name]=b
        with xr.open_dataset(io.BytesIO(data["active/initial.nc"]),engine="scipy") as ds:ini=ds.load()
        req(digest(ini)==r["cases"]["active"]["initial"]["owner_sha"],"source owner")
        with np.load(io.BytesIO(data["active/series.npz"]),allow_pickle=False) as s:x={k:s[k].copy() for k in s.files}
    req(np.array_equal(x["step"],np.arange(1,241)),"source time")
    return ini,x
def one(RBD,ini,x,k,mode,out):
    g,c,act,t,hd,b=grid()
    req(np.array_equal(t,ini["tail"].values) and np.array_equal(hd,ini["head"].values) and np.array_equal(c,ini.core.values.astype(bool)),"topology")
    j=k-1
    V=x["trace_V"][j-1].copy();z=x["trace_z"][j-1].copy()
    h=x["trace_flow_h"][j].copy();q=x["trace_q"][j].copy();hl=x["trace_hl"][j].copy()
    u=np.zeros_like(q);u[act]=q[act]/hl[act];prev=np.zeros_like(q);prev[act]=x["trace_q"][j-1][act]/x["trace_hl"][j-1][act]
    fine=np.where(c,np.where(g.node_x<=6,.25,.75),.5)
    if mode=="uniform":fine=np.full(len(c),.5)
    if mode=="mirror":fine=1.-fine
    frac=np.column_stack((1.-fine,fine));M=V[:,None]*frac
    before=input_hash(V,z,h,q,hl,u,prev,M,fine)
    g.at_node["topographic__elevation"]=z.copy()
    g.at_node["surface_water__depth"]=h.copy();g.at_link["surface_water__depth"]=hl.copy();g.at_link["surface_water__velocity"]=u.copy()
    gsd=np.vstack((np.r_[4.,np.full(len(c),100.)],np.r_[2.,100.*fine],np.r_[1.,np.zeros(len(c))]))
    r=RBD(g,gsd=gsd,bed_surf__gsd_loc_node=np.arange(len(c)),rho=1000,rho_s=RHO_S,bedload_equation="WilcockAndCrowe",lambda_p=PHI,variable_critical_shear_stress=False,outlet_boundary_condition="fixedValue",dt=DT,surface_water__velocity_prev_time_link=prev.copy(),track_stratigraphy=False)
    req(np.allclose(r._gs,[np.sqrt(8),np.sqrt(2)],rtol=0,atol=1e-12),"size mapping")
    req(np.allclose(r._bed_surf__gsd_node,frac,rtol=0,atol=1e-12),"node fraction mapping")
    req(np.allclose(r._bed_surf__gsd_link,.5*(frac[t]+frac[hd]),rtol=0,atol=1e-12),"link mapping")
    # Use the native proposal operators, not native ownership of stratigraphy.
    with np.errstate(all="warn"):
        r.shear_stress();r.bedload_equation();r.calculate_net_bedload();r.update_bed_elevation()
    qb=np.asarray(r._sed_transp__bedload_rate_link).copy();p=np.asarray(r._sed_transp__bedload_gsd_link).copy();zn=g.at_node["topographic__elevation"].copy()
    req(np.isfinite(qb).all() and np.all(qb[~act]==0) and np.all(qb[b!=0]==0) and np.array_equal(zn[~c],z[~c]),"native boundary")
    for a0,a1 in ((g.at_node["surface_water__depth"],h),(g.at_link["surface_water__depth"],hl),(g.at_link["surface_water__velocity"],u),(r._surface_water__velocity_prev_time_link,prev)):req(np.array_equal(a0,a1),"hydraulic mutation")
    f=qb*DX*DT;Mn,d,amount,gross,res=move(M,f,p,t,hd,c,act)
    prop=frac[np.where(f>=0,t,hd)]
    ctrl,cd,ca,cg,cr=move(M,f,prop,t,hd,c,act)
    dv=(zn-z)*(1-PHI)*DX*DX
    err=float(np.max(np.abs(d.sum(axis=1)-dv)))
    req(err<=1e-10 and np.max(np.abs(Mn.sum(axis=1)-(V+dv)))<=1e-10,"Exner/aggregate")
    req(input_hash(V,z,h,q,hl,u,prev,M,fine)==before,"input mutation")
    diff=float(np.max(np.abs(Mn-ctrl)))
    np.savez_compressed(out/f"{k}_{mode}.npz",V=V,z=z,z_new=zn,h=h,hl=hl,u=u,prev=prev,stock=M,stock_new=Mn,stock_proportional=ctrl,qb=qb,p=p,flow_classes=amount,gross=gross,core=c,tail=t,head=hd)
    return {"step":k,"mode":mode,"max_qb":float(np.max(np.abs(qb))),"gross_m3":float(np.abs(f).sum()),"class_residual_m3":res,"exner_local_m3":err,"selection_diff_m3":diff},qb
def main(dest):
    out=pathlib.Path(dest);out.mkdir();start=time.monotonic()
    rep={"experiment":"P5.1","cases":{},"comparisons":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"]
    rep["fixture"]=fixture();ini,x=retained();RBD=load_provider(out)
    for k in STEPS:
        qbs={}
        for mode in MODES:
            try:
                row,qb=one(RBD,ini,x,k,mode,out);rep["cases"][f"{k}_{mode}"]=row;qbs[mode]=qb
            except Exception as e:traceback.print_exc();fail.append(f"{k}_{mode}:{e!r}")
        if len(qbs)==3:rep["comparisons"][str(k)]={"heterogeneous_mirror_qb_diff":float(np.max(np.abs(qbs["heterogeneous"]-qbs["mirror"]))),"uniform_selection_diff":rep["cases"][f"{k}_uniform"]["selection_diff_m3"]}
    if len(rep["cases"])!=9:fail.append("incomplete")
    if sum(v["uniform_selection_diff"]>1e-8 for v in rep["comparisons"].values())<2:fail.append("selection")
    if sum(v["heterogeneous_mirror_qb_diff"]>1e-8 for v in rep["comparisons"].values())<2:fail.append("composition_feedback")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_GRAIN_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
