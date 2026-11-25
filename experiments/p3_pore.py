"""P3.2b: pore water bookkeeping during native sediment transfers.
Well-mixed pore storage, release at donor, dry-solid deposition closure.
Not a wetting front reconstruction or new sediment equation.
"""
import hashlib,io,json,pathlib,resource,sys,time,traceback,zipfile
import numpy as np
from landlab.components import SoilInfiltrationGreenAmpt
from p3_exner import grid,req,DX,DT,PHI,FILM

ARCH=pathlib.Path("evidence/P3/P32a/15021698022.zip")
ARCH_SHA="8e75c09ac396d2635dd79e36d67d5a5c13c19f7d427f2f6d18e2d82961d2cd37"
THETA=.15
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def input_hash(*arrays):
    d=hashlib.sha256()
    for a in arrays:
        a=np.asarray(a);d.update(str(a.shape).encode());d.update(a.dtype.str.encode());d.update(a.tobytes())
    return d.hexdigest()
def geometry(V,A):
    T=V/((1-PHI)*A)
    C=PHI*V/(1-PHI)
    return T,C
def domain(V,W,A,core):
    T,C=geometry(V,A)
    F=W/A-THETA*T
    return bool(np.all(T[core]>0) and np.all(F[core]>FILM)),F
def remap(V,W,H,f,t,head,core,active,A):
    args=(V,W,H,f,t,head,core,active);finger=input_hash(*args)
    V,W,H,f,t,head,core,active=map(np.asarray,args)
    n=len(V);m=len(f)
    req(V.shape==W.shape==H.shape==core.shape==(n,) and f.shape==t.shape==head.shape==active.shape==(m,),"shape")
    req(core.dtype==active.dtype==bool and t.dtype.kind in "iu" and head.dtype.kind in "iu","index types")
    req(np.isfinite(A) and A>0 and np.all((t>=0)&(t<n)&(head>=0)&(head<n)&(t!=head)),"area/topology")
    req(all(np.isfinite(x).all() for x in (V,W,H,f)) and np.all(V>=0) and np.all(W>=0) and np.all(H>=0),"stock bounds")
    T,C=geometry(V,A)
    req(np.all(W<=C+1e-12) and all(np.all(x[~core]==0) for x in (V,W,H)),"capacity/halo")
    req(np.all(f[~active]==0),"inactive flux")
    source=np.where(f>=0,t,head);target=np.where(f>=0,head,t)
    absf=np.abs(f)
    req(np.all(absf[~core[source]]==0),"undeclared solid inflow")
    gross=np.zeros(n);incoming=np.zeros(n)
    np.add.at(gross,source,absf);np.add.at(incoming,target,absf)
    out=float(incoming[~core].sum())
    req(np.all(gross<=V),"overdraw")
    frac=np.divide(gross,V,out=np.zeros_like(V),where=V>0)
    release=W*frac
    vp=V-gross+incoming;vp[~core]=0.
    wp=W-release;hp=H+release
    tp,cp=geometry(vp,A)
    req(np.all(vp>=0) and np.all(wp>=0) and np.all(hp>=0) and np.all(wp<=cp+1e-12),"candidate bounds")
    solres=float(vp.sum()+out-V.sum());waterres=float((wp+hp).sum()-(W+H).sum())
    local=float(np.max(np.abs(((wp+hp)-(W+H))[core])))
    req(abs(solres)<=1e-9 and abs(waterres)<=1e-10 and local<=1e-12,"conservation")
    req(input_hash(*args)==finger,"input mutation")
    valid,F=domain(vp,wp,A,core)
    result={"solid_out_m3":out,"solid_residual_m3":solres,"water_residual_m3":waterres,"max_local_water_residual_m3":local,"released_water_m3":float(release.sum()),"min_F_equiv_m":float(np.min(F[core])),"GA_domain_valid":valid,"source_unchanged":True}
    arrays={"V_before":V,"W_before":W,"H_before":H,"f":f,"tail":t,"head":head,"core":core,"active":active,"V_after":vp,"W_after":wp,"H_after":hp,"released":release,"gross":gross,"incoming":incoming,"thickness_after":tp,"capacity_after":cp,"F_equiv_m":F}
    return result,arrays
def analytic(V,W,H,f,t,head,core=None,active=None):
    V,W,H,f=map(lambda a:np.asarray(a,dtype=float),(V,W,H,f))
    t,head=map(lambda a:np.asarray(a,dtype=np.int64),(t,head))
    core=np.ones(len(V),bool) if core is None else np.asarray(core,bool)
    active=np.ones(len(f),bool) if active is None else np.asarray(active,bool)
    return (V,W,H,f,t,head,core,active,1.)
def ga(a,halo_h,K):
    v,w,h=(a[k].copy() for k in ("V_after","W_after","H_after"))
    core=a["core"];A=DX*DX
    valid,F=domain(v,w,A,core);req(valid,"GA equivalent domain")
    g,mask,active,t,head,bsign=grid()
    req(np.array_equal(mask,core) and np.array_equal(t,a["tail"]) and np.array_equal(head,a["head"]),"GA grid")
    h0=np.where(core,h/A,halo_h)
    f0=np.where(core,F,.01)
    g.at_node["surface_water__depth"]=h0.copy()
    g.at_node["soil_water_infiltration__depth"]=f0.copy()
    native=SoilInfiltrationGreenAmpt(g,hydraulic_conductivity=np.where(core,K,0.),soil_bulk_density=1590.,rock_density=2650.,initial_soil_moisture_content=THETA,soil_type="sandy loam",surface_water_minimum_depth=FILM)
    native.run_one_step(DT)
    hn=g.at_node["surface_water__depth"];fn=g.at_node["soil_water_infiltration__depth"]
    raw=(fn-f0)*A
    req(np.isfinite(raw).all() and np.all(raw>=0) and np.all(raw[~core]==0) and np.all(hn[~core]==h0[~core]),"GA native proposal")
    req(np.max(np.abs((h0-hn)*A-raw))<=1e-12,"GA native balance")
    _,cap=geometry(v,A)
    # Explicit finite-capacity boundary for infiltration, not residual correction.
    wp=np.minimum(w+raw,cap)
    transfer=wp-w
    req(np.all(transfer>=0) and np.all(transfer<=raw+1e-12),"GA bounded transfer")
    hp=h-transfer
    req(np.all(hp>=0) and np.all(wp>=0) and np.all(wp<=cap),"GA candidate")
    residual=float((wp+hp).sum()-(w+h).sum())
    local=float(np.max(np.abs(wp+hp-w-h)))
    req(abs(residual)<=1e-10 and local<=1e-12,"GA water conservation")
    if K==0:req(np.array_equal(wp,w) and np.array_equal(hp,h),"GA K0")
    return {"infiltration_m3":float(transfer.sum()),"water_residual_m3":residual,"max_local_water_residual_m3":local,"max_native_transfer_error_m3":float(np.max(np.abs((h0-hn)*A-raw)))},{"GA_W_after":wp,"GA_H_after":hp,"GA_raw":raw,"GA_transfer":transfer}
def main(dest):
    out=pathlib.Path(dest);out.mkdir()
    start=time.monotonic();rep={"experiment":"P3.2b","analytic":{},"retained":{},"controls":{},"failed_preregistered_gates":[]}
    fail=rep["failed_preregistered_gates"]
    def check(name,fun):
        try:fun()
        except Exception as e:
            traceback.print_exc();fail.append(name+":"+repr(e))
    def fixture(name,args,expected,domain_valid):
        r,a=remap(*args)
        rep["analytic"][name]=r;np.savez_compressed(out/(name+".npz"),**a)
        for key,val in expected.items():req(np.allclose(a[key],val,rtol=0,atol=1e-12),name+" "+key)
        req(r["GA_domain_valid"]==domain_valid,"GA domain classification")
        if not domain_valid:
            # Enforce same guard used by ga() before any provider call.
            before=input_hash(a["V_after"],a["W_after"],a["H_after"])
            try:ga(a,np.zeros_like(a["V_after"]),1e-6)
            except ValueError as e:req(str(e)=="GA equivalent domain","domain guard cause")
            else:raise ValueError("invalid GA domain accepted")
            req(before==input_hash(a["V_after"],a["W_after"],a["H_after"]),"domain guard mutated")
            rep["controls"][name+"_GA_rejected"]=True
        if name=="through":req(abs(r["solid_out_m3"]-.06)<=1e-12 and abs(r["released_water_m3"]-.1)<=1e-12,"through/export")
        if name=="zero":req(all(np.array_equal(a[k+"_before"],a[k+"_after"]) for k in ("V","W","H")),"zero exact")
    b=analytic([.6,.6,.6,0],[.2,.2,.2,0],[.1,.1,.1,0],[.12,.12,.06],[0,1,2],[1,2,3],[1,1,1,0])
    check("through",lambda:fixture("through",b,{"V_after":[.48,.6,.66,0],"W_after":[.16,.16,.18,0],"H_after":[.14,.14,.12,0]},True))
    check("saturated",lambda:fixture("saturated",analytic([.6,.6],[.4,.4],[.1,.1],[.3],[0],[1]),{"V_after":[.3,.9],"W_after":[.2,.4],"H_after":[.3,.1],"capacity_after":[.2,.6]},True))
    check("deposition_domain",lambda:fixture("deposition_domain",analytic([.6,.6],[.2,.16],[.1,.1],[.3],[0],[1]),{"V_after":[.3,.9],"W_after":[.1,.16],"H_after":[.2,.1]},False))
    check("exhaustion",lambda:fixture("exhaustion",analytic([.6,.6],[.2,.2],[.1,.1],[.6],[0],[1]),{"V_after":[0.,1.2],"W_after":[0.,.2],"H_after":[.3,.1]},False))
    zero=list(b);zero[3]=np.zeros_like(b[3])
    check("zero",lambda:fixture("zero",tuple(zero),{},True))
    def reject(name,args,reason):
        finger=input_hash(*args[:-1])
        try:remap(*args)
        except ValueError as e:req(str(e)==reason,"wrong rejection "+str(e))
        else:raise ValueError("accepted corruption "+name)
        req(input_hash(*args[:-1])==finger,"rejection mutated")
        rep["controls"][name]=True
    check("overdraw",lambda:reject("overdraw",analytic([.6,.6],[.2,.2],[.1,.1],[.61],[0],[1]),"overdraw"))
    check("undeclared_inflow",lambda:reject("undeclared_inflow",analytic([.6,0],[.2,0],[.1,0],[-.1],[0],[1],[1,0]),"undeclared solid inflow"))
    check("inactive",lambda:reject("inactive",analytic([.6,.6],[.2,.2],[.1,.1],[.1],[0],[1],active=[0]),"inactive flux"))
    check("overcapacity",lambda:reject("overcapacity",analytic([.6,.6],[.41,.2],[.1,.1],[.1],[0],[1]),"capacity/halo"))
    check("nan",lambda:reject("nan",analytic([.6,.6],[.2,.2],[.1,.1],[float("nan")],[0],[1]),"stock bounds"))
    def retained():
        req(sha(ARCH)==ARCH_SHA,"archive SHA")
        with zipfile.ZipFile(ARCH) as z:
            req(z.testzip() is None,"archive CRC")
            old=json.loads(z.read("result.json"));req(not old["failed_preregistered_gates"],"reference gates")
            for name in ("live_40","live_80","live_120","imposed_export","dry"):
                def sample(name=name):
                    n=name+".npz";payload=z.read(n);req(hashlib.sha256(payload).hexdigest()==old["outputs"][n],"reference hash")
                    with np.load(io.BytesIO(payload),allow_pickle=False) as p:x={k:p[k].copy() for k in p.files}
                    core=x["core"];A=DX*DX
                    args=(x["stock_before"],np.where(core,.04*A,0.),np.where(core,x["h_node"]*A,0.),x["qb"]*DX*DT,x["tail"],x["head"],core,x["active"],A)
                    r,a=remap(*args)
                    candidate_core_error=float(np.max(np.abs(a["V_after"][core]-x["stock_candidate"][core])))
                    owner_error=float(np.max(np.abs(a["V_after"]-x["stock_accepted"])))
                    halo_flux_error=float(np.max(np.abs(x["stock_candidate"][~core]-a["incoming"][~core])))
                    halo_sum_error=abs(float(x["stock_candidate"][~core].sum())-r["solid_out_m3"])
                    reference_out_error=abs(old["cases"][name]["out_solid_m3"]-r["solid_out_m3"])
                    req(candidate_core_error<=1e-12 and owner_error<=1e-12,"stock candidate/accepted owner mismatch")
                    req(np.all(a["V_after"][~core]==0) and np.all(x["stock_accepted"][~core]==0),"halo owner stock")
                    req(max(halo_flux_error,halo_sum_error,reference_out_error)<=1e-12,"halo boundary flux mismatch")
                    r["reference_comparison"]={"candidate_core_max_error_m3":candidate_core_error,"accepted_owner_max_error_m3":owner_error,"halo_flux_max_error_m3":halo_flux_error,"halo_sum_error_m3":halo_sum_error,"reference_export_error_m3":reference_out_error,"halo_candidate_m3":float(x["stock_candidate"][~core].sum())}
                    req(np.max(np.abs(a["gross"]-x["gross_out"]))<=1e-12,"gross mismatch")
                    req(np.max(np.abs((x["bedrock"]+a["thickness_after"]-x["z_after"])[core]))<=1e-12,"geometry mismatch")
                    req(r["GA_domain_valid"],"retained GA domain")
                    if name.startswith("live_"):req(r["released_water_m3"]>0,"live pore water release")
                    if name=="imposed_export":req(abs(r["solid_out_m3"]-.00002)<=1e-12,"imposed export")
                    if name=="dry":req(all(np.array_equal(a[k+"_before"],a[k+"_after"]) for k in ("V","W","H")),"dry exact")
                    ga_result,ga_arrays=ga(a,x["h_node"],1e-6)
                    r["GA"]=ga_result;r["source_npz_sha256"]=old["outputs"][n]
                    if name=="live_120":
                        req(ga_result["infiltration_m3"]>0,"GA live infiltration")
                        k0,ga0=ga(a,x["h_node"],0.)
                        req(k0["infiltration_m3"]==0,"GA K0");rep["controls"]["GA_K0"]=True
                        ga_arrays.update({"K0_W_after":ga0["GA_W_after"],"K0_H_after":ga0["GA_H_after"]})
                    np.savez_compressed(out/(name+".npz"),**a,**ga_arrays)
                    rep["retained"][name]=r
                check(name,sample)
    check("retained_archive",retained)
    req_names={"deposition_domain_GA_rejected","exhaustion_GA_rejected","overdraw","undeclared_inflow","inactive","overcapacity","nan","GA_K0"}
    if len(rep["analytic"])!=5 or len(rep["retained"])!=5 or set(rep["controls"])!=req_names:fail.append("incomplete cases or controls")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={p.name:sha(p) for p in out.iterdir() if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n")
    print("WORLD_PORE_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
