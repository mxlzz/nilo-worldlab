"""S2 retrospective portability, not prospective confirmation. S1 exact comparator unchanged."""
import hashlib,io,json,pathlib,sys,zipfile
import numpy as np
import xarray as xr
P_SHA="fd0ba07f97a2e0f9121bc5cdf6b19ba37d2dcd06a1c751a193fa60d40b271d7d"
A_SHA="20d1856b2e76aa5ea4b8a9837a1833e87bca1d498a04d170178bd4fd99d83686"
LOCAL=1e-7;L1=1e-4;LEDGER=.01
def req(ok,msg):
    if not bool(ok):raise ValueError(msg)
def sha(x):return hashlib.sha256(x).hexdigest()
def exact(a,b):
    x,y=np.asarray(a),np.asarray(b)
    return x.dtype==y.dtype and x.shape==y.shape and x.tobytes()==y.tobytes()
def nc(data):
    with xr.open_dataset(io.BytesIO(data),engine="scipy") as x:return x.load()
def compare(x,y):
    req(set(x.variables)==set(y.variables) and dict(x.sizes)==dict(y.sizes),"schema")
    req(set(x.attrs)==set(y.attrs),"attrs")
    for k in x.attrs:
        if k!="p1c_E_export_J":req(x.attrs[k]==y.attrs[k],"attr "+k)
    delta={}
    for k in x.variables:
        a,b=x[k],y[k]
        req(a.dims==b.dims and a.dtype==b.dtype and a.attrs==b.attrs,"metadata "+k)
        if k in ("E_roi","E_coarse_rest"):
            req(np.isfinite(a).all() and np.isfinite(b).all(),"nonfinite "+k)
            d=np.abs(a.values-b.values)
            delta[k]={"max":float(d.max()),"L1":float(d.sum())}
            req(delta[k]["max"]<=LOCAL,"local "+k)
        else:req(exact(a.values,b.values),"exact "+k)
    export=abs(float(x.attrs["p1c_E_export_J"])-float(y.attrs["p1c_E_export_J"]))
    req(np.isfinite(export),"nonfinite export")
    total=export+sum(d["L1"] for d in delta.values())
    req(total<=L1,"combined L1")
    balances=[]
    for z in (x,y):
        e=sum(float(z[k].sum()) for k in ("E_roi","E_coarse_rest","E_grain"))+float(z.attrs["p1c_E_export_J"])
        residual=e-float(z.attrs["p1c_E_incident_J"])
        req(abs(residual)<=LEDGER,"ledger")
        # Coarse halo and empty grain stocks cannot hold energy.
        req(np.all(z.E_coarse_rest.values[z.soil_mass_kg.values==0]==0),"empty soil energy")
        req(np.all(z.E_grain.values[z.grain_mass_kg.values==0]==0),"empty grain energy")
        balances.append(residual)
    return {"fields":delta,"export_delta_J":export,"combined_L1_J":total,"ledger_residuals_J":balances}
def metrics(a,b):
    req(set(a)==set(b),"metric schema")
    for k,v in a.items():
        w=b[k]
        if k in ("incident_J","grain_added_J","grain_J","intercepted_cells"):req(v==w,"exact metric "+k)
        elif k.endswith("_J"):req(np.isfinite(v) and np.isfinite(w) and abs(v-w)<=L1,"energy metric "+k)
        elif k.endswith("_K"):req(np.isfinite(v) and np.isfinite(w) and abs(v-w)<=LOCAL,"temperature "+k)
        else:req(v==w,"metric "+k)
def rejected(a,b):
    try:compare(a,b)
    except ValueError:return True
    return False
def main(p,a,s1,out):
    out=pathlib.Path(out);out.mkdir()
    report={"experiment":"P1-C.2b-S2","retrospective":True,"cases":{},"failed_preregistered_gates":[]}
    try:
        req(sha(pathlib.Path(p).read_bytes())==P_SHA and sha(pathlib.Path(a).read_bytes())==A_SHA,"input archive digest")
        exact_report=json.loads((pathlib.Path(s1)/"result.json").read_text())
        req(not exact_report["failed_preregistered_gates"] and len(exact_report["cases"])==4,"S1 gate")
        with zipfile.ZipFile(p) as pz,zipfile.ZipFile(a) as az:
            x=json.loads(pz.read("result.json"));y=json.loads(az.read("result.json"))
            req(not x["failed_preregistered_gates"] and not y["failed_preregistered_gates"],"prior internal gates")
            req(set(x["cases"])==set(y["cases"]) and len(y["cases"])==4,"cases")
            req(len(y["outputs"])==36,"output count")
            for n,h in y["outputs"].items():req(sha(az.read(n))==h,"consumer digest")
            for key,v in y["cases"].items():
                ref=x["cases"][key]
                req(v["history"]==ref["history"] and v["formation_steps"]==ref["formation_steps"] and v["checkpoint_owner_sha"]==ref["epochs"][3]["owner_sha"],"checkpoint/history")
                req(v["new_bullet_steps"]==2160 and v["new_soil_steps"]==9 and v["new_thermal_steps"]==18 and len(v["epochs"])==len(v["requests"])==3,"suffix")
                req(len(v["controls"])==7 and all(v["controls"].values()),"corruption")
                rr=[]
                for m,n,ev,refev in zip(v["epochs"],ref["epochs"][4:],v["requests"],ref["requests"][4:]):
                    epoch=m["epoch"]
                    req(set(m)==set(n)-{"flux_comparison"},"epoch keys")
                    for k in m:
                        if k not in ("owner_sha","radiation","transport"):req(m[k]==n[k],"epoch exact "+k)
                    req({k:v for k,v in ev.items() if k!="files"}=={k:v for k,v in refev.items() if k!="files"},"worker exact")
                    metrics(m["radiation"],n["radiation"])
                    req(len(m["transport"])==len(n["transport"])==3,"transport count")
                    for u,w in zip(m["transport"],n["transport"]):metrics(u,w)
                    for kind,ext in (("epoch","nc"),("flux","npz"),("ray","npz")):
                        name=f"{key}/{kind}_{epoch}.{ext}"
                        src,dst=pz.read(name),az.read(name)
                        req(sha(src)==x["outputs"][name],"producer digest")
                        if ext=="npz":
                            with np.load(io.BytesIO(src),allow_pickle=False) as u,np.load(io.BytesIO(dst),allow_pickle=False) as w:
                                req(set(u.files)==set(w.files) and all(exact(u[k],w[k]) for k in u.files),"arrays exact")
                        else:
                            orig,port=nc(src),nc(dst)
                            result=compare(orig,port)
                            plus=port.copy(deep=True);plus.E_roi.values.flat[0]+=1.
                            comp=port.copy(deep=True);i=int(np.argmax(comp.E_roi.values));j=(i+1)%comp.E_roi.size
                            req(comp.E_roi.values.flat[i]>1,"control stock")
                            comp.E_roi.values.flat[i]-=1.;comp.E_roi.values.flat[j]+=1.
                            req(rejected(orig,plus) and rejected(orig,comp),"portability corruption control")
                    rr.append({"epoch":epoch,**result,"one_J_rejected":True,"compensated_rejected":True})
                report["cases"][key]=rr
    except Exception as e:report["failed_preregistered_gates"].append(repr(e))
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print("WORLD_PORTABILITY_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if report["failed_preregistered_gates"]:raise SystemExit(1)
if __name__=="__main__":main(*sys.argv[1:])
