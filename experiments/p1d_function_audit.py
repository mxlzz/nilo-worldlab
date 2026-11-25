"""P1-D.1 read-only functional audit with analytical bare control."""
import hashlib,io,json,pathlib,sys,time,traceback,zipfile
import numpy as np
import p1c_live_world as live
from p1c_live_world import ker,thermal
ROOT=live.ROOT
P=ROOT/"evidence/P1-C/C2a/14977708755.zip"
B=ROOT/"evidence/P1-B/B2c/result.json"
P_SHA="fd0ba07f97a2e0f9121bc5cdf6b19ba37d2dcd06a1c751a193fa60d40b271d7d"
DOSE=324000.;CELL=90.;N=60
REQ=ker.req
def sha(x):return hashlib.sha256(x).hexdigest()
def read(z,report,name):
    data=z.read(name);REQ(sha(data)==report["outputs"][name],"input digest "+name);return data
def nc(z,report,name):
    read(z,report,name)
    return thermal.read_nc(z,name)
def stats(red,bare):
    dt=thermal.deltaT(red);bt=thermal.deltaT(bare);delta=bt-dt
    x=thermal.cold(red,bare)
    x.update(min_deficit_K=float(delta.min()),negative_cells=int((delta < -1e-7).sum()),actual_max_K=float(300+dt.max()),bare_max_K=float(300+bt.max()),roi_J=float(red.E_roi.sum()),bare_roi_J=float(bare.E_roi.sum()),grain_J=float(red.E_grain.sum()),rest_J=float(red.E_rest.sum()),export_J=float(red.attrs["E_export_J"]),bare_rest_J=float(bare.E_rest.sum()),bare_export_J=float(bare.attrs["E_export_J"]),bare_ledger_residual_J=ker.total(bare)-float(bare.attrs["E_incident_J"]))
    REQ(np.isfinite(list(x.values())).all(),"nonfinite metrics")
    return x
def main(out):
    out=pathlib.Path(out);out.mkdir()
    start=time.monotonic()
    result={"experiment":"P1-D.1","retrospective":True,"cases":{},"failed_preregistered_gates":[]}
    try:
        REQ(sha(P.read_bytes())==P_SHA,"producer archive")
        bref=json.loads(B.read_text());REQ(not bref["failed_preregistered_gates"],"B2c gates")
        with zipfile.ZipFile(P) as z:
            REQ(z.testzip() is None,"archive CRC")
            a=json.loads(z.read("result.json"));REQ(not a["failed_preregistered_gates"] and set(a["cases"])==set(bref["cases"]) and len(a["cases"])==4,"reference gates")
            for key,v in a["cases"].items():
                seed=int(key.split("_")[0]);records=[];bare=None
                REQ(len(v["epochs"])==7,"epoch count")
                for epoch,rec in enumerate(v["epochs"]):
                    owner=nc(z,a,f"{key}/epoch_{epoch}.nc")
                    REQ(live.digest(owner)==rec["owner_sha"],"owner digest")
                    actual=live.unpack(owner,seed,epoch,"actual")
                    if not epoch:bare=thermal.reduced(live.project(owner),seed,N)
                    else:
                        with np.load(io.BytesIO(read(z,a,f"{key}/flux_{epoch}.npz")),allow_pickle=False) as f:
                            REQ(live.flux_digest({n:f[n] for n in f.files})==rec["flux_sha"],"flux digest")
                            for step in range(3*epoch-2,3*epoch+1):
                                flows={n:f[f"s{step}_{n}"] for n in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
                                bare,_=ker.advance(bare,flows,step)
                                REQ(abs(ker.total(bare)-bare.attrs["E_incident_J"])<=.01,"transport balance")
                    REQ(np.array_equal(bare.soil_mass_kg,actual.soil_mass_kg) and np.array_equal(bare.grain_mass_kg,actual.grain_mass_kg) and np.array_equal(bare.grain,actual.grain),"same stocks")
                    REQ(abs(bare.attrs["mass_export_kg"]-actual.attrs["mass_export_kg"])<=.09216,"same material export")
                    bare.E_roi.values[:]+=CELL
                    bare.attrs["E_incident_J"]=DOSE*(epoch+1)
                    ker.validate(bare,*ker.contract(bare,bare.soil_mass_kg.values))
                    REQ(np.ptp(bare.E_roi.values)<=1e-7 and np.all(bare.E_grain.values==0),"uniform unoccluded control")
                    REQ(actual.attrs["E_incident_J"]==bare.attrs["E_incident_J"] and abs(ker.total(bare)-bare.attrs["E_incident_J"])<=.01,"dose/ledger")
                    x=stats(actual,bare);x["epoch"]=epoch
                    with np.load(io.BytesIO(read(z,a,f"{key}/ray_{epoch}.npz")),allow_pickle=False) as f:
                        x["intercepted_cells"]=int((f["source_id"]>=0).sum())
                        REQ(x["intercepted_cells"]==rec["radiation"]["intercepted_cells"],"ray count")
                    if epoch==0:
                        REQ(np.max(np.abs(thermal.deltaT(bare)-.25))<=1e-7,"B2c bare analytic")
                        ref=bref["cases"][key]["before"]
                        for name in ("cold_cells","largest"):REQ(x[name]==ref[name],"B2c "+name)
                        for name in ("mean_deficit_K","max_deficit_K"):REQ(abs(x[name]-ref[name])<=1e-7,"B2c "+name)
                    x["field_max_difference_J"]=float(np.max(np.abs(bare.E_roi.values-actual.E_roi.values)))
                    if key=="7321_0.1":
                        fr=nc(z,a,f"frozen/epoch_{epoch}.nc")
                        frec=a["controls"]["frozen"]["epochs"][epoch]
                        REQ(live.digest(fr)==frec["owner_sha"] and rec["projection_sha"]==frec["projection_sha"],"frozen lineage")
                        frozen=live.unpack(fr,seed,epoch,"frozen")
                        REQ(frozen.attrs["E_incident_J"]==bare.attrs["E_incident_J"],"frozen dose")
                        x["frozen"]=stats(frozen,bare)
                        x["actual_frozen_roi_max_J"]=float(np.max(np.abs(actual.E_roi.values-frozen.E_roi.values)))
                    records.append(x)
                REQ(max(x["field_max_difference_J"] for x in records[1:])>1e-7,"no detectable contrast")
                result["cases"][key]=records
        REQ(sha(P.read_bytes())==P_SHA,"input mutated")
    except Exception as e:traceback.print_exc();result["failed_preregistered_gates"].append(repr(e))
    result["wall_seconds"]=time.monotonic()-start
    (out/"result.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("WORLD_FUNCTION_JSON "+json.dumps(result,sort_keys=True),flush=True)
    if result["failed_preregistered_gates"]:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
