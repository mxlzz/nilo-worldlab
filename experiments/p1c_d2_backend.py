"""C.2b-D2: identical-input thermal diagnostic; no live mechanics or Landlab."""
import hashlib,io,json,os,pathlib,platform,sys,zipfile
from importlib.metadata import version
import numpy as np
from threadpoolctl import threadpool_info
import p1c_live_world as live
from p1c_live_world import ker,thermal,b
ROOT=live.ROOT
P=ROOT/"evidence/P1-C/C2a/14977708755.zip"
Q=ROOT/"evidence/P1-C/C2a/14977883611.zip"
PS="fd0ba07f97a2e0f9121bc5cdf6b19ba37d2dcd06a1c751a193fa60d40b271d7d"
QS="a81ecf3bb6091b2e4a0578baaa3b84702a1e40f2ee8db93e7a91f1577e1f5be6"
def sha(data):return hashlib.sha256(data).hexdigest()
def req(ok,msg):
    if not ok:raise ValueError(msg)
def main(out,mode):
    out=pathlib.Path(out);out.mkdir(parents=True)
    req(sha(P.read_bytes())==PS and sha(Q.read_bytes())==QS,"archive digest")
    req(os.environ.get("OPENBLAS_CORETYPE")==("Nehalem" if mode=="fixed" else None),"core env")
    req(all(os.environ.get(k)=="1" for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")),"threads")
    result={"mode":mode,"cpu":next((x.split(":",1)[1].strip() for x in pathlib.Path("/proc/cpuinfo").read_text().splitlines() if x.startswith("model name")),"unknown"),"python":platform.python_version(),"versions":{n:version(n) for n in ("numpy","scipy","xarray","fipy","threadpoolctl","pybullet")},"cases":{}}
    with zipfile.ZipFile(P) as z,zipfile.ZipFile(Q) as q:
        for key in ("7321_0.0","7321_0.1"):
            seed=7321
            a=thermal.read_nc(z,key+"/epoch_3.nc");ck=thermal.read_nc(q,key+"/owner.nc")
            req(live.digest(a)==live.digest(ck),"checkpoint != epoch3")
            base=live.project(thermal.read_nc(z,key+"/epoch_4.nc"))
            ref=thermal.read_nc(z,key+"/epoch_4.nc")
            red=live.unpack(ck,seed,3,"actual")
            hashes={"checkpoint":live.digest(ck),"base":live.digest(base),"input":ker.digest(red)}
            with np.load(io.BytesIO(z.read(key+"/flux_4.npz")),allow_pickle=False) as f:
                hashes["flows"]=live.flux_digest({n:f[n] for n in f.files})
                for step in (10,11,12):
                    flow={n:f[f"s{step}_{n}"] for n in ("before","after","erosion","deposition","influx","outflux","receiver","order","core")}
                    red,_=ker.advance(red,flow,step)
                    hashes[f"transport_{step}"]=ker.digest(red)
            red,hits,metrics=live.irradiate(red,base,seed,"actual")
            owner=live.publish(base,red,seed,4,"actual")
            hashes.update(irradiation=ker.digest(red),owner=live.digest(owner),reference=live.digest(ref),rays=sha(hits.tobytes()))
            owner.to_netcdf(out/(key+".nc"),engine="scipy")
            req(live.digest(b.read_owner(out/(key+".nc")))==hashes["owner"],"roundtrip")
            result["cases"][key]={"hashes":hashes,"metrics":metrics,"max_roi_delta_J":float(np.abs(owner.E_roi.values-ref.E_roi.values).max()),"max_rest_delta_J":float(np.abs(owner.E_coarse_rest.values-ref.E_coarse_rest.values).max())}
    pools=threadpool_info()
    result["pools"]=pools
    bl=[p for p in pools if p.get("internal_api")=="openblas"]
    result["profile_verified"]=bool(bl) and all(p.get("num_threads")==1 and (mode!="fixed" or str(p.get("architecture","")).lower()=="nehalem") for p in bl)
    result["outputs"]={p.name:sha(p.read_bytes()) for p in out.glob("*.nc")}
    (out/"result.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("WORLD_D2_JSON "+json.dumps(result,sort_keys=True),flush=True)
if __name__=="__main__":main(*sys.argv[1:])
