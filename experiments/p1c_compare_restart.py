"""Compare C.2b only after consumer finishes; no evolution in this job."""
import hashlib,io,json,pathlib,sys,zipfile
import numpy as np

def req(v,msg):
    if not bool(v):raise ValueError(msg)
def sha(data):return hashlib.sha256(data).hexdigest()
def main(producer,consumer,out):
    root=pathlib.Path(consumer);target=pathlib.Path(out);target.mkdir()
    result={"experiment":"P1-C.2b-compare","cases":{},"failed_preregistered_gates":[]}
    try:
        p=pathlib.Path(producer);req(sha(p.read_bytes())=="fd0ba07f97a2e0f9121bc5cdf6b19ba37d2dcd06a1c751a193fa60d40b271d7d","producer digest")
        with zipfile.ZipFile(p) as z:
            a=json.loads(z.read("result.json"));b=json.loads((root/"result.json").read_text())
            req(not a["failed_preregistered_gates"] and not b["failed_preregistered_gates"],"producer/consumer gates")
            req(set(a["cases"])==set(b["cases"]) and len(b["cases"])==4,"case set")
            req(len(b["outputs"])==36,"suffix output count")
            for name,h in b["outputs"].items():req(sha((root/name).read_bytes())==h,"consumer digest "+name)
            for key,v in b["cases"].items():
                ref=a["cases"][key]
                req(v["history"]==ref["history"] and v["formation_steps"]==ref["formation_steps"],"history "+key)
                req(v["checkpoint_owner_sha"]==ref["epochs"][3]["owner_sha"],"checkpoint "+key)
                req(v["new_bullet_steps"]==2160 and v["new_thermal_steps"]==18 and v["new_soil_steps"]==9 and len(v["epochs"])==3,"suffix steps")
                req(all(v["controls"].values()) and len(v["controls"])==7,"corruption controls")
                records=[]
                for x,y,ev,expected_ev in zip(v["epochs"],ref["epochs"][4:],v["requests"],ref["requests"][4:]):
                    epoch=x["epoch"];expected={k:w for k,w in y.items() if k!="flux_comparison"}
                    req(x==expected,"epoch logical mismatch "+key+":"+str(epoch))
                    req({k:w for k,w in ev.items() if k!="files"}=={k:w for k,w in expected_ev.items() if k!="files"},"worker event")
                    byte={}
                    for kind,ext in [("epoch","nc"),("flux","npz"),("ray","npz")]:
                        name=f"{key}/{kind}_{epoch}.{ext}";src=z.read(name);dst=(root/name).read_bytes()
                        req(sha(src)==a["outputs"][name],"producer output digest")
                        byte[kind]=src==dst
                        if ext=="npz":
                            with np.load(io.BytesIO(src),allow_pickle=False) as s,np.load(io.BytesIO(dst),allow_pickle=False) as t:
                                req(set(s.files)==set(t.files),"array names")
                                for n in s.files:req(s[n].dtype==t[n].dtype and s[n].shape==t[n].shape and s[n].tobytes()==t[n].tobytes(),"array bit mismatch "+n)
                        else:
                            # The logical hash was recomputed and validated independently by each runner.
                            req(x["owner_sha"]==y["owner_sha"],"owner digest")
                    records.append({"epoch":epoch,"owner_sha":x["owner_sha"],"trace_sha":x["trace_sha"],"byte_exact_files":byte,"exact":True})
                result["cases"][key]={"suffix":records,"corruption_controls":v["controls"],"new_bullet_steps":v["new_bullet_steps"],"new_thermal_steps":v["new_thermal_steps"],"new_soil_steps":v["new_soil_steps"]}
    except Exception as e:result["failed_preregistered_gates"].append(repr(e))
    (target/"result.json").write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
    print("WORLD_COMPARE_JSON "+json.dumps(result,sort_keys=True),flush=True)
    if result["failed_preregistered_gates"]:raise SystemExit(1)
if __name__=="__main__":main(*sys.argv[1:])
