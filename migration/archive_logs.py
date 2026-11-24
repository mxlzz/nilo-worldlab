#!/usr/bin/env python3
"""Archive original WorldLab GitHub Actions run logs, not execute experiments."""
import concurrent.futures as cf, hashlib, json, os, pathlib, threading, time, urllib.request, urllib.error, zipfile
ROOT=pathlib.Path.cwd()
MANIFEST=json.loads((ROOT/"evidence/ARTIFACT_MANIFEST.json").read_text())
SOURCE=MANIFEST["source"]
TOKEN=os.environ["GH_TOKEN"]
API="https://api.github.com/repos/"+SOURCE
OUT=ROOT/"evidence"/"logs"
BUDGET=250*1024*1024
SIZE_LIMIT=90*1024*1024
total=0
lock=threading.Lock()
critical={45389348480,45399096372,45399630729,45400499614,45407165828,45413992020,45445405808,45447643891,45449986939,45488606315}
class Redirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        r=super().redirect_request(req,fp,code,msg,headers,newurl)
        if r is not None:r.remove_header("Authorization")
        return r
def archive(run):
    global total
    rid=run["id"];dest=OUT/(str(rid)+".zip");tmp=dest.with_suffix(".part")
    info={"run_id":rid,"name":run.get("name"),"head_sha":run.get("head_sha"),"conclusion":run.get("conclusion"),"url":run.get("html_url")}
    for attempt in range(2):
        try:
            req=urllib.request.Request(API+"/actions/runs/"+str(rid)+"/logs",headers={"Authorization":"Bearer "+TOKEN,"Accept":"application/vnd.github+json","User-Agent":"WorldLab-approved-evidence-archive"})
            with urllib.request.build_opener(Redirect).open(req,timeout=60) as resp,tmp.open("wb") as f:
                h=hashlib.sha256();n=0
                while True:
                    chunk=resp.read(1024*1024)
                    if not chunk:break
                    n+=len(chunk)
                    if n>SIZE_LIMIT:raise ValueError("individual ZIP exceeds safe git limit")
                    f.write(chunk);h.update(chunk)
            with zipfile.ZipFile(tmp) as z:
                if not z.namelist() or z.testzip() is not None:raise ValueError("invalid log ZIP")
            with lock:
                if total+n>BUDGET:raise ValueError("global archival budget exceeded")
                total+=n
                tmp.replace(dest)
            return {**info,"status":"archived_local_sha256","path":dest.relative_to(ROOT).as_posix(),"bytes":n,"sha256":h.hexdigest()}
        except urllib.error.HTTPError as e:
            tmp.unlink(missing_ok=True)
            if e.code in (403,404,410,429) or attempt:
                return {**info,"status":"not_archived","error":"HTTP "+str(e.code)}
            time.sleep(2)
        except Exception as e:
            tmp.unlink(missing_ok=True)
            if attempt:return {**info,"status":"not_archived","error":str(e)}
    return {**info,"status":"not_archived","error":"retry limit"}
OUT.mkdir(parents=True,exist_ok=True)
if any(OUT.glob("*.zip")):raise RuntimeError("Archive exists; reconcile before rerunning")
runs=sorted(MANIFEST["runs"],key=lambda x:(x["id"] not in critical,"isolated matter tests" in (x.get("name") or ""),-x["id"]))
with cf.ThreadPoolExecutor(max_workers=5) as executor: records=list(executor.map(archive,runs))
result={"source":SOURCE,"source_commit":MANIFEST["source_commit"],"run_count":len(records),"archived_count":sum(r["status"]=="archived_local_sha256" for r in records),"archive_bytes":total,"note":"Log ZIPs have locally calculated SHA256 and ZIP CRC validation; GitHub provides no expected digest for these log downloads.","records":records}
(ROOT/"evidence"/"LOG_MANIFEST.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print("WORLDLAB_LOG_ARCHIVE",json.dumps({k:v for k,v in result.items() if k!="records"}))
