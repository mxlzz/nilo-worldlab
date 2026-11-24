#!/usr/bin/env python3
"""One-time transfer of approved WorldLab paths; no scientific execution."""
import hashlib, json, os, pathlib, shutil, subprocess, tempfile, urllib.request, urllib.error, time
ROOT=pathlib.Path.cwd()
SOURCE="mxlzz/Nilo_public"
BRANCH="research/world-foundations-2026-11"
SHA="0448ea6e5e02d77142cec270ea99ac0773493fc1"
API="https://api.github.com/repos/"+SOURCE
TOKEN=os.environ.get("GH_TOKEN","")
def cmd(*a, cwd=None):
    return subprocess.check_output(a,cwd=cwd).decode().strip()
def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
def get_json(url,auth=True):
    headers={"Accept":"application/vnd.github+json","User-Agent":"WorldLab-approved-migration","X-GitHub-Api-Version":"2022-11-28"}
    if auth and TOKEN: headers["Authorization"]="Bearer "+TOKEN
    req=urllib.request.Request(url,headers=headers)
    try:
        with urllib.request.urlopen(req,timeout=60) as f:return json.load(f)
    except urllib.error.HTTPError:
        if auth:return get_json(url,False)
        raise
def download(url, dest, expected_digest=None):
    # Authenticate only to GitHub API; redirects to storage must not receive the token.
    req=urllib.request.Request(url,headers={"Accept":"application/vnd.github+json","Authorization":"Bearer "+TOKEN,"User-Agent":"WorldLab-migration"})
    class SafeRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self,req,fp,code,msg,headers,newurl):
            r=super().redirect_request(req,fp,code,msg,headers,newurl)
            if r is not None:r.remove_header("Authorization")
            return r
    opener=urllib.request.build_opener(SafeRedirect)
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_suffix(".part")
    try:
        with opener.open(req,timeout=120) as response, tmp.open("wb") as out:
            h=hashlib.sha256();size=0
            while True:
                block=response.read(1024*1024)
                if not block:break
                size+=len(block)
                if size>90*1024*1024:raise ValueError("individual file exceeds Git safe limit")
                out.write(block);h.update(block)
        digest=h.hexdigest()
        if expected_digest and expected_digest.removeprefix("sha256:")!=digest:raise ValueError("archive digest mismatch")
        tmp.replace(dest)
        return {"status":"archived_verified_digest" if expected_digest else "archived_local_sha256_only","sha256":digest,"bytes":size,"path":dest.relative_to(ROOT).as_posix()}
    except Exception:
        tmp.unlink(missing_ok=True);raise
def pages(url,key):
    out=[]
    for page in range(1,101):
        data=get_json(url+("&" if "?" in url else "?")+"per_page=100&page="+str(page))
        rows=data[key];out.extend(rows)
        if len(rows)<100:return out
    raise RuntimeError("pagination cap reached; inventory incomplete")
def source_ref():
    return get_json(API+"/git/refs/heads/"+BRANCH)["object"]["sha"]
if source_ref()!=SHA:raise RuntimeError("source HEAD changed; reconcile before copying")
if (ROOT/"inherited").exists():raise RuntimeError("inherited directory already exists; do not overwrite")
with tempfile.TemporaryDirectory() as tmp:
    repo=pathlib.Path(tmp)
    cmd("git","init","-q",cwd=repo)
    cmd("git","remote","add","origin","https://github.com/"+SOURCE+".git",cwd=repo)
    cmd("git","fetch","--depth=1","origin",SHA,cwd=repo)
    cmd("git","checkout","-q","--detach","FETCH_HEAD",cwd=repo)
    assert cmd("git","rev-parse","HEAD",cwd=repo)==SHA
    tree=cmd("git","rev-parse","HEAD^{tree}",cwd=repo)
    raw=subprocess.check_output(["git","ls-tree","-r","-z",SHA],cwd=repo)
    files=[]
    for line in raw.split(b"\0"):
        if not line:continue
        meta,path=line.split(b"\t",1)
        mode,typ,blob=meta.decode().split()
        p=path.decode()
        if not (p.startswith("WORLD_RESEARCH/") or p.startswith(".github/workflows/world-")):continue
        if typ!="blob" or mode not in ("100644","100755"):raise ValueError("unexpected nonregular path: "+p)
        source=repo/p
        dest=ROOT/"inherited"/p
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,dest)
        actual=cmd("git","hash-object","--no-filters",str(dest))
        if actual!=blob:raise ValueError("blob mismatch: "+p)
        files.append({"source":p,"destination":"inherited/"+p,"git_blob_sha":blob,"mode":mode,"size":dest.stat().st_size})
    if len(files)!=314:raise ValueError("unexpected selected file count "+str(len(files)))
    dump(ROOT/"MIGRATION_MANIFEST.json",{"source":SOURCE,"branch":BRANCH,"source_commit":SHA,"source_tree":tree,"count":len(files),"world_research":sum(x["source"].startswith("WORLD_RESEARCH/") for x in files),"historical_workflows":sum(x["source"].startswith(".github/workflows/") for x in files),"files":files,"status":"all_selected_blobs_verified_locally","history":"Original DAG retained in source; not imported; historical workflows inactive under inherited/.github/workflows"})
# Enumerate relevant public workflow runs and artifacts. No scientific tests are launched.
records=[]; errors=[]; total=0; budget=250*1024*1024
try:
    runs=pages(API+"/actions/runs?branch=research%2Fworld-foundations-2026-11","workflow_runs")
except Exception as e:
    runs=[];errors.append({"scope":"run enumeration","error":str(e)})
# Known critical runs supplement list if API paging fails.
critical=[45389348480,45399096372,45399630729,45400499614,45407165828,45413992020,45445405808,45447643891,45449986939,45488606315]
seen={x["id"] for x in runs}
for rid in critical:
    if rid not in seen:
        try:runs.append(get_json(API+"/actions/runs/"+str(rid)))
        except Exception as e:errors.append({"run_id":rid,"error":str(e)})
runs=[r for r in runs if r.get("head_branch")==BRANCH and ("world-" in r.get("path","") or "World research" in r.get("name",""))]
# Save all available run metadata (no credentials), and artifacts for every run.
# Prioritize expiring critical inputs, then most recent experiments.
runs.sort(key=lambda r:(r["id"] not in critical,critical.index(r["id"]) if r["id"] in critical else -r["id"]))
for run in runs:
    rid=run["id"]
    rec={k:run.get(k) for k in ["id","name","head_sha","head_branch","status","conclusion","html_url","created_at","updated_at","path"]}
    rec["artifacts"]=[]
    try: arts=pages(API+"/actions/runs/"+str(rid)+"/artifacts","artifacts")
    except Exception as e:
        rec["artifact_listing_error"]=str(e);records.append(rec);continue
    for a in arts:
        entry={k:a.get(k) for k in ["id","name","size_in_bytes","expired","expires_at","digest","archive_download_url"]}
        if a.get("expired"):entry["archive"]={"status":"expired_not_recovered"}
        elif not a.get("size_in_bytes") or a["size_in_bytes"]>90*1024*1024 or total+a["size_in_bytes"]>budget:
            entry["archive"]={"status":"size_budget_not_archived"}
        else:
            try:
                result=download(API+"/actions/artifacts/"+str(a["id"])+"/zip",ROOT/"evidence"/"archives"/str(rid)/(str(a["id"])+".zip"),a.get("digest"))
                total+=result["bytes"];entry["archive"]=result
            except Exception as e:entry["archive"]={"status":"download_failed","error":str(e)}
        rec["artifacts"].append(entry)
    records.append(rec)
dump(ROOT/"evidence"/"ARTIFACT_MANIFEST.json",{"source":SOURCE,"source_commit":SHA,"scope":"world research workflow artifacts on original branch; metadata plus ZIP bytes if available","runs":records,"errors":errors,"archived_bytes":total,"budget_bytes":budget})
if source_ref()!=SHA:raise RuntimeError("source HEAD changed during migration; do not publish snapshot as current")
print("WORLDLAB_MIGRATION",json.dumps({"files":len(files),"world_research":254,"historical_workflows":60,"runs":len(records),"archived":sum(a["archive"]["status"].startswith("archived") for r in records for a in r["artifacts"]),"artifact_failures":sum(not a["archive"]["status"].startswith("archived") for r in records for a in r["artifacts"]),"bytes":total,"listing_errors":len(errors)}))
