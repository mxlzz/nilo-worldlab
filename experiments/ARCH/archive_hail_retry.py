import hashlib,json,pathlib,zipfile
p=pathlib.Path("evidence/ARCH/HAIL-retry")
z=p/"retry.zip"
sha=lambda b:hashlib.sha256(b).hexdigest()
expected="5ed01c2b1e8a161dda0fc1f14a6cc17b013c20d5b188c946b670e716ef952e6b"
assert z.stat().st_size==1452320 and sha(z.read_bytes())==expected
with zipfile.ZipFile(z) as f:
    assert f.testzip() is None
    names=[i.filename for i in f.infolist() if not i.is_dir()]
    assert len(names)==len(set(names))
    assert all(not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts for n in names)
    members={n:{"sha256":sha(f.read(n)),"bytes":len(f.read(n))} for n in names}
    assert f.read("gate.txt").strip()==b"PASS"
    assert f.read("source.txt").strip()==b"5913d54b9bf5cf97a417ae30c98b39e624bac013"
    h1=json.loads(f.read("hash1.json"));h2=json.loads(f.read("hash2.json"))
    assert h1==h2 and h1
    for i in (1,2):
        prefix=f"run{i}/"
        paths={n[len(prefix):]:m["sha256"] for n,m in members.items() if n.startswith(prefix)}
        assert paths==h1
        assert b"THE SIMULATION IS FINISHED!" in f.read(f"run{i}.log")
    old=json.loads(f.read("original-hashes.json"))
    comparison=f.read("original-comparison.txt").decode().strip()
    assert (comparison=="MATCH")==(h1==old)
    oldmanifest=json.loads(pathlib.Path("evidence/ARCH/HAIL-first/archive_manifest.json").read_text())
    assert old=={n[5:]:m["sha256"] for n,m in oldmanifest["members"].items() if n.startswith("run1/")}
    time={n:f.read(n).decode() for n in ["build.time","run1.time","run2.time"]}
    summary={"source":f.read("source.txt").decode().strip(),"same_host_exact":True,"original_comparison":comparison,"files_per_run":len(h1),"waterdepth_rasters":sum("waterdepth" in n.lower() and n.endswith(".asc") for n in h1),"times":time,"environment":f.read("environment.txt").decode()}
    for n in ("run1.time","run2.time"):
        assert "Exit status: 0" in time[n]
manifest={"artifact_id":15088406154,"run_id":45756683979,"zip_sha256":expected,"zip_bytes":z.stat().st_size,"members":members,"summary":summary}
(p/"archive_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
print("VERIFIED",len(members),"members",json.dumps(summary))