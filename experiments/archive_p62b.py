import hashlib,json,pathlib,zipfile
p=pathlib.Path("evidence/P6/P62B")
z=p/"p62-territory.zip"
sha=lambda b:hashlib.sha256(b).hexdigest()
expected="e106f2a471f455d5824be3f8b74bd70b3f456af1577017e15e6cc6a573443460"
assert z.stat().st_size==437836 and sha(z.read_bytes())==expected
log=json.loads((p/"readback.json").read_text())
with zipfile.ZipFile(z) as f:
    assert f.testzip() is None
    names=[i.filename for i in f.infolist() if not i.is_dir()]
    assert len(names)==len(set(names))
    assert all(not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts for n in names)
    members={n:{"sha256":sha(f.read(n)),"bytes":len(f.read(n))} for n in names}
    a=json.loads(f.read("reference/result.json"));b=json.loads(f.read("consumer/seed26061104/result.json"));c=json.loads(f.read("compare.json"))
    assert a==log["WORLD_TERRITORY_JSON"] and b==log["WORLD_HORIZON_R_RESTART_JSON"] and c==log["WORLD_HORIZON_R_COMPARE_JSON"]
    assert not a["failed_preregistered_gates"] and all(v["steps"]==240 for v in a["cases"].values())
    assert b["steps"]==120 and len(b["controls"])==7 and all(b["controls"].values()) and c["seed26061104"]["suffix_exact"]
    for prefix,report in (("reference/",a),("consumer/seed26061104/",b)):
        for n,h in report["outputs"].items():assert members[prefix+n]["sha256"]==h,(prefix,n)
    m=json.loads(f.read("checkpoint/seed26061104/manifest.json"))
    assert m["step"]==120 and m["owner_sha"]==b["input_owner_sha"]
    assert members["checkpoint/seed26061104/owner.nc"]["sha256"]==m["sha256"]
    assert m["owner_sha"]==a["cases"]["seed26061104"]["snapshots"]["120"]["owner_sha"]
    assert b["final"]["owner_sha"]==a["cases"]["seed26061104"]["snapshots"]["240"]["owner_sha"]
    assert members["consumer/seed26061104/step0240.nc"]["sha256"]==b["final"]["sha256"]
manifest={"artifact_id":15125444846,"run_id":45838934630,"zip_sha256":expected,"zip_bytes":z.stat().st_size,"members":members,"reports_equal":True,"checkpoint_verified":True,"output_paths_verified":True}
(p/"archive_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
print("VERIFIED",len(members),"members; reports, checkpoint and output paths")
