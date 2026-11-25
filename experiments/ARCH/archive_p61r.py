import hashlib, json, pathlib, zipfile
p=pathlib.Path("evidence/P6/P61R")
z=p/"p6-horizon-r.zip"
sha=lambda b:hashlib.sha256(b).hexdigest()
expect="a1309b755e908093755f3f29bbc83a5ab9820f741195c1eb2f88c6efa8e3adf9"
assert z.stat().st_size==1369091 and sha(z.read_bytes())==expect
log=json.loads((p/"readback.json").read_text())
with zipfile.ZipFile(z) as f:
    assert f.testzip() is None
    names=[i.filename for i in f.infolist() if not i.is_dir()]
    assert len(names)==len(set(names))
    assert all(not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts for n in names)
    members={n:{"sha256":sha(f.read(n)),"bytes":len(f.read(n))} for n in names}
    expected=list(log["producer"]["outputs"].values())
    expected += [v for c in log["consumers"] for v in c["outputs"].values()]
    actual={m["sha256"] for m in members.values()}
    missing=sorted(set(expected)-actual)
    assert not missing,missing
    matches=[n for n in names if n.endswith("compare.json")]
    assert len(matches)==1 and json.loads(f.read(matches[0]))==log["compare"], matches
manifest={"artifact_id":15031613248,"run_id":45633814018,"zip_sha256":expect,"zip_bytes":z.stat().st_size,"members":members,"logged_output_hashes_all_present":True,"compare_json_matches":True}
(p/"archive_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
print("VERIFIED",len(members),"members")
