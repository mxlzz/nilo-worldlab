import hashlib,json,os,pathlib,zipfile
p=pathlib.Path("evidence/ARCH/HAIL-first")
z=p/"original.zip"
sha=lambda x:hashlib.sha256(x).hexdigest()
assert z.stat().st_size==724943
assert sha(z.read_bytes())=="8e54f03cd932b696f0ba6033958a4d17a399157a7fd8742960f652b6ca2d16b0"
with zipfile.ZipFile(z) as f:
    assert f.testzip() is None
    names=[i.filename for i in f.infolist() if not i.is_dir()]
    assert len(names)==len(set(names))
    assert all(not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts for n in names)
    members={n:{"sha256":sha(f.read(n)),"bytes":len(f.read(n))} for n in names}
    grains=[n for n in names if pathlib.PurePosixPath(n).name.lower().startswith("grainz") and n.endswith(".asc")]
    assert grains
    assert any("nrows" not in f.read(n).decode().lower().splitlines()[:6] for n in grains)
    f.extractall(pathlib.Path(os.environ["RUNNER_TEMP"])/"hail-original")
    m={"artifact_id":15044938432,"run_id":45662003426,"zip_sha256":sha(z.read_bytes()),"zip_bytes":z.stat().st_size,"members":members,"grain_files":len(grains)}
(p/"archive_manifest.json").write_text(json.dumps(m,indent=2)+"\n")
print("VERIFIED",len(members),"members; grain tables",len(grains))
