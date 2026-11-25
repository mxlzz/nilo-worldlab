import hashlib, json, math, pathlib, sys
root=pathlib.Path(sys.argv[1]); result={}
for p in sorted(root.rglob("*")):
    if not p.is_file(): continue
    data=p.read_bytes()
    result[str(p.relative_to(root))]=hashlib.sha256(data).hexdigest()
    if p.suffix.lower()==".asc":
        lines=p.read_text().splitlines()
        hdr={x.split()[0].lower():float(x.split()[1]) for x in lines[:6]}
        vals=[float(x) for line in lines[6:] for x in line.split()]
        assert len(vals)==int(hdr["nrows"]*hdr["ncols"]), str(p)
        assert all(math.isfinite(v) for v in vals), str(p)
        if "waterdepth" in p.name.lower():
            assert all(v>=0 or v==hdr.get("nodata_value") for v in vals),str(p)
    if p.suffix.lower()==".dat":
        rows=[l.split() for l in p.read_text().splitlines() if l.strip()]
        assert rows, str(p)
        # upstream header optional; verify at least one numeric row
        parsed=[]
        for row in rows:
            try: parsed.append([float(x) for x in row])
            except ValueError: continue
        assert parsed and all(math.isfinite(v) for row in parsed for v in row), str(p)
assert any("waterdepth" in k.lower() and k.lower().endswith(".asc") for k in result),result.keys()
assert any(k.lower().endswith(".dat") for k in result),result.keys()
print(json.dumps(result,sort_keys=True,indent=2))
