import json,pathlib
p=pathlib.Path("arch04-out")
a=json.loads((p/"producer/active-bmi/rows.json").read_text())
b=json.loads((p/"consumer/suffix.json").read_text())
assert len(b)==480 and a[480:]==b
(p/"compare.json").write_text(json.dumps({"suffix_exact":True,"steps":len(b)}))
print("ARCH04_COMPARE", (p/"compare.json").read_text())
