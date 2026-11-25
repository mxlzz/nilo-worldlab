"""Read-only postmortem: compare retained restart/continuous arrays and metadata, no physics."""
import io,json,hashlib,pathlib,zipfile
import numpy as np
import xarray as xr

D=pathlib.Path("evidence/P1-C/C2b")
P=pathlib.Path("evidence/P1-C/C2a/14977708755.zip")
def sha(b):return hashlib.sha256(b).hexdigest()
def normal(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    return x
def nc(z,n):
    with xr.open_dataset(io.BytesIO(z.read(n)),engine="scipy") as d:return d.load()
def compare(a,b):
    out={}
    for name in sorted(set(a.variables)|set(b.variables)):
        if name not in a.variables or name not in b.variables:
            out[name]={"missing":True};continue
        x,y=a[name],b[name];u,v=np.asarray(x),np.asarray(y)
        m={}
        if x.dims!=y.dims:m["dims"]=[x.dims,y.dims]
        if u.dtype!=v.dtype:m["dtype"]=[str(u.dtype),str(v.dtype)]
        if x.attrs!=y.attrs:m["attrs"]=[x.attrs,y.attrs]
        if u.shape!=v.shape:m["shape"]=[u.shape,v.shape]
        elif u.dtype.kind in "fiu" and v.dtype.kind in "fiu":
            bits=u.astype(u.dtype.newbyteorder("<"),copy=False).tobytes()==v.astype(v.dtype.newbyteorder("<"),copy=False).tobytes()
            if not bits:
                m["different_values"]=int(np.count_nonzero(u!=v))
                m["max_abs"]=float(np.max(np.abs(u.astype(float)-v.astype(float))))
                m["different_bits"]=int(np.count_nonzero(u.view(np.uint8)!=v.view(np.uint8))) if u.flags.c_contiguous and v.flags.c_contiguous and u.dtype==v.dtype else None
                if u.dtype==v.dtype==np.float64 and np.all(u>=0) and np.all(v>=0):
                    uu=u.view(np.uint64);vv=v.view(np.uint64)
                    m["max_ulp"]=int(np.max(np.maximum(uu,vv)-np.minimum(uu,vv)))
        elif not np.array_equal(u,v):m["different_values"]=int(np.count_nonzero(u!=v))
        if m:out[name]=m
    attrs={k:[normal(a.attrs.get(k)),normal(b.attrs.get(k))] for k in sorted(set(a.attrs)|set(b.attrs)) if normal(a.attrs.get(k))!=normal(b.attrs.get(k))}
    return {"fields":out,"attrs":attrs}
def main():
    out={}
    with zipfile.ZipFile(P) as p,zipfile.ZipFile(D/"14977979956.zip") as q:
        a=json.loads(p.read("result.json"));b=json.loads(q.read("result.json"))
        for key,v in b["cases"].items():
            ref=a["cases"][key];rows=[]
            for epoch in (4,5,6):
                name=f"{key}/epoch_{epoch}.nc";x,y=nc(p,name),nc(q,name)
                d=compare(x,y)
                e=ref["epochs"][epoch];f=v["epochs"][epoch-4]
                d["logical_equal"]={k:e[k]==f[k] for k in f}
                d["file_bytes_equal"]=p.read(name)==q.read(name)
                for prefix in ("flux","ray"):
                    name=f"{key}/{prefix}_{epoch}.npz"
                    with np.load(io.BytesIO(p.read(name)),allow_pickle=False) as u,np.load(io.BytesIO(q.read(name)),allow_pickle=False) as w:
                        d[prefix+"_arrays_exact"]=set(u.files)==set(w.files) and all(u[n].dtype==w[n].dtype and u[n].shape==w[n].shape and u[n].tobytes()==w[n].tobytes() for n in u.files)
                rows.append({"epoch":epoch,**d})
            out[key]=rows
    result={"diagnostic_only":True,"cases":out}
    (D/"diagnostic.json").write_text(json.dumps(result,indent=2,sort_keys=True,default=normal)+"\n")
    print("WORLD_DIAGNOSTIC_JSON "+json.dumps(result,sort_keys=True,default=normal),flush=True)
if __name__=="__main__":main()
