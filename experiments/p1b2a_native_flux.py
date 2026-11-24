"""P1-B.2a: observe native SPACE gross fluxes; do not modify inherited physics."""
import base64, gzip, hashlib, json, os, pathlib, platform, sys, time, traceback, zipfile
from importlib.metadata import version
import numpy as np
import xarray as xr
from numpy._core import _multiarray_umath as cpu
import probe_e6_h33_resume_landlab as h33
import probe_e5_h31_owner_contract as c

ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/"p1b2a-output"
ARCHIVES={
 "raw":("45389348480/14917807888.zip","f2b204f42f7d055600587230093ad36b1a89f0f85ae758bae114a4b6b94a5170"),
 "continuous":("45447643891/14946621830.zip","34849933ca9b1c4fd1b2b1b5c9538f4093caef62950a71859f8583a58c7d9c44"),
 "results":("45447643891/14948035431.zip","0cb082e656d34aef4d2b5387d1e4d2cba81b0dd9337f100dde4a715956db5fad"),
 "thermal":("45488606315/14965401716.zip","7f919ae46e68c8da7d775f297f00a662c98f37752f7b78658b419eddcf79d512")}
EXTRA={"grain_state","grain_mass_kg","grain_radius_m_nominal","source_id","release_step","native_id","grain_reservoir_kg"}
ATTRS={"p0d_seed","p0d_formation_steps","p0d_bullet_steps","p0d_epoch","p0d_plan_sha","p0d_geometry"}
PROFILE={"NPY_DISABLE_CPU_FEATURES":"X86_V3,X86_V4,AVX512_ICL,AVX512_SPR","OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1","MKL_NUM_THREADS":"1"}

def read_inputs():
    out={}
    for name,(relative,sha) in ARCHIVES.items():
        p=ROOT/"evidence/archives"/relative
        assert hashlib.sha256(p.read_bytes()).hexdigest()==sha, name+" digest"
        dest=OUT/"input"/name
        dest.mkdir(parents=True)
        with zipfile.ZipFile(p) as z:
            assert z.testzip() is None
            for n in z.namelist():
                assert not pathlib.PurePosixPath(n).is_absolute() and ".." not in pathlib.PurePosixPath(n).parts
            z.extractall(dest)
        out[name]=dest
    return out

def one(root,name):
    paths=list(root.rglob(name))
    assert len(paths)==1,(root,name,paths)
    return paths[0]

def reference(inputs,raw,seed,k,epoch,report):
    key=f"{seed}_{k}"
    path=one(inputs["continuous"],key+f"/epoch_{epoch}.nc")
    with xr.open_dataset(path,engine="scipy") as ds:
        soil=ds.load().drop_vars(list(EXTRA)).drop_dims(["grain","motion"])
    soil.attrs={k:v for k,v in soil.attrs.items() if k not in ATTRS}
    c.validate(soil,raw,seed,k,max_step=18)
    expected=report["cases"][key]["epochs"][epoch]["soil_sha"]
    assert c.hash_owner(soil)==expected, "archived soil reference mismatch"
    return expected

def compare_case(inputs,raw,seed,k,report):
    owner=h33.initial_owner(raw,seed,k)
    key=f"{seed}_{k}"
    errors=[]; epochs=[]; steps=[]; arrays={}
    def gate(label,ok):
        if not bool(ok):errors.append(label)
    def check_epoch(epoch):
        expected=reference(inputs,raw,seed,k,epoch,report)
        actual=c.hash_owner(owner)
        epochs.append({"epoch":epoch,"expected":expected,"actual":actual,"exact":actual==expected})
        gate(f"epoch_{epoch}_hash",actual==expected)
    check_epoch(0)
    for block in range(1,19,3):
        grid,flow,solver,outlet=h33.model_from_owner(owner,k)
        for step in range(block,block+3):
            before=owner.soil_mass_kg.values.ravel().copy()
            np.copyto(grid.at_node["soil__depth"],before/(c.RHO*c.AREA))
            np.copyto(grid.at_node["topographic__elevation"],c.topography(owner).ravel())
            flow.run_one_step()
            suspended,bedload=solver.run_one_step(c.LANDLAB_DT)
            dt=c.LANDLAB_DT; area=grid.cell_area_at_node
            node=grid.at_node
            e=node["sediment__erosion_flux"].copy()*area*dt*c.RHO
            d=node["sediment__deposition_flux"].copy()*area*dt*c.RHO
            qin=node["sediment__influx"].copy()*dt*c.RHO
            qout=node["sediment__outflux"].copy()*dt*c.RHO
            br=node["bedrock__erosion_flux"].copy()*area*dt*c.RHO
            rec=node["flow__receiver_node"].copy()
            order=node["flow__upstream_node_order"][::-1].copy()
            active=(grid.status_at_node==0)&(node["surface_water__discharge"]>0)
            core=owner.is_core_node.values.ravel()!=0
            after=grid.at_node["soil__depth"].copy()*c.RHO*c.AREA
            export=float(suspended+bedload)*c.RHO
            candidate=owner.copy(deep=True)
            candidate.soil_mass_kg.values[:]=after.reshape(c.SHAPE)
            candidate.attrs["exported_soil_kg"]=float(owner.attrs["exported_soil_kg"])+export
            candidate.attrs["landlab_step"]=step
            candidate.attrs["landlab_t_bench"]=step*c.LANDLAB_DT
            check=c.validate(candidate,raw,seed,k,max_step=18)
            h33.h29.validate(c.soil_only(candidate),raw,c.SOURCES[seed],seed)
            tol=np.maximum(1e-6,1e-9*np.maximum.reduce([before,after,e,d,qin,qout,np.ones_like(e)]))
            mass=after-before+e-d
            balance=qin+e-qout-d
            incoming=np.zeros_like(qin)
            pos=np.empty(len(order),dtype=int);pos[order]=np.arange(len(order))
            valid=np.all((rec>=0)&(rec<len(rec)))
            edges=np.flatnonzero(qout>0)
            gate(f"{step}:receiver",valid and all(rec[i]!=i and pos[i]<pos[rec[i]] for i in edges))
            if valid:np.add.at(incoming,rec,qout)
            route=incoming-qin
            boundary_export=float(qin[~core].sum())
            for name,condition in {
                "finite":np.isfinite(np.stack([before,after,e,d,qin,qout,br])).all(),
                "nonnegative":all(np.all(a>=0) for a in [before,after,e,d,qin,qout,br]),
                "mass":np.all(np.abs(mass)<=tol),
                "flow":np.all(np.abs(balance[core])<=tol[core]),
                "route":np.all(np.abs(route)<=tol),
                "export":abs(boundary_export-export)<=max(1e-6,1e-9*max(export,1)),
                "no_rock_wash":float(np.abs(br).sum())+abs(float(suspended)*c.RHO)<=1e-6,
                "erosion_available":np.all(e<=before),
                "inactive_zero_out":np.all(qout[~active]==0),
                "global":abs(check["mass_residual_kg"])<=.09216,
                "K0":k!=0 or (np.array_equal(before,after) and export==0 and not any(np.any(a) for a in [e,d,qin,qout])),
            }.items():gate(f"{step}:{name}",condition)
            steps.append({"step":step,"mass_abs_max_kg":float(np.abs(mass).max()),"flow_abs_max_kg":float(np.abs(balance[core]).max()),"route_abs_max_kg":float(np.abs(route).max()),"export_kg":export,"boundary_export_kg":boundary_export,"export_residual_kg":boundary_export-export,"global_residual_kg":float(check["mass_residual_kg"]),"erosion_kg":float(e.sum()),"deposition_kg":float(d.sum()),"inactive_out_kg":float(qout[~active].sum()),"max_erosion_minus_stock_kg":float((e-before).max())})
            for n,a in {"before":before,"after":after,"erosion":e,"deposition":d,"influx":qin,"outflux":qout,"bedrock":br,"receiver":rec,"order":order,"active":active,"core":core}.items():
                arrays[f"s{step}_{n}"]=a
            owner=candidate
        # Mirror original worker serialization at each transition.
        p=OUT/f"{key}_soil.nc";owner.to_netcdf(p,engine="scipy")
        with xr.open_dataset(p,engine="scipy") as ds: owner=ds.load()
        p.unlink()
        check_epoch((block+2)//3)
    if seed==7321 and k==.1:
        control=h33.initial_owner(raw,seed,k)
        for block in range(1,19,3):
            control,_=h33.advance(control,raw,seed,k,block,block+2,max_step=18)
            p=OUT/"control.nc";control.to_netcdf(p,engine="scipy")
            with xr.open_dataset(p,engine="scipy") as ds:control=ds.load()
            p.unlink()
            gate(f"uninstrumented_epoch_{(block+2)//3}",c.hash_owner(control)==epochs[(block+2)//3]["actual"])
    np.savez_compressed(OUT/f"{key}_fluxes.npz",**arrays)
    return {"epochs":epochs,"steps":steps,"failed":errors}

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    report={"experiment":"P1-B.2a","cases":{},"failed_preregistered_gates":[],"runtime":{},"input_archives":ARCHIVES}
    try:
        assert all(os.environ.get(k)==v for k,v in PROFILE.items())
        assert not any(cpu.__cpu_features__.get(k,False) for k in PROFILE["NPY_DISABLE_CPU_FEATURES"].split(",") if k in cpu.__cpu_dispatch__)
        report["runtime"]={"python":platform.python_version(),"versions":{p:version(p) for p in ("landlab","numpy","scipy","xarray")},"cpu":platform.processor(),"cpu_features":cpu.__cpu_features__,"profile":PROFILE}
        inputs=read_inputs()
        refs=json.loads(one(inputs["continuous"],"result.json").read_text())
        assert refs["experiment"]=="P0-D" and refs["mode"]=="produce" and not refs["failed_preregistered_gates"]
        raws=h33.read_raw(one(inputs["raw"],"worldengine_two_seeds.npz"))
        for seed in (7321,8128):
            for k in (0.,.1):
                key=f"{seed}_{k}"
                try:
                    case=compare_case(inputs,raws[seed],seed,k,refs)
                    report["cases"][key]=case
                    report["failed_preregistered_gates"] += [key+":"+x for x in case["failed"]]
                except Exception as e:
                    traceback.print_exc()
                    report["cases"][key]={"error":repr(e)}
                    report["failed_preregistered_gates"].append(key+":instrument:"+repr(e))
    except Exception as e:
        traceback.print_exc()
        report["failed_preregistered_gates"].append("instrument:"+repr(e))
    report["elapsed_s"]=time.monotonic()-start
    report["outputs"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob("*.npz")}
    (OUT/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    # Durable readback of metrics through CI logs; raw flux NPZ uploaded separately.
    print("WORLD_EVAL_JSON "+json.dumps(report,sort_keys=True),flush=True)
    if report["failed_preregistered_gates"]:sys.exit(1)

if __name__=="__main__":main()
