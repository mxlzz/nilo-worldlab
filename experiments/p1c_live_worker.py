"""P1-C live SPACE worker, B.2a step audit unchanged, no archived future."""
import hashlib,json,os,pathlib,shutil,sys
import numpy as np
import xarray as xr
from numpy._core import _multiarray_umath as cpu
import probe_e6_h33_resume_landlab as h33
import probe_e5_h31_owner_contract as c
from p1b2a_native_flux import PROFILE

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def req(v,m):
    if not bool(v):raise ValueError(m)
def main(rawpath,seed,k,previous,target,folder):
    seed=int(seed);k=float(k);target=int(target);out=pathlib.Path(folder)
    req(all(os.environ.get(key)==val for key,val in PROFILE.items()),"profile")
    req(not any(cpu.__cpu_features__.get(key,False) for key in PROFILE["NPY_DISABLE_CPU_FEATURES"].split(",") if key in cpu.__cpu_dispatch__),"CPU dispatch")
    req(not out.exists() and target in range(0,19,3),"duplicate or invalid request")
    raw=h33.read_raw(rawpath)[seed]
    if target==0:
        req(previous=="-","unexpected predecessor")
        owner=h33.initial_owner(raw,seed,k);input_sha=None
    else:
        with xr.open_dataset(previous,engine="scipy") as ds:owner=ds.load()
        c.validate(owner,raw,seed,k,max_step=18)
        req(owner.attrs["landlab_step"]==target-3,"wrong predecessor")
        input_sha=c.hash_owner(owner)
    errors=[];steps=[];arrays={}
    def gate(label,ok):
        if not bool(ok):errors.append(label)
    if target:
        block=target-2
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

    req(not errors,"B.2a gates: "+repr(errors))
    c.validate(owner,raw,seed,k,max_step=18)
    stage=out.with_name(out.name+".staging")
    req(not stage.exists(),"staging exists")
    stage.mkdir()
    try:
        owner.to_netcdf(stage/"soil.nc",engine="scipy")
        with xr.open_dataset(stage/"soil.nc",engine="scipy") as ds:back=ds.load()
        req(c.hash_owner(back)==c.hash_owner(owner),"worker roundtrip")
        np.savez_compressed(stage/"flux.npz",**arrays)
        event={"seed":seed,"K":k,"from":target-2 if target else None,"to":target,"input_sha":input_sha,"output_sha":c.hash_owner(owner),"steps":steps,"files":{p.name:sha(p) for p in stage.iterdir()}}
        (stage/"event.json").write_text(json.dumps(event,sort_keys=True)+"\n")
        stage.replace(out)
        print("P1C_SOIL_JSON "+json.dumps(event,sort_keys=True),flush=True)
    finally:
        if stage.exists():shutil.rmtree(stage)
if __name__=="__main__":main(*sys.argv[1:])
