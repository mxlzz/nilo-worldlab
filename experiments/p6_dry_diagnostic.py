"""P6.1-D: read-only one-interval dry-boundary diagnosis."""
import io,json,pathlib,hashlib,zipfile,time,sys
import numpy as np
import xarray as xr
from landlab.components import OverlandFlow
from p6_horizon import grid,req,DX,DT,A,FILM,zphys,rate,validate
from p3_water import digest,sha
ARCH=pathlib.Path("evidence/P6/P61/15031640270.zip")
ASHA="a37481e1a0b72a2cca2bb826cb5ae4d92f14c0bac2fff26819e6be8e47cc5c6c"
MEM="reference/relaxation/last_valid.nc"
MSHA="5466887bef58bc68b137ee2cf143468a9a5d1bfdf03668d4eb11a59f3aade2d1"
OWNER="772e0b24764ea7751868038f42ddc244f4d97d34300afc6824a8e8c6c5f64611"
def load():
    req(sha(ARCH)==ASHA,"archive hash")
    with zipfile.ZipFile(ARCH) as z:
        req(z.testzip() is None,"archive CRC")
        b=z.read(MEM);req(hashlib.sha256(b).hexdigest()==MSHA,"member hash")
    with xr.open_dataset(io.BytesIO(b),engine="scipy") as f:ds=f.load()
    req(digest(ds)==OWNER,"owner hash");validate(ds,2078)
    return ds
def probe(ds,n,steep):
    g,c,act,t,hd,b=grid()
    g.at_node["topographic__elevation"]=zphys(ds).copy()
    g.at_node["surface_water__depth"]=np.zeros(len(c))
    of=OverlandFlow(g,h_init=FILM,alpha=.7,theta=.8,mannings_n=.03,steep_slopes=steep,rainfall_intensity=0.)
    h=g.at_node["surface_water__depth"];q=g.at_link["surface_water__discharge"]
    np.copyto(h,np.where(c,ds.H.values/A,FILM));np.copyto(q,ds.q.values);np.copyto(of._h_links,ds.hl.values)
    req(np.array_equal(c,ds.core.values.astype(bool)) and np.array_equal(act,ds.active.values.astype(bool)),"topology")
    dt=DT/n;rr=rate("relaxation",2079);of.rainfall_intensity=rr
    initial=h.copy();rows=[];arrays={"initial_h":initial,"initial_q":q.copy(),"initial_hl":of._h_links.copy(),"z":zphys(ds),"core":c,"active":act,"tail":t,"head":hd,"boundary":b}
    elapsed=0.;flow=np.zeros_like(q);out=inc=0.;maxlocal=0.
    for j in range(n):
        old=h.copy();cfl=float(of.calc_time_step());req(dt<=cfl,"CFL")
        req(of.run_one_step(dt)==dt,"duration");elapsed+=dt
        hh=h.copy();qq=q.copy();hl=of._h_links.copy();div=g.calc_flux_div_at_node(qq)
        local=hh-old-(rr-div)*dt
        f=qq*DX*dt;flow+=f;bf=f*b;out+=float(bf[bf>0].sum());inc+=float(-bf[bf<0].sum())
        gross=np.zeros(len(c));incoming=np.zeros(len(c))
        np.add.at(gross,t,np.maximum(f,0));np.add.at(gross,hd,np.maximum(-f,0))
        np.add.at(incoming,hd,np.maximum(f,0));np.add.at(incoming,t,np.maximum(-f,0))
        guards={"nonnegative_h":bool(np.all(hh>=0)),"halo_film":bool(np.all(hh[~c]==FILM)),"inactive_q":bool(np.all(qq[~act]==0)),"active_hl":bool(np.all(hl[act]>0))}
        loc=float(np.max(np.abs(local[c])));maxlocal=max(maxlocal,loc)
        neg=np.where(hh<0)[0];over=np.where(c & (gross>old*A+1e-12))[0]
        row={"substep":j+1,"dt":dt,"cfl":cfl,"guards":guards,"negative_nodes":neg.tolist(),"negative_volume_m3":float(-hh[neg].sum()*A),"min_h_m":float(hh[c].min()),"local_residual_m":loc,"global_residual_m3":float(local[c].sum()*A),"gross_overdraw_nodes":over.tolist(),"max_gross_overdraw_m3":float(max(0,np.max((gross-old*A)[c]))),"negative_details":[{"node":int(i),"old_h_m":float(old[i]),"new_h_m":float(hh[i]),"gross_out_m3":float(gross[i]),"gross_in_m3":float(incoming[i])} for i in neg]}
        rows.append(row)
        arrays.update({f"{j+1}_{k}":v for k,v in {"old_h":old,"h":hh,"q":qq,"hl":hl,"div":div,"local":local,"gross":gross,"incoming":incoming}.items()})
        if not all(guards.values()):break
    residual=float((h-initial)[c].sum()*A+out-inc-rr*elapsed*A*int(c.sum()))
    arrays["integrated_flux_m3"]=flow
    return {"n":n,"steep_slopes":steep,"elapsed":elapsed,"rows":rows,"out_m3":out,"in_m3":inc,"global_residual_m3":residual,"max_local_residual_m":maxlocal,"eligible_hydraulic_candidate":bool(len(rows)==n and all(all(v["guards"].values()) for v in rows) and maxlocal<=1e-10 and abs(residual)<=1e-8),"gross_overdraw_free":all(not v["gross_overdraw_nodes"] for v in rows)},arrays
def main():
    start=time.perf_counter();d=pathlib.Path(sys.argv[1]);d.mkdir();ds=load();cases={};arrs={}
    for label,n,steep in [("original",1,False),("repeat",1,False),("half",2,False),("quarter",4,False),("steep",1,True)]:
        row,arrays=probe(ds,n,steep);cases[label]=row;arrs[label]=arrays;np.savez_compressed(d/(label+".npz"),**arrays)
    exact=cases["original"]==cases["repeat"] and all(np.array_equal(v,arrs["repeat"][k]) for k,v in arrs["original"].items())
    fail=[]
    if not exact:fail.append("same_host_repeat")
    if all(cases["original"]["rows"][0]["guards"].values()):fail.append("original_guard_not_reproduced")
    if digest(ds)!=OWNER:fail.append("owner_mutation")
    rep={"experiment":"P6.1-D","input_owner_sha":OWNER,"repeat_exact":exact,"cases":cases,"failed_diagnostic_gates":fail,"wall_s":time.perf_counter()-start}
    rep["outputs"]={p.name:sha(p) for p in d.glob("*.npz")}
    (d/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_DRY_DIAG_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    req(not fail,"diagnostic")
if __name__=="__main__":main()
