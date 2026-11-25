"""P6.1-L: candidate mass-flux draining limiter; no state accepted."""
import io,json,pathlib,hashlib,zipfile,time,sys
import numpy as np
import xarray as xr
from landlab.components import OverlandFlow
from p6_horizon import grid,req,DX,DT,A,FILM,zphys,rate,validate,total
from p4_energy import temp,capacity,transfer,CW,TR
from p3_water import digest,sha
from p6_dry_diagnostic import ARCH,ASHA
def limit(H,f,t,head,c,act):
    H,f,t,head,c,act=map(np.asarray,(H,f,t,head,c,act))
    req(np.isfinite(H).all() and np.isfinite(f).all() and np.all(H>=0) and np.all(H[~c]==0),"bounds")
    req(np.all(f[~act]==0),"inactive")
    src=np.where(f>=0,t,head);mag=np.abs(f)
    req(np.all(mag[~c[src]]==0),"undeclared inflow")
    G=np.zeros_like(H);np.add.at(G,src,mag)
    s=np.ones_like(H);m=G>H
    s[m]=np.nextafter(H[m]/G[m],0.)
    s[m & (H==0)]=0.
    return f*s[src],s,G
def fixture():
    H=np.array([1.,0.,0.]);f=np.array([.8,.7]);t=np.array([0,0]);hd=np.array([1,2]);c=np.ones(3,dtype=bool);act=np.ones(2,dtype=bool)
    a,s,G=limit(H,f,t,hd,c,act)
    req(np.allclose(a,[8/15,7/15],rtol=0,atol=1e-12) and a.sum()<=1.,"analytic")
    z,_,_=limit(H*0,f,t,hd,c,act);req(np.all(z==0),"dry")
    z,_,_=limit(H,f*0,t,hd,c,act);req(np.array_equal(z,f*0),"zero")
    try:limit(H,np.array([-.8,.7]),t,hd,np.array([True,False,True]),act)
    except ValueError as e:req(str(e)=="undeclared inflow","inflow reason")
    else:raise ValueError("inflow accepted")
    return {"analytic_flux_m3":a.tolist(),"dry":True,"zero":True,"reject_inflow":True}
def states():
    req(sha(ARCH)==ASHA,"archive hash")
    out={}
    with zipfile.ZipFile(ARCH) as z:
        req(z.testzip() is None,"CRC");rep=json.loads(z.read("reference/result.json"))
        for mode,step,member,key in [("periodic",240,"step0240.nc","240"),("periodic",2400,"step2400.nc","2400"),("relaxation",2078,"last_valid.nc","last_valid")]:
            name=f"{mode}/{member}";b=z.read("reference/"+name);req(hashlib.sha256(b).hexdigest()==rep["outputs"][name],"member hash")
            with xr.open_dataset(io.BytesIO(b),engine="scipy") as f:ds=f.load()
            req(digest(ds)==rep["cases"][mode]["snapshots"][key]["owner_sha"],"owner hash");validate(ds,step);out[f"{mode}_{step}"]=ds
    return out
def probe(ds):
    d=digest(ds);g,c,act,t,hd,b=grid()
    g.at_node["topographic__elevation"]=zphys(ds).copy();g.at_node["surface_water__depth"]=np.zeros(len(c))
    of=OverlandFlow(g,h_init=FILM,alpha=.7,theta=.8,mannings_n=.03,steep_slopes=False,rainfall_intensity=0.)
    h=g.at_node["surface_water__depth"];q=g.at_link["surface_water__discharge"];H=ds.H.values.copy()
    np.copyto(h,np.where(c,H/A,FILM));np.copyto(q,ds.q.values);np.copyto(of._h_links,ds.hl.values)
    old=h.copy();rr=rate(ds.attrs["thermal_mode"],int(ds.attrs["step"])+1);of.rainfall_intensity=rr
    req(DT<=of.calc_time_step() and of.run_one_step(DT)==DT,"native duration/CFL")
    proposed=q.copy();hh=h.copy();hl=of._h_links.copy();f=proposed*DX*DT
    accepted,s,G=limit(H,f,t,hd,c,act);changed=np.where(s<1)[0]
    qq=proposed.copy()
    src=np.where(f>=0,t,hd);affected=s[src]<1
    qq[affected]=accepted[affected]/(DX*DT)
    if len(changed):
        # Use accepted flux exactly once per interface; evaluate continuity in volume.
        ff=qq*DX*DT;mag=np.abs(ff);src2=np.where(ff>=0,t,hd);dst2=np.where(ff>=0,hd,t)
        gross=np.zeros_like(H);incoming=np.zeros_like(H)
        np.add.at(gross,src2,mag);np.add.at(incoming,dst2,mag)
        Hn=H-gross+incoming+np.where(c,rr*DT*A,0.);hn=np.where(c,Hn/A,FILM)
    else:
        hn=hh.copy();ff=f.copy()
    req(np.array_equal(qq[~affected],proposed[~affected]),"unaffected changed")
    req(np.all(hn>=0) and np.all(hn[~c]==FILM) and np.all(qq[~act]==0) and np.all(hl[act]>0),"candidate bounds")
    div=g.calc_flux_div_at_node(qq);local=float(np.max(np.abs((hn-old-(rr-div)*DT)[c])))
    bf=ff*b;out=float(bf[bf>0].sum());inc=float(-bf[bf<0].sum());res=float((hn-old)[c].sum()*A+out-inc-rr*DT*A*int(c.sum()))
    req(local<=1e-10 and abs(res)<=1e-8 and inc==0,"candidate continuity")
    V=total(ds);W=ds.W.values;E=ds.E_sensible_J.values;T=temp(E,capacity(V,W,H),c)
    ev,dv,exp,vout,ef,gross=transfer(H,T,ff,t,hd,c,act,CW)
    rv=np.where(c,rr*DT*A,0.);re=rv*CW*(300.-TR);En=E+re+ev
    er=abs(float(En.sum()+exp-re.sum()-E.sum()))
    req(np.all(gross<=H+1e-12) and np.max(np.abs(H+rv+dv-np.where(c,hn*A,0.)))<=1e-10 and abs(vout-out)<=1e-10 and er<=.001,"advection")
    req(digest(ds)==d,"owner mutation")
    rep={"input_owner_sha":d,"limited_nodes":changed.tolist(),"limited_links":int(affected.sum()),"native_min_h":float(hh[c].min()),"min_h":float(hn[c].min()),"max_local_residual_m":local,"global_residual_m3":res,"energy_residual_J":er,"max_overdraw_m3":float(max(0,np.max(gross-H))),"native_identity":bool(np.array_equal(hn,hh) and np.array_equal(qq,proposed)),"flux_removed_m3":float(np.abs(f-ff).sum()),"out_m3":out}
    arrays={"H":H,"old_h":old,"native_h":hh,"native_q":proposed,"h":hn,"q":qq,"hl":hl,"f":f,"accepted_f":ff,"s":s,"G":G,"gross":gross,"energy_delta":ev,"energy_link":ef,"core":c,"active":act,"tail":t,"head":hd}
    return rep,arrays
def main():
    start=time.perf_counter();d=pathlib.Path(sys.argv[1]);d.mkdir();rep={"experiment":"P6.1-L","fixture":fixture(),"cases":{},"failed_gates":[]}
    for name,ds in states().items():
        r,a=probe(ds);r2,a2=probe(ds);r["repeat_exact"]=r==r2 and all(np.array_equal(v,a2[k]) for k,v in a.items())
        if not r["repeat_exact"]:rep["failed_gates"].append(name+"_repeat")
        if name.startswith("periodic"):
            if r["limited_nodes"] or not r["native_identity"]:rep["failed_gates"].append(name+"_wet_identity")
        elif 17 not in r["limited_nodes"] or r["native_identity"]:rep["failed_gates"].append(name+"_dry")
        rep["cases"][name]=r;np.savez_compressed(d/(name+".npz"),**a)
    rep["wall_s"]=time.perf_counter()-start;rep["outputs"]={p.name:sha(p) for p in d.glob("*.npz")}
    (d/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_LIMITER_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    req(not rep["failed_gates"],"limiter")
if __name__=="__main__":main()
