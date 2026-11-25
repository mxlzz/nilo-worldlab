"""P6.1: same P5.2 operators, extended horizon and explicit repeated forcing."""
import hashlib,json,pathlib,sys,time,resource,traceback,zipfile,io
import numpy as np
import xarray as xr
from p3_reciprocal import Model as Hydro,initial as hydro_initial,zphys as hydro_zphys,rain,load_provider,DX,DT,A,PHI,RHO_S,FILM,MID,STEPS,LEDGER as HYDRO_LEDGER
from p3_exner import grid,req
from p3_pore import remap,ga,geometry,domain,input_hash
from p3_water import digest,sha
from p4_energy import capacity,temp,transfer,total_balance,TR,CS,CW,SIGMA,LEDGER as ENERGY_LEDGER
from p5_grain import move
MODES=("periodic","relaxation")
STEPS=2400;MID=1200;CYCLE=240
ARCH=pathlib.Path("evidence/P5/P52/15030611016.zip")
ARCH_SHA="d9716143daad2d8d331dbdbc3526ea88294f677da307c9f6f116733c9f16fea1"
def phase(k):return (k-1)%CYCLE+1
def rate(mode,k):return 0. if mode=="relaxation" and k>CYCLE else (72. if phase(k)<=60 else 36. if 120<phase(k)<=180 else 0.)/3.6e6
def sun(mode,k):return not (mode=="relaxation" and k>CYCLE) and phase(k)<=120
def forcing_total(mode,k):
    n=min(k,CYCLE) if mode=="relaxation" else k
    full,rem=divmod(n,CYCLE)
    rain_m=(full*(72.*60+36.*60)+72.*min(rem,60)+36.*max(0,min(rem-120,60)))/3.6e6*DT*A*25
    sun_n=full*120+min(rem,120)
    return rain_m,sun_n*.8*600*A*DT*25
L=2.45e6;CPW=4180.;RHOW=1000.;CE=.0013;AIR=295.;RH=.5;PRESSURE=101.325
VAPOR=("vapor_kg","vapor_latent_J","vapor_sensible_J")
LEDGERS=HYDRO_LEDGER+ENERGY_LEDGER+VAPOR
DYNAMIC=("V_class","W","H","q","hl","E_sensible_J")
from p4_evap import evaporate

def frac_initial(ds):
    return np.where(ds.core.values.astype(bool),np.where(ds.x.values<=6,.25,.75),.5)
def total(ds):return ds.V_class.values.sum(axis=1)
def view(ds):
    v=ds.copy(deep=True);v["V"]=("node",total(ds));return v
def zphys(ds):return hydro_zphys(view(ds))
def fraction(ds):
    V=total(ds);c=ds.core.values.astype(bool);req(np.all(V[c]>0),"depleted column GSD domain")
    f=np.full(len(c),.5);f[c]=ds.V_class.values[c,1]/V[c];return f

def proposal(RBD,ds,h,hl,u,prev,fine):
    z=zphys(ds);g,c,act,t,hd,b=grid()
    g.at_node["topographic__elevation"]=z.copy()
    g.at_node["surface_water__depth"]=h.copy();g.at_link["surface_water__depth"]=hl.copy();g.at_link["surface_water__velocity"]=u.copy()
    gsd=np.vstack((np.r_[4.,np.full(len(c),100.)],np.r_[2.,100.*fine],np.r_[1.,np.zeros(len(c))]))
    r=RBD(g,gsd=gsd,bed_surf__gsd_loc_node=np.arange(len(c)),rho=1000,rho_s=RHO_S,bedload_equation="WilcockAndCrowe",lambda_p=PHI,variable_critical_shear_stress=False,outlet_boundary_condition="fixedValue",dt=DT,surface_water__velocity_prev_time_link=prev.copy(),track_stratigraphy=False)
    f=np.column_stack((1.-fine,fine));req(np.allclose(r._bed_surf__gsd_node,f,rtol=0,atol=1e-12),"GSD mapping")
    req(np.allclose(r._bed_surf__gsd_link,.5*(f[t]+f[hd]),rtol=0,atol=1e-12),"link GSD")
    req(np.array_equal(r._surface_water__velocity_prev_time_link,prev),"previous velocity")
    with np.errstate(all="warn"):
        r.shear_stress();r.bedload_equation();r.calculate_net_bedload();r.update_bed_elevation()
    qb=np.asarray(r._sed_transp__bedload_rate_link).copy();p=np.asarray(r._sed_transp__bedload_gsd_link).copy();nativez=g.at_node["topographic__elevation"].copy()
    req(np.isfinite(qb).all() and np.isfinite(p).all() and np.all(p>=0) and np.all(np.abs(p[qb!=0].sum(axis=1)-1)<=1e-12),"native fractions")
    for x,y in ((g.at_node["surface_water__depth"],h),(g.at_link["surface_water__depth"],hl),(g.at_link["surface_water__velocity"],u),(r._surface_water__velocity_prev_time_link,prev)):req(np.array_equal(x,y),"RBD mutation")
    req(np.all(qb[~act]==0) and np.all(qb[b!=0]==0) and np.array_equal(nativez[~c],z[~c]),"sediment boundary")
    return qb,p,nativez

def initial(mode):
    ds=hydro_initial("active");E=capacity(ds.V.values,ds.W.values,ds.H.values)*(300.-TR)
    f=frac_initial(ds);M=ds.V.values[:,None]*np.column_stack((1.-f,f))
    ds=ds.drop_vars("V").assign_coords(grain=np.array([0,1],dtype=np.int32));ds["V_class"]=(("node","grain"),M,{"units":"m3 solid","grain_intervals_mm":"2-4;1-2"})
    ds["E_sensible_J"]=("node",E,{"units":"J"})
    ds.attrs.update(thermal_schema="P6.1-v1",thermal_mode=mode,thermal_initial_J=float(E.sum()),temperature_reference_K=TR,latent_J_kg=L,CE=CE,air_T_K=AIR,air_RH=RH,**{n:0. for n in ENERGY_LEDGER+VAPOR})
    return ds
def validate(ds,expected=None):
    a=ds.attrs;mode=a["thermal_mode"];req(mode in MODES,"mode");ref=initial(mode)
    req(set(ds.variables)==set(ref.variables) and set(a)==set(ref.attrs),"schema/duplicate")
    req(all(np.array_equal(ds[k].values,ref[k].values) for k in ref.variables if k not in DYNAMIC),"fixed geometry")
    req(all(ds[k].dims==ref[k].dims and ds[k].shape==ref[k].shape and ds[k].attrs==ref[k].attrs for k in ref.variables),"dimensions/units")
    req(all(a[n]==ref.attrs[n] for n in ref.attrs if n not in LEDGERS+("step","time_s")),"metadata")
    c=ds.core.values.astype(bool);act=ds.active.values.astype(bool)
    M,W,H,q,hl,E=(ds[n].values for n in DYNAMIC);V=M.sum(axis=1)
    req(M.shape==(len(c),2) and np.isfinite(M).all() and np.all(M>=0) and np.all(M[~c]==0),"class stocks")
    req(np.max(np.abs(M.sum(axis=0)-ref.V_class.values.sum(axis=0)))<=1e-9,"class balance")
    req(all(np.isfinite(v).all() for v in (V,W,H,q,hl,E)) and all(np.all(v>=0) and np.all(v[~c]==0) for v in (V,W,H)),"stock bounds")
    req(np.all(W<=geometry(V,A)[1]+1e-12) and domain(V,W,A,c)[0],"GA/capacity")
    req(np.all(q[~act]==0) and np.all(hl[act]>0) and np.all(hl[~act]==FILM),"hydro cache")
    T=temp(E,capacity(V,W,H),c)
    k=int(a["step"]);req(a["step"]==k and 0<=k<=STEPS and a["time_s"]==k*DT and (expected is None or k==expected),"clock")
    req(all(np.isfinite(a[n]) and a[n]>=0 for n in LEDGERS),"ledger")
    rr,sw_expected=forcing_total(mode,k)
    req(abs(a["rain"]-rr)<=1e-10 and a["solid_out"]==a["solid_out_J"]==a["inc"]==0,"boundary/forcing")
    req(abs(a["rain_J"]-rr*CW*(300.-TR))<=.001,"rain energy")
    req(abs(a["shortwave_J"]-sw_expected)<=1e-6,"pulse")
    req(abs(a["vapor_latent_J"]-a["vapor_kg"]*L)<=.001,"latent ledger")
    sr=float(V.sum()+a["solid_out"]-total(ref).sum());wr=float((W+H).sum()+a["out"]+a["vapor_kg"]/RHOW-a["inc"]-a["rain"]-(ref.W.values+ref.H.values).sum())
    er=total_balance(E,a,a["thermal_initial_J"])+a["vapor_latent_J"]+a["vapor_sensible_J"]
    req(abs(sr)<=1e-9 and abs(wr)<=1e-8 and abs(er)<=.01,"global material/energy budget")
    return sr,wr,er
def snap(ds,p):
    validate(ds);d=digest(ds);ds.to_netcdf(p,engine="scipy")
    with xr.open_dataset(p,engine="scipy") as f:b=f.load()
    req(digest(b)==d,"NC roundtrip");validate(b)
    return {"sha256":sha(p),"owner_sha":d}
def read(p,meta=None,step=None):
    if meta:req(sha(p)==meta["sha256"],"file digest")
    with xr.open_dataset(p,engine="scipy") as f:ds=f.load()
    validate(ds,step)
    if meta:req(digest(ds)==meta["owner_sha"],"owner digest")
    return ds
def observe(ds):
    d=digest(ds);v=(float(ds.V_class.sum()),float(ds.H.sum()),float(ds.E_sensible_J.sum()))
    req(digest(ds)==d,"observation mutation");return v
class Model(Hydro):
    def __init__(self,ds,RBD):
        super().__init__(view(ds),RBD);self.ds=ds.copy(deep=True)
    def flow(self,k,freeze=False):
        ds=self.ds;c=self.c;g=self.g
        req(k==ds.attrs["step"]+1,"step")
        z=ds.z0.values if freeze else zphys(ds)
        np.copyto(g.at_node["topographic__elevation"],z)
        np.copyto(self.h,np.where(c,ds.H.values/A,FILM));req(np.array_equal(self.q,ds.q.values) and np.array_equal(self.of._h_links,ds.hl.values),"hydraulic cache")
        old=self.h.copy();prev=np.zeros_like(self.q);prev[self.act]=self.q[self.act]/ds.hl.values[self.act]
        rr=rate(ds.attrs["thermal_mode"],k);self.of.rainfall_intensity=rr
        req(DT<=self.of.calc_time_step(),"CFL");req(self.of.run_one_step(DT)==DT,"duration")
        q=self.q.copy();hl=np.asarray(self.of._h_links).copy();hh=self.h.copy()
        div=g.calc_flux_div_at_node(q);local=float(np.max(np.abs((hh-old-(rr-div)*DT)[c])))
        bf=q*self.b*DX*DT;out=float(bf[bf>0].sum());inc=float(-bf[bf<0].sum())
        req(local<=1e-10 and abs(float(div[c].sum()*A*DT)-(out-inc))<=1e-10,"hydraulic local")
        req(np.all(hh>=0) and np.all(hh[~c]==FILM) and np.all(q[~self.act]==0) and np.all(hl[self.act]>0),"hydraulic bounds")
        u=np.zeros_like(q);u[self.act]=q[self.act]/hl[self.act]
        req(np.isfinite(u).all(),"velocity")
        return hh,q,hl,u,prev,rr,out,inc,local
    def advance(self,k):
        ds=self.ds;c=self.c;mode=ds.attrs["thermal_mode"]
        M=ds.V_class.values;V=total(ds);W,H,E=(ds[n].values for n in ("W","H","E_sensible_J"))
        T=temp(E,capacity(V,W,H),c)
        h,q,hl,u,prev,rr,out,inc,flowlocal=self.flow(k,False)
        req(inc==0,"undeclared liquid inflow")
        Hflow=np.where(c,h*A,0.);rv=np.where(c,rr*DT*A,0.);re=rv*CW*(300.-TR)
        ew,dw,ewout,wout,linkw,_=transfer(H,T,q*DX*DT,self.t,self.hd,c,self.act,CW)
        req(np.max(np.abs(H+rv+dw-Hflow))<=1e-10 and abs(wout-out)<=1e-10,"water material transfer")
        Eflow=E+re+ew;op1=abs(float(Eflow.sum()+ewout-re.sum()-E.sum()))
        Tf=temp(Eflow,capacity(V,W,Hflow),c)
        z=zphys(ds);fine=fraction(ds);used=fine
        qb,p,nativez=proposal(self.RBD,ds,h,hl,u,prev,used)
        fs=qb*DX*DT;frac=np.column_stack((1.-fine,fine))
        chosen=p
        Mn,dc,amount,gross,classres=move(M,fs,chosen,self.t,self.hd,c,self.act)
        fs=qb*DX*DT;es,dv,esout,sout,links,_=transfer(V,Tf,fs,self.t,self.hd,c,self.act,CS)
        rem,a=remap(V,W,Hflow,fs,self.t,self.hd,c,self.act,A)
        Vn=Mn.sum(axis=1);req(np.max(np.abs(Vn-a["V_after"]))<=1e-10,"aggregate class/remap")
        dV=Vn-V;newz=np.where(c,ds.bedrock.values+Vn/((1-PHI)*A),ds.z0.values)
        exlocal=float(np.max(np.abs(((nativez-z)*(1-PHI)*A-dV)[c])))
        req(exlocal<=1e-10 and np.max(np.abs(newz-nativez))<=1e-12 and np.max(np.abs(V+dv-Vn))<=1e-12,"Exner/geometry")
        req(rem["GA_domain_valid"] and rem["solid_out_m3"]==sout==esout==0,"pore/solid boundary")
        gr,ga_a=ga(a,h,1e-6);Wn=ga_a["GA_W_after"];Hp=ga_a["GA_H_after"]
        Eadv=Eflow+es;op2=abs(float(Eadv.sum()+esout-Eflow.sum()));Cp=capacity(Vn,Wn,Hp);Ta=temp(Eadv,Cp,c)
        sw=np.zeros_like(E);li=sw.copy();lo=sw.copy();si=sw.copy();so=sw.copy()
        if sun(mode,k):sw[c]=.8*600.*A*DT
        li[c]=.95*SIGMA*285.**4*A*DT;lo[c]=.95*SIGMA*Ta[c]**4*A*DT
        sensible=10.*(AIR-Ta[c])*A*DT;si[c]=np.maximum(sensible,0.);so[c]=np.maximum(-sensible,0.)
        Ep=Eadv+sw+li-lo+si-so;op3=abs(float(Ep.sum()-Eadv.sum()-(sw+li-lo+si-so).sum()))
        Hn,En,ev=evaporate(Vn,Wn,Hp,Ep,c,ds.attrs["CE"])
        req(max(op1,op2,op3,ev["energy_error"])<=.001,"operator energy")
        nxt=ds.copy(deep=True)
        for name,val in (("V_class",Mn),("W",Wn),("H",Hn),("q",q),("hl",hl),("E_sensible_J",En)):np.copyto(nxt[name].values,val)
        at=nxt.attrs;at["step"]=k;at["time_s"]=k*DT
        add={"rain":float(rv.sum()),"out":out,"inc":inc,"solid_out":sout,"erosion":float(-dV[dV<0].sum()),"deposition":float(dV[dV>0].sum()),"release":rem["released_water_m3"],"infiltration":gr["infiltration_m3"],"rain_J":float(re.sum()),"water_out_J":ewout,"solid_out_J":esout,"shortwave_J":float(sw.sum()),"longwave_in_J":float(li.sum()),"longwave_out_J":float(lo.sum()),"sensible_in_J":float(si.sum()),"sensible_out_J":float(so.sum()),"vapor_kg":float(ev["m"].sum()),"vapor_latent_J":float(ev["latent"].sum()),"vapor_sensible_J":float(ev["sensible"].sum())}
        for name,val in add.items():at[name]+=val
        sr,wr,er=validate(nxt);self.ds=nxt
        Cn=capacity(Vn,Wn,Hn);Tn=temp(En,Cn,c)
        metrics={"step":k,"class_res":float(np.max(np.abs(Mn.sum(axis=0)-initial(mode).V_class.values.sum(axis=0)))),"owner_sha":digest(nxt),**{n:at[n] for n in LEDGERS},"solid_res":sr,"water_res":wr,"energy_res":er,"flow_local":flowlocal,"exner_local":exlocal,"pore_local":rem["max_local_water_residual_m3"],"ga_local":gr["max_local_water_residual_m3"],"evap_local_K":ev["temperature_error"],"max_energy_op":max(op1,op2,op3,ev["energy_error"]),"T_mean_K":float(TR+En.sum()/Cn.sum()),"T_min_K":float(Tn[c].min()),"T_max_K":float(Tn[c].max()),"min_H_m3":float(Hn[c].min()),"min_class_m3":float(Mn[c].min()),"min_F":rem["min_F_equiv_m"],"max_dz":float(np.max(np.abs(newz-ds.z0.values)))}
        trace={"V":Vn.copy(),"V_class":Mn.copy(),"p":p.copy(),"chosen_p":chosen.copy(),"W":Wn.copy(),"H":Hn.copy(),"q":q,"hl":hl,"flow_h":h,"qb":qb,"z":newz,"E":En.copy(),"T":Tn,"H_before_evap":Hp.copy(),"E_before_evap":Ep.copy(),"vapor_mass":ev["m"],"vapor_demand":ev["demand"],"vapor_latent":ev["latent"],"vapor_sensible":ev["sensible"],"water_link_energy":linkw,"solid_link_energy":links}
        return metrics,trace
def series(rows,p):
    if rows:np.savez_compressed(p,**{n:np.asarray([r[n] for r in rows]) for n in rows[0]})
def physical_fields(ds):return {n:ds[n].values.copy() for n in DYNAMIC}
def reference():
    req(sha(ARCH)==ARCH_SHA,"P5 archive SHA")
    with zipfile.ZipFile(ARCH) as z:
        req(z.testzip() is None,"P5 CRC");r=json.loads(z.read("reference/result.json"));req(not r["failed_preregistered_gates"],"P5 gate")
        name="active/final.nc";b=z.read("reference/"+name);req(hashlib.sha256(b).hexdigest()==r["outputs"][name],"P5 snapshot SHA")
        with xr.open_dataset(io.BytesIO(b),engine="scipy") as d:ds=d.load()
        req(digest(ds)==r["cases"]["active"]["final"]["owner_sha"],"P5 logical SHA")
        return ds
def producer(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir();start=time.monotonic()
    rep={"experiment":"P6.1","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"];finals={};prefix={}
    RBD=load_provider(out);ref=reference()
    for mode in MODES:
        folder=out/mode;folder.mkdir();rows=[];snaps={};m=None;k=0
        try:
            m=Model(initial(mode),RBD);snaps["0"]=snap(m.ds,folder/"step0000.nc")
            for k in range(1,STEPS+1):
                row,_=m.advance(k);rows.append(row)
                if k%CYCLE==0:
                    sm=snap(m.ds,folder/f"step{k:04d}.nc");snaps[str(k)]=sm;observe(m.ds)
                    if k==240:prefix[mode]=m.ds.copy(deep=True)
                    if k==MID and mode=="periodic":
                        meta=snap(m.ds,cp/"owner.nc");req(meta["owner_sha"]==sm["owner_sha"],"checkpoint")
                        (cp/"manifest.json").write_text(json.dumps({"step":MID,"mode":mode,**meta},sort_keys=True)+"\n")
            finals[mode]=m.ds.copy(deep=True)
        except Exception as e:
            traceback.print_exc();fail.append(mode+":"+repr(e))
            if m is not None:
                last=int(m.ds.attrs["step"]);snaps["last_valid"]=snap(m.ds,folder/"last_valid.nc")
                rep.setdefault("stops",{})[mode]={"attempted_step":k,"last_valid_step":last,"error":repr(e)}
        series(rows,folder/"series.npz")
        rep["cases"][mode]={"steps":len(rows),"snapshots":snaps,"last":rows[-1] if rows else None,"max_residuals":{n:max(abs(r[n]) for r in rows) for n in ("class_res","solid_res","water_res","energy_res","max_energy_op")} if rows else {},"min_H":min((r["min_H_m3"] for r in rows),default=None),"min_class":min((r["min_class_m3"] for r in rows),default=None),"min_F":min((r["min_F"] for r in rows),default=None),"T_range":[min(r["T_min_K"] for r in rows),max(r["T_max_K"] for r in rows)] if rows else None}
    if len(prefix)==2:
        a,b=(prefix[n] for n in MODES)
        exact=all(np.array_equal(a[n].values,b[n].values) for n in DYNAMIC) and all(a.attrs[n]==b.attrs[n] for n in LEDGERS)
        diffs={n:float(np.max(np.abs(a[n].values-ref[n].values))) for n in DYNAMIC}
        tol={"V_class":1e-10,"W":1e-10,"H":1e-10,"q":1e-9,"hl":1e-10,"E_sensible_J":.001}
        rep["prefix"]={"exact_between_modes":exact,"P52_max_diffs":diffs}
        if not exact or not all(diffs[n]<=tol[n] for n in diffs):fail.append("prefix")
    else:fail.append("prefix_incomplete")
    if len(finals)==2:
        a,b=(finals[n] for n in MODES)
        v={"H_diff_m3":float(np.max(np.abs(a.H.values-b.H.values))),"E_diff_J":float(np.max(np.abs(a.E_sensible_J.values-b.E_sensible_J.values))),"fine_change":float(np.max(np.abs(fraction(a)-frac_initial(a)))),"max_dz_m":float(np.max(np.abs(zphys(a)-a.z0.values)))}
        rep["comparisons"]=v
        if not (v["H_diff_m3"]>1e-6 or v["E_diff_J"]>100):fail.append("forcing_effect")
        if not (v["fine_change"]>1e-5 and v["max_dz_m"]>1e-4 and a.attrs["vapor_kg"]>0 and a.attrs["out"]>0):fail.append("evolution")
        for mode,ds in finals.items():
            rr,sw=forcing_total(mode,STEPS)
            if abs(ds.attrs["rain"]-rr)>1e-10 or abs(ds.attrs["shortwave_J"]-sw)>1e-6:fail.append("forcing_ledger_"+mode)
    else:fail.append("incomplete horizon")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_HORIZON_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
def consumer(checkpoint,dest):
    cp=pathlib.Path(checkpoint);out=pathlib.Path(dest);out.mkdir();meta=json.loads((cp/"manifest.json").read_text())
    req(meta["step"]==MID and meta["mode"]=="periodic","checkpoint metadata");ds=read(cp/"owner.nc",meta,MID)
    controls={}
    def reject(name,fn):
        try:fn()
        except ValueError:controls[name]=True
        else:raise ValueError("accepted "+name)
    bad=out/"bad.nc";b=bytearray((cp/"owner.nc").read_bytes());b[-1]^=1;bad.write_bytes(b);reject("file",lambda:read(bad,meta,MID));bad.unlink()
    def verify(x):
        validate(x,MID);req(digest(x)==meta["owner_sha"],"logical digest")
    alt=ds.copy(deep=True);alt.q.values[np.flatnonzero(alt.active.values)[0]]+=1e-4;reject("q",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.E_sensible_J.values[np.flatnonzero(alt.core.values)[0]]+=.001;reject("E",lambda:verify(alt))
    alt=ds.copy(deep=True);alt["V"]=("node",total(alt));reject("duplicate",lambda:verify(alt))
    alt=ds.copy(deep=True);i=np.flatnonzero(alt.core.values)[0];alt.W.values[i]=geometry(total(alt),A)[1][i]+.01;reject("capacity",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.attrs["step"]+=1;reject("clock",lambda:verify(alt))
    alt=ds.copy(deep=True);alt.V_class.values[np.flatnonzero(alt.core.values)[0],0]=-.01;reject("negative_class",lambda:verify(alt))
    RBD=load_provider(out);m=Model(ds,RBD);req(digest(m.ds)==meta["owner_sha"],"restored")
    rows=[];snaps={}
    for k in range(MID+1,STEPS+1):
        a,_=m.advance(k);rows.append(a)
        if k%CYCLE==0:snaps[str(k)]=snap(m.ds,out/f"step{k:04d}.nc");observe(m.ds)
    fin=snaps[str(STEPS)];series(rows,out/"series.npz")
    rep={"steps":len(rows),"first_step":rows[0]["step"],"input_owner_sha":meta["owner_sha"],"final":fin,"snapshots":snaps,"controls":controls,"outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_HORIZON_RESTART_JSON "+json.dumps(rep,sort_keys=True),flush=True)
def compare(reference,consumer_path,dest):
    p=pathlib.Path(reference);c=pathlib.Path(consumer_path);a=json.loads((p/"result.json").read_text());b=json.loads((c/"result.json").read_text())
    req(not a["failed_preregistered_gates"] and b["steps"]==MID and b["first_step"]==MID+1 and len(b["controls"])==7 and all(b["controls"].values()),"source/controls")
    req(a["cases"]["periodic"]["snapshots"][str(MID)]["owner_sha"]==b["input_owner_sha"] and a["cases"]["periodic"]["snapshots"][str(STEPS)]["owner_sha"]==b["final"]["owner_sha"],"owner")
    with np.load(p/"periodic/series.npz",allow_pickle=False) as x,np.load(c/"series.npz",allow_pickle=False) as y:
        req(set(x.files)==set(y.files) and all(np.array_equal(x[k][MID:],y[k]) for k in x.files),"suffix arrays")
    for k in range(MID+CYCLE,STEPS+1,CYCLE):req(digest(read(p/"periodic"/f"step{k:04d}.nc"))==digest(read(c/f"step{k:04d}.nc")),"snapshot suffix")
    rep={"suffix_exact":True,"steps":MID,"controls":b["controls"],"owner_sha":b["final"]["owner_sha"]}
    pathlib.Path(dest).write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_HORIZON_COMPARE_JSON "+json.dumps(rep,sort_keys=True),flush=True)
if __name__=="__main__":{"produce":producer,"consume":consumer,"compare":compare}[sys.argv[1]](*sys.argv[2:])
