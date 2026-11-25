"""P5.2: live mixed-size sediment with exclusive class ownership; P4.2 operators retained."""
import hashlib,json,pathlib,sys,time,resource,traceback
import numpy as np
import xarray as xr
from p3_reciprocal import Model as Hydro,initial as hydro_initial,zphys as hydro_zphys,rain,load_provider,DX,DT,A,PHI,RHO_S,FILM,MID,STEPS,LEDGER as HYDRO_LEDGER
from p3_exner import grid,req
from p3_pore import remap,ga,geometry,domain,input_hash
from p3_water import digest,sha
from p4_energy import capacity,temp,transfer,total_balance,TR,CS,CW,SIGMA,LEDGER as ENERGY_LEDGER
from p5_grain import move
MODES=("active","passive","frozen")
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
    ds.attrs.update(thermal_schema="P5.2-v1",thermal_mode=mode,thermal_initial_J=float(E.sum()),temperature_reference_K=TR,latent_J_kg=L,CE=CE,air_T_K=AIR,air_RH=RH,**{n:0. for n in ENERGY_LEDGER+VAPOR})
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
    rr=sum(rain("active",i)*DT*A*int(c.sum()) for i in range(1,k+1))
    req(abs(a["rain"]-rr)<=1e-10 and a["solid_out"]==a["solid_out_J"]==a["inc"]==0,"boundary/forcing")
    req(abs(a["rain_J"]-rr*CW*(300.-TR))<=.001,"rain energy")
    req(abs(a["shortwave_J"]-(min(k,MID)*.8*600*A*DT*int(c.sum())))<=1e-6,"pulse")
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
        owner=self.ds;self.ds=view(owner)
        try:return super().flow(k,freeze)
        finally:self.ds=owner
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
        z=zphys(ds);fine=fraction(ds);used=frac_initial(ds) if mode=="frozen" else fine
        qb,p,nativez=proposal(self.RBD,ds,h,hl,u,prev,used)
        fs=qb*DX*DT;frac=np.column_stack((1.-fine,fine))
        chosen=frac[np.where(fs>=0,self.t,self.hd)] if mode=="passive" else p
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
        if k<=MID:sw[c]=.8*600.*A*DT
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
        metrics={"step":k,"class_res":float(np.max(np.abs(Mn.sum(axis=0)-initial(mode).V_class.values.sum(axis=0)))),"owner_sha":digest(nxt),**{n:at[n] for n in LEDGERS},"solid_res":sr,"water_res":wr,"energy_res":er,"flow_local":flowlocal,"exner_local":exlocal,"pore_local":rem["max_local_water_residual_m3"],"ga_local":gr["max_local_water_residual_m3"],"evap_local_K":ev["temperature_error"],"max_energy_op":max(op1,op2,op3,ev["energy_error"]),"T_mean_K":float(TR+En.sum()/Cn.sum()),"min_F":rem["min_F_equiv_m"],"max_dz":float(np.max(np.abs(newz-ds.z0.values)))}
        trace={"V":Vn.copy(),"V_class":Mn.copy(),"p":p.copy(),"chosen_p":chosen.copy(),"W":Wn.copy(),"H":Hn.copy(),"q":q,"hl":hl,"flow_h":h,"qb":qb,"z":newz,"E":En.copy(),"T":Tn,"H_before_evap":Hp.copy(),"E_before_evap":Ep.copy(),"vapor_mass":ev["m"],"vapor_demand":ev["demand"],"vapor_latent":ev["latent"],"vapor_sensible":ev["sensible"],"water_link_energy":linkw,"solid_link_energy":links}
        return metrics,trace
def series(rows,traces,p):
    np.savez_compressed(p,**{n:np.asarray([r[n] for r in rows]) for n in rows[0]},**{"trace_"+n:np.asarray([t[n] for t in traces]) for n in traces[0]})
def producer(dest,checkpoint):
    out=pathlib.Path(dest);out.mkdir();cp=pathlib.Path(checkpoint);cp.mkdir();start=time.monotonic()
    rep={"experiment":"P5.2","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"];finals={}
    RBD=load_provider(out)
    for mode in MODES:
        folder=out/mode;folder.mkdir();rows=[];traces=[]
        try:
            m=Model(initial(mode),RBD);s0=snap(m.ds,folder/"initial.nc")
            for k in range(1,STEPS+1):
                if mode=="active" and k==MID+1:
                    d=digest(m.ds);branch=Model(m.ds,RBD);h,q,hl,u,prev,*_=branch.flow(k,False)
                    qb0,p0,z0=proposal(RBD,m.ds,h,hl,u,prev,fraction(m.ds))
                    qbf,pf,zf=proposal(RBD,m.ds,h,hl,u,prev,frac_initial(m.ds))
                    fork={"input_owner_sha":d,"qb_diff":float(np.max(np.abs(qb0-qbf))),"p_diff":float(np.max(np.abs(p0-pf)))}
                    req(digest(m.ds)==d,"fork mutated owner")
                    np.savez_compressed(folder/"fork.npz",real_qb=qb0,frozen_qb=qbf,real_p=p0,frozen_p=pf)
                row,tr=m.advance(k);rows.append(row);traces.append(tr)
                if mode=="active" and k==MID+1:
                    fork["real_exact"]=bool(np.array_equal(qb0,tr["qb"]) and np.array_equal(p0,tr["p"]));rep["fork"]=fork
                if k==MID:
                    sm=snap(m.ds,folder/"mid.nc");observe(m.ds)
                    if mode=="active":
                        meta=snap(m.ds,cp/"owner.nc");req(meta["owner_sha"]==sm["owner_sha"],"checkpoint")
                        (cp/"manifest.json").write_text(json.dumps({"step":MID,"mode":mode,**meta},sort_keys=True)+"\n")
            sf=snap(m.ds,folder/"final.nc");observe(m.ds);series(rows,traces,folder/"series.npz")
            finals[mode]=m.ds.copy(deep=True)
            rep["cases"][mode]={"initial":s0,"mid":sm,"final":sf,"steps":len(rows),"last":rows[-1],"max_residuals":{n:max(abs(r[n]) for r in rows) for n in ("class_res","solid_res","water_res","energy_res","flow_local","exner_local","pore_local","ga_local","evap_local_K","max_energy_op")},"min_F":min(r["min_F"] for r in rows)}
        except Exception as e:
            traceback.print_exc();fail.append(mode+":"+repr(e))
            if rows:series(rows,traces,folder/"partial.npz")
    def gate(n,v):
        if not bool(v):fail.append(n)
    if len(finals)==3:
        a=finals["active"];b=finals["passive"];f=finals["frozen"]
        af=folder.parent/"active"/"series.npz";ff=folder.parent/"frozen"/"series.npz"
        with np.load(af,allow_pickle=False) as x,np.load(ff,allow_pickle=False) as y:qb_diff=float(np.max(np.abs(x["trace_qb"][-1]-y["trace_qb"][-1])))
        rep["comparisons"]={"fine_change":float(np.max(np.abs(fraction(a)-frac_initial(a)))),"active_passive_stock_diff_m3":float(np.max(np.abs(a.V_class.values-b.V_class.values))),"active_frozen_qb_diff":qb_diff,"active_frozen_q_diff":float(np.max(np.abs(a.q.values-f.q.values)))}
        v=rep["comparisons"];gate("composition_change",v["fine_change"]>1e-5);gate("sorting",v["active_passive_stock_diff_m3"]>1e-6)
        gate("feedback",v["active_frozen_qb_diff"]>1e-8 and v["active_frozen_q_diff"]>1e-9)
        f=rep["fork"];gate("fork",f["qb_diff"]>1e-8 and f["p_diff"]>1e-5 and f["real_exact"])
    else:fail.append("incomplete cases")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_COMPOSITION_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
def consumer(checkpoint,dest):
    cp=pathlib.Path(checkpoint);out=pathlib.Path(dest);out.mkdir();meta=json.loads((cp/"manifest.json").read_text())
    req(meta["step"]==MID and meta["mode"]=="active","checkpoint metadata");ds=read(cp/"owner.nc",meta,MID)
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
    rows=[];traces=[]
    for k in range(MID+1,STEPS+1):
        a,b=m.advance(k);rows.append(a);traces.append(b)
    fin=snap(m.ds,out/"final.nc");series(rows,traces,out/"series.npz")
    rep={"steps":len(rows),"first_step":rows[0]["step"],"input_owner_sha":meta["owner_sha"],"final":fin,"controls":controls,"outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_COMPOSITION_RESTART_JSON "+json.dumps(rep,sort_keys=True),flush=True)
def compare(reference,consumer_path,dest):
    p=pathlib.Path(reference);c=pathlib.Path(consumer_path);a=json.loads((p/"result.json").read_text());b=json.loads((c/"result.json").read_text())
    req(not a["failed_preregistered_gates"] and b["steps"]==MID and b["first_step"]==MID+1 and len(b["controls"])==7 and all(b["controls"].values()),"source/controls")
    req(a["cases"]["active"]["mid"]["owner_sha"]==b["input_owner_sha"] and a["cases"]["active"]["final"]["owner_sha"]==b["final"]["owner_sha"],"owner")
    with np.load(p/"active/series.npz",allow_pickle=False) as x,np.load(c/"series.npz",allow_pickle=False) as y:
        req(set(x.files)==set(y.files) and all(np.array_equal(x[k][MID:],y[k]) for k in x.files),"suffix arrays")
    req(digest(read(p/"active/final.nc"))==digest(read(c/"final.nc")),"final state")
    rep={"suffix_exact":True,"steps":MID,"controls":b["controls"],"owner_sha":b["final"]["owner_sha"]}
    pathlib.Path(dest).write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n");print("WORLD_COMPOSITION_COMPARE_JSON "+json.dumps(rep,sort_keys=True),flush=True)
if __name__=="__main__":{"produce":producer,"consume":consumer,"compare":compare}[sys.argv[1]](*sys.argv[2:])
