"""P4.1: sensible energy on retained hydromorphodynamic trajectory.
Single column temperature, published radiative/sensible terms, no evaporation.
"""
import hashlib,io,json,pathlib,resource,sys,time,traceback,zipfile
import numpy as np
import xarray as xr
from p3_exner import req,DX,DT,PHI,RHO_S
from p3_pore import remap,geometry,input_hash
from p3_water import digest,sha
ARCH=pathlib.Path("evidence/P3/P33/15028222831.zip")
ARCH_SHA="36814182ed6d793d1d0a40a0b2745601605e1171f08b45ace5a47df29bc2c16c"
A=DX*DX;STEPS=240;MID=120;TR=273.15
CS=800.*RHO_S;CW=4180.*1000.;SIGMA=5.670374419e-8
MODES=("forced","dark","transport_only","heterogeneous")
LEDGER=("rain_J","water_out_J","solid_out_J","shortwave_J","longwave_in_J","longwave_out_J","sensible_in_J","sensible_out_J")
def capacity(V,W,H):return CS*V+CW*(W+H)
def temp(E,C,core):
    req(np.isfinite(E).all() and np.isfinite(C).all() and np.all(C[core]>0) and np.all(E[core]>0) and np.all(C[~core]==0) and np.all(E[~core]==0),"energy/capacity")
    T=TR+np.divide(E,C,out=np.zeros_like(E),where=C>0)
    req(np.all((T[core]>TR)&(T[core]<373.15)),"temperature domain")
    return T
def transfer(stock,T,f,t,head,core,active,volcp):
    args=(stock,T,f,t,head,core,active);finger=input_hash(*args)
    stock,T,f,t,head,core,active=map(np.asarray,args)
    n=len(stock)
    req(stock.shape==T.shape==core.shape==(n,) and f.shape==t.shape==head.shape==active.shape,"transfer shape")
    req(all(np.isfinite(x).all() for x in (stock,T,f)) and np.all(stock>=0) and np.all(stock[~core]==0) and np.all(T[core]>TR),"transfer bounds")
    req(np.all(f[~active]==0),"inactive")
    source=np.where(f>=0,t,head);target=np.where(f>=0,head,t);mag=np.abs(f)
    req(np.all(mag[~core[source]]==0),"undeclared inflow")
    gross=np.zeros(n);np.add.at(gross,source,mag);req(np.all(gross<=stock+1e-12),"overdraw")
    ent=mag*volcp*(T[source]-TR)
    ev=np.zeros(n);dv=np.zeros(n)
    np.add.at(ev,source,-ent);np.add.at(ev,target,ent)
    np.add.at(dv,source,-mag);np.add.at(dv,target,mag)
    exp=float(ev[~core].sum());out=float(dv[~core].sum())
    ev[~core]=0.;dv[~core]=0.
    req(abs(float(ev.sum())+exp)<=.001 and abs(float(dv.sum())+out)<=1e-10,"transfer conservation")
    req(input_hash(*args)==finger,"transfer mutated input")
    return ev,dv,exp,out,np.sign(f)*ent,gross
def fixture(out):
    s=np.array([1.,1.,0.]);T=np.array([303.15,283.15,TR]);f=np.array([.1,.2]);t=np.array([0,1]);hd=np.array([1,2]);c=np.array([True,True,False]);act=np.array([True,True])
    ev,dv,exp,vout,ef,g=transfer(s,T,f,t,hd,c,act,CW)
    req(np.allclose(ev,[-12540000.,4180000.,0.],rtol=0,atol=.001) and abs(exp-8360000.)<=.001 and np.allclose(dv,[-.1,-.1,0.],rtol=0,atol=1e-12) and abs(vout-.2)<=1e-12,"analytic transfer")
    zero=transfer(s,T,np.zeros(2),t,hd,c,act,CW);req(all(np.all(np.asarray(x)==0) for x in zero),"zero transfer")
    controls={}
    def reject(name,fn):
        try:fn()
        except ValueError:controls[name]=True
        else:raise ValueError("accepted "+name)
    reject("overdraw",lambda:transfer(s,T,np.array([1.1,.2]),t,hd,c,act,CW))
    reject("external_inflow",lambda:transfer(s,T,np.array([.1,-.2]),t,hd,c,act,CW))
    reject("inactive",lambda:transfer(s,T,f,t,hd,c,np.array([False,True]),CW))
    reject("nan",lambda:transfer(s,T,np.array([float("nan"),.2]),t,hd,c,act,CW))
    reject("energy_capacity",lambda:temp(np.array([1.,1.,1.]),s,c))
    np.savez_compressed(out/"fixture.npz",energy_delta=ev,volume_delta=dv,energy_export=exp,volume_export=vout,energy_link=ef,gross=g)
    return {"energy_delta_J":ev.tolist(),"export_J":exp,"controls":controls}
def reference():
    req(sha(ARCH)==ARCH_SHA,"P3 archive digest")
    with zipfile.ZipFile(ARCH) as z:
        req(z.testzip() is None,"archive CRC");a=json.loads(z.read("reference/result.json"))
        req(not a["failed_preregistered_gates"] and a["cases"]["active"]["steps"]==STEPS,"P3 gates")
        ds={}
        for name in ("initial","mid","final"):
            p="active/"+name+".nc";b=z.read("reference/"+p)
            req(hashlib.sha256(b).hexdigest()==a["outputs"][p],"NC digest")
            with xr.open_dataset(io.BytesIO(b),engine="scipy") as f:ds[name]=f.load()
            req(digest(ds[name])==a["cases"]["active"][name]["owner_sha"],"NC owner digest")
        p="active/series.npz";b=z.read("reference/"+p);req(hashlib.sha256(b).hexdigest()==a["outputs"][p],"series digest")
        with np.load(io.BytesIO(b),allow_pickle=False) as f:x={k:f[k].copy() for k in f.files}
    req(np.array_equal(x["step"],np.arange(1,STEPS+1)),"steps")
    for name,k in (("mid",MID-1),("final",STEPS-1)):
        for field in ("V","W","H","q","hl"):req(np.array_equal(ds[name][field].values,x["trace_"+field][k]),"trace/snapshot "+name+field)
        req(x["owner_sha"][k]==digest(ds[name]),"trace snapshot logical hash")
    return ds,x,a
def snapshot(hydro,E,mode,ledger,initial_energy,p):
    out=hydro.copy(deep=True);out["E_sensible_J"]=("node",E.copy(),{"units":"J"})
    out.attrs.update(thermal_schema="P4.1-v1",thermal_mode=mode,temperature_reference_K=TR,solid_volumetric_heat_capacity_J_m3_K=CS,water_volumetric_heat_capacity_J_m3_K=CW,thermal_initial_J=initial_energy,**ledger)
    req(set(out.data_vars)==set(hydro.data_vars)|{"E_sensible_J"},"owner schema")
    core=out.core.values.astype(bool);d=digest(out);temp(out.E_sensible_J.values,capacity(out.V.values,out.W.values,out.H.values),core)
    req(digest(out)==d,"observation mutated")
    out.to_netcdf(p,engine="scipy")
    with xr.open_dataset(p,engine="scipy") as f:back=f.load()
    req(digest(back)==d,"snapshot roundtrip")
    return {"sha256":sha(p),"owner_sha":d,"hydro_owner_sha":digest(hydro)}
def total_balance(E,ledger,E0):
    return float(E.sum())+ledger["water_out_J"]+ledger["solid_out_J"]+ledger["longwave_out_J"]+ledger["sensible_out_J"]-ledger["rain_J"]-ledger["shortwave_J"]-ledger["longwave_in_J"]-ledger["sensible_in_J"]-E0
def run(mode,ds,x,out):
    ini=ds["initial"];c=ini.core.values.astype(bool);act=ini.active.values.astype(bool);t=ini["tail"].values;head=ini["head"].values
    V,W,H=(ini[k].values.copy() for k in ("V","W","H"))
    C=capacity(V,W,H);T0=np.where(c,297.+.5*ini.x.values if mode=="heterogeneous" else 300.,TR)
    E=C*(T0-TR);E0=float(E.sum());ledger={k:0. for k in LEDGER}
    folder=out/mode;folder.mkdir();snaps={"initial":snapshot(ini,E,mode,ledger,E0,folder/"initial.nc")}
    rows=[];arrays=[]
    raincum=outcum=relcum=infcum=0.
    for k in range(1,STEPS+1):
        j=k-1;old_E=float(E.sum());T=temp(E,capacity(V,W,H),c)
        rain_depth=(72. if k<=60 else 36. if 120<k<=180 else 0.)/3.6e6*DT
        rain=np.where(c,rain_depth*A,0.);rain_E=rain*CW*(300.-TR)
        e_w,d_w,ewout,wout,linkw,_=transfer(H,T,x["trace_q"][j]*DX*DT,t,head,c,act,CW)
        Hflow=np.where(c,x["trace_flow_h"][j]*A,0.)
        req(np.max(np.abs(H+rain+d_w-Hflow))<=1e-10,"native water material")
        Eflow=E+rain_E+e_w
        op1=abs(float(Eflow.sum())+ewout-float(rain_E.sum())-old_E)
        Tf=temp(Eflow,capacity(V,W,Hflow),c)
        fs=x["trace_qb"][j]*DX*DT
        e_s,d_s,esout,sout,links,_=transfer(V,Tf,fs,t,head,c,act,CS)
        r,a=remap(V,W,Hflow,fs,t,head,c,act,A)
        Vn,Wn,Hn=(x["trace_"+n][j].copy() for n in ("V","W","H"))
        req(np.max(np.abs(V+d_s-Vn))<=1e-12 and np.max(np.abs(a["V_after"]-Vn))<=1e-12,"native solid stock")
        zn=np.where(c,ini.bedrock.values+geometry(Vn,A)[0],ini.z0.values)
        req(np.max(np.abs(zn-x["trace_z"][j]))<=1e-12,"native geometry")
        infiltration=Wn-a["W_after"]
        req(np.all(infiltration>=0) and np.max(np.abs(Hn+infiltration-a["H_after"]))<=1e-12 and np.all(Wn<=geometry(Vn,A)[1]+1e-12),"native pore/GA")
        req(sout==esout==r["solid_out_m3"]==0,"solid boundary")
        Eadv=Eflow+e_s;op2=abs(float(Eadv.sum())+esout-float(Eflow.sum()))
        Cn=capacity(Vn,Wn,Hn);Ta=temp(Eadv,Cn,c)
        sw=np.zeros_like(E);li=sw.copy();lo=sw.copy();si=sw.copy();so=sw.copy()
        if mode in ("forced","dark"):
            if mode=="forced" and k<=MID:sw[c]=.8*600.*A*DT
            li[c]=.95*SIGMA*285.**4*A*DT;lo[c]=.95*SIGMA*Ta[c]**4*A*DT
            sens=10.*(295.-Ta[c])*A*DT
            si[c]=np.maximum(sens,0.);so[c]=np.maximum(-sens,0.)
        Enew=Eadv+sw+li-lo+si-so;op3=abs(float(Enew.sum())-float(Eadv.sum())-float((sw+li-lo+si-so).sum()))
        req(max(op1,op2,op3)<=.001,"operator energy")
        updates=(float(rain_E.sum()),ewout,esout,float(sw.sum()),float(li.sum()),float(lo.sum()),float(si.sum()),float(so.sum()))
        for name,val in zip(LEDGER,updates):ledger[name]+=val
        res=total_balance(Enew,ledger,E0);req(abs(res)<=.01,"global energy balance")
        raincum+=float(rain.sum());outcum+=wout;relcum+=r["released_water_m3"];infcum+=float(infiltration.sum())
        for val,key in ((raincum,"rain"),(outcum,"out"),(relcum,"release"),(infcum,"infiltration")):req(abs(val-x[key][j])<=1e-8,"P3 ledger "+key)
        req(x["inc"][j]==x["solid_out"][j]==0,"P3 boundary")
        Tnew=temp(Enew,Cn,c)
        rows.append({"step":k,"T_mean_K":float(TR+Enew.sum()/Cn.sum()),"T_min_K":float(Tnew[c].min()),"T_max_K":float(Tnew[c].max()),"T_uniform_error_K":float(np.max(np.abs(Tnew[c]-300.))),"energy_residual_J":res,"max_operator_residual_J":max(op1,op2,op3),**ledger})
        arrays.append({"E_J":Enew.copy(),"T_K":Tnew,"water_link_energy_J":linkw,"solid_link_energy_J":links,"rain_energy_J":rain_E,"shortwave_J_node":sw,"longwave_in_J_node":li,"longwave_out_J_node":lo,"sensible_in_J_node":si,"sensible_out_J_node":so})
        V,W,H,E=Vn,Wn,Hn,Enew
        if k in (MID,STEPS):
            name="mid" if k==MID else "final"
            for key,val in (("V",V),("W",W),("H",H)):req(np.array_equal(ds[name][key].values,val),"snapshot stock")
            snaps[name]=snapshot(ds[name],E,mode,ledger,E0,folder/(name+".nc"))
    np.savez_compressed(folder/"series.npz",**{k:np.asarray([r[k] for r in rows]) for k in rows[0]},**{k:np.asarray([a[k] for a in arrays]) for k in arrays[0]})
    last=rows[-1];mid=rows[MID-1];Tlast=arrays[-1]["T_K"]
    req(last["water_out_J"]>0 and last["solid_out_J"]==0 and last["rain_J"]>0,"export/rain energy")
    if mode=="transport_only":req(max(r["T_uniform_error_K"] for r in rows)<=1e-9,"uniform advection")
    if mode=="heterogeneous":
        req(float(np.max(np.abs(Tlast[c]-T0[c])))>1e-4,"heterogeneous transfer")
        req(all(r["T_min_K"]>=min(T0[c].min(),300.)-1e-9 and r["T_max_K"]<=max(T0[c].max(),300.)+1e-9 for r in rows),"advective extrema")
    if mode=="forced":
        req(mid["T_mean_K"]>300.001 and last["T_mean_K"]<mid["T_mean_K"]-.001,"heating/cooling")
        req(abs(last["shortwave_J"]-576000.)<=1e-6,"absorbed pulse")
    else:req(last["shortwave_J"]==0,"shortwave control")
    if mode=="dark":req(last["T_mean_K"]<299.999,"dark cooling")
    if mode in ("forced","dark"):req(last["longwave_out_J"]>0 and last["sensible_out_J"]>0,"dissipation")
    else:req(all(last[n]==0 for n in LEDGER[3:]),"transport-only forcing")
    return {"steps":STEPS,"snapshots":snaps,"mid":mid,"last":last,"max_energy_residual_J":max(abs(r["energy_residual_J"]) for r in rows),"max_operator_residual_J":max(r["max_operator_residual_J"] for r in rows),"max_uniform_error_K":max(r["T_uniform_error_K"] for r in rows),"max_final_temperature_change_K":float(np.max(np.abs(Tlast[c]-T0[c])))}
def main(dest):
    out=pathlib.Path(dest);out.mkdir();start=time.monotonic();rep={"experiment":"P4.1","cases":{},"failed_preregistered_gates":[]};fail=rep["failed_preregistered_gates"]
    try:
        rep["fixture"]=fixture(out);ds,x,a=reference();finger=input_hash(*[x[k] for k in sorted(x)])
        before={k:digest(v) for k,v in ds.items()}
        for mode in MODES:
            try:rep["cases"][mode]=run(mode,ds,x,out)
            except Exception as e:traceback.print_exc();fail.append(mode+":"+repr(e))
        req(finger==input_hash(*[x[k] for k in sorted(x)]) and before=={k:digest(v) for k,v in ds.items()},"retained data mutation")
        rep["source_unchanged"]=True;rep["source_archive_sha256"]=ARCH_SHA
    except Exception as e:traceback.print_exc();fail.append("input/fixture:"+repr(e))
    if len(rep["cases"])!=4:fail.append("incomplete cases")
    rep["wall_seconds"]=time.monotonic()-start;rep["peak_rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    rep["outputs"]={str(p.relative_to(out)):sha(p) for p in out.rglob("*") if p.is_file()}
    (out/"result.json").write_text(json.dumps(rep,indent=2,sort_keys=True)+"\n")
    print("WORLD_ENERGY_JSON "+json.dumps(rep,sort_keys=True),flush=True)
    if fail:raise SystemExit(1)
if __name__=="__main__":main(sys.argv[1])
