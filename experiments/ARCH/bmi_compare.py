import json,pathlib,resource,sys,time
import numpy as np
from landlab.bmi import wrap_as_bmi
from landlab.components import SoilInfiltrationGreenAmpt
from p3_water import Model,initial,snapshot,read,DT,DX,SHAPE,N,STEPS,MID,req,digest
H="surface_water__depth"
F="soil_water_infiltration__depth"
class BMIInfiltration:
    def __init__(self,m,folder):
        self.h=m.h;self.f=m.F
        cfg={"grid":{"RasterModelGrid":[list(SHAPE),{"xy_spacing":DX},{"fields":{"node":{H:{"constant":[{"value":0.0}]},F:{"constant":[{"value":0.01}]}}}}]},
             "clock":{"start":m.owner.attrs["time_s"],"stop":STEPS*DT,"step":DT,"units":"s"},
             "soil_infiltration_green_ampt":{"hydraulic_conductivity":0.0 if m.case=="noinf" else 1e-6,"soil_bulk_density":1590.0,"rock_density":2650.0,"initial_soil_moisture_content":0.15,"soil_type":"sandy loam","surface_water_minimum_depth":1e-8}}
        p=folder/"bmi-config.json";p.write_text(json.dumps(cfg))
        self.b=wrap_as_bmi(SoilInfiltrationGreenAmpt)();self.b.initialize(str(p))
        self.t=m.owner.attrs["time_s"]
        req(self.b.get_current_time()==self.t and self.b.get_time_step()==DT and self.b.get_time_units()=="s","BMI clock")
        for name in (H,F):
            req(name in self.b.get_input_var_names() and name in self.b.get_output_var_names(),"BMI inout")
            req(self.b.get_var_units(name)=="m" and self.b.get_var_location(name)=="node","BMI metadata")
            grid=self.b.get_var_grid(name)
            req(self.b.get_grid_size(grid)==N,"BMI grid size")
            req(np.array_equal(self.b.get_grid_shape(grid,np.zeros(2,dtype=int)),SHAPE),"BMI grid shape")
            req(np.array_equal(self.b.get_grid_spacing(grid,np.zeros(2)),[DX,DX]),"BMI spacing")
        self.b.set_value("boundary_condition_flag",m.g.status_at_node.copy())
        self.b.set_value(H,self.h);self.b.set_value(F,self.f)
        temp=self.b.get_value(H,np.empty_like(self.h));temp[:]=-123.
        req(np.array_equal(self.b.get_value(H,np.empty_like(self.h)),self.h),"copy observation mutated")
    def run_one_step(self,dt):
        req(dt==DT and self.b.get_current_time()==self.t,"BMI clock before")
        self.b.set_value(H,self.h);self.b.set_value(F,self.f)
        self.b.update();self.t+=dt
        req(self.b.get_current_time()==self.t,"BMI clock after")
        self.b.get_value(H,self.h);self.b.get_value(F,self.f)
def run(dest,checkpoint=None):
    p=pathlib.Path(dest);p.mkdir(parents=True,exist_ok=True)
    if checkpoint:
        m=Model(read(pathlib.Path(checkpoint),step=MID));m.inf=BMIInfiltration(m,p)
        start=time.monotonic();rows=[m.advance(i) for i in range(MID+1,STEPS+1)]
        (p/"suffix.json").write_text(json.dumps(rows));snapshot(m.owner,p/"final.nc")
        print("ARCH04_CONSUMER",json.dumps({"seconds":time.monotonic()-start,"final":digest(m.owner),"steps":len(rows)}));return
    report={}
    for case in ("active","noinf"):
        modes={}
        for mode in ("direct","bmi"):
            folder=p/(case+"-"+mode);folder.mkdir()
            m=Model(initial(case))
            if mode=="bmi":m.inf=BMIInfiltration(m,folder)
            t=time.monotonic();rows=[]
            for i in range(1,STEPS+1):
                rows.append(m.advance(i))
                if case=="active" and mode=="bmi" and i==MID:
                    snapshot(m.owner,p/"checkpoint.nc")
            elapsed=time.monotonic()-t
            modes[mode]={"seconds":elapsed,"final":digest(m.owner),"rows":rows}
            snapshot(m.owner,folder/"final.nc")
            (folder/"rows.json").write_text(json.dumps(rows))
        req(modes["direct"]["rows"]==modes["bmi"]["rows"],"direct/BMI rows mismatch "+case)
        report[case]={k:{a:b for a,b in v.items() if a!="rows"} for k,v in modes.items()}
    report["copy_bytes_min_per_step"]=4*N*8
    report["rss_kib"]=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (p/"result.json").write_text(json.dumps(report,indent=2));print("ARCH04_PRODUCER",json.dumps(report))
if __name__=="__main__":run(*sys.argv[1:])
