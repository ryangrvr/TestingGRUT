"use strict";
const fs=require("node:fs"),path=require("node:path"),K=require("./kernel_modes.cjs"),D=require("./diagnostics.cjs"),base=require("../protocol.json").base,profile=require("../profile.json");
function execute(spec,protocol,onProgress){
  const cfg={...base,...spec,modes:protocol.modes,objects:K.caseObjects(spec.name,base)},r=K.makeRun(profile,cfg),trace=[],profiles=[];
  const steps=Math.round(spec.T/spec.dt),period=Math.round(protocol.diagnostic_period/spec.dt),frames=protocol.profile_times;
  let first=null,maxEnergy=0,maxCharge=0,firstLinearWarning=null;
  const started=Date.now();
  for(let k=0;k<=steps;k++){
    if(k%period===0||k===steps){
      const frame=r.diagnostics(),modes=r.modeDiagnostics(),local=D.ballDiagnostics(r,frame,cfg),s=frame.summary;
      if(!first)first=s;
      maxEnergy=Math.max(maxEnergy,Math.abs(s.Ebar-first.Ebar)/first.Ebar);
      maxCharge=Math.max(maxCharge,Math.abs(s.Qbar-first.Qbar)/first.absQbar);
      if(firstLinearWarning===null&&modes.some(x=>x.relative_to_background>=protocol.linear_warning))firstLinearWarning=s.time;
      trace.push({...s,modes,...local});
      if(frames.some(t=>Math.abs(t-s.time)<1e-7))profiles.push({time:s.time,lineDensity:frame.lineDensity,lineCharge:frame.lineCharge,axisF2:frame.axisF2,axisS:frame.axisS});
      if(!Object.values(s).every(Number.isFinite)||modes.some(x=>!Number.isFinite(x.norm_squared)))throw Error("Nonfinite state at "+s.time);
      if(onProgress&&Math.abs(s.time/20-Math.round(s.time/20))<1e-8)onProgress(trace.at(-1));
    }
    if(k<steps)r.step();
  }
  return {id:spec.id,configuration:cfg,grid:r.grid,wall_seconds:(Date.now()-started)/1000,trace,profiles,
    maximum_relative_energy_drift:maxEnergy,maximum_charge_error:maxCharge,first_linear_warning_time:firstLinearWarning,
    note:"Perturbations solve a linear tangent problem. Their amplitudes do not deform the evolved background; large tangent growth is not a completed nonlinear instability demonstration."};
}
if(require.main===module){
  const protocol=require("./protocol.json"),id=process.argv[2],spec=protocol.runs.find(x=>x.id===id);if(!spec)throw Error("Pass a prescribed run id");
  const out=path.join(__dirname,"results");fs.mkdirSync(out,{recursive:true});
  const result=execute(spec,protocol,row=>console.log(JSON.stringify({id,time:row.time,energy:row.Ebar,charge:row.Qbar,regions:row.centers,modes:row.modes})));
  fs.writeFileSync(path.join(out,id+".json"),JSON.stringify(result)+"\n");
  console.log(JSON.stringify({completed:id,seconds:result.wall_seconds,energy_drift:result.maximum_relative_energy_drift,charge_error:result.maximum_charge_error,linear_warning:result.first_linear_warning_time}));
}
module.exports={execute};
