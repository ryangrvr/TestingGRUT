"use strict";
const fs=require("node:fs"),path=require("node:path"),full=require("../kernel_modes.cjs"),rigid=require("./rigid_kernel.cjs"),D=require("../diagnostics.cjs"),base=require("../../protocol.json").base,profile=require("../../profile.json");
function execute(spec,protocol,onProgress){
  const cfg={...base,...spec,modes:[],objects:[]},K=spec.control==='rigid_phi'?rigid:full,r=K.makeRun(profile,cfg),[U,W,S,P,Q,T]=r.state;
  const A=Math.sqrt(protocol.Qbar/(protocol.preparation_frequency*Math.pow(Math.PI,1.5)*Math.pow(spec.sigma,3)));
  for(let i=0;i<r.grid.Nr;i++)for(let j=0;j<r.grid.Nz;j++){
    const rho=(i+.5)*spec.h,z=-spec.Z+(j+.5)*spec.h,k=i*r.grid.Nz+j;
    U[k]=A*Math.exp(-(rho*rho+z*z)/(2*spec.sigma*spec.sigma));W[k]=0;S[k]=1;
    P[k]=0;Q[k]=-protocol.preparation_frequency*U[k];T[k]=0;
  }
  r.recompute();const steps=Math.round(spec.T/spec.dt),period=Math.round(protocol.diagnostic_period/spec.dt),trace=[],profiles=[];
  let first=null,maxEnergy=0,maxCharge=0;const start=Date.now();
  for(let k=0;k<=steps;k++){
    if(k%period===0||k===steps){
      const f=r.diagnostics(),s=f.summary,local=D.ballDiagnostics(r,f,cfg);if(!first)first=s;
      maxEnergy=Math.max(maxEnergy,Math.abs(s.Ebar-first.Ebar)/first.Ebar);maxCharge=Math.max(maxCharge,Math.abs(s.Qbar-first.Qbar)/first.absQbar);
      if(!Object.values(s).every(Number.isFinite))throw Error("Nonfinite background at "+s.time);
      trace.push({...s,...local});
      if(protocol.profile_times.some(t=>Math.abs(t-s.time)<1e-7))profiles.push({time:s.time,lineDensity:f.lineDensity,lineCharge:f.lineCharge,axisF2:f.axisF2,axisS:f.axisS});
      if(onProgress&&Math.abs(s.time/20-Math.round(s.time/20))<1e-8)onProgress(trace.at(-1));
    }
    if(k<steps)r.step();
  }
  return {id:spec.id,configuration:cfg,initial_amplitude:A,prepared_Qbar:protocol.Qbar,wall_seconds:(Date.now()-start)/1000,
    trace,profiles,maximum_relative_energy_drift:maxEnergy,maximum_charge_error:maxCharge,
    note:spec.control==='rigid_phi'?"Counterfactual feedback ablation: S is held identically 1. The scalar X energy remains conserved; this is a comparison law, not the full coupled model.":"Full coupled scalar equations; only the initial Gaussian charged cloud is prepared. No preformed depleted-phi object, external damping, or later drive."};
}
if(require.main===module){
  const protocol=require("./protocol.json"),id=process.argv[2],spec=protocol.runs.find(x=>x.id===id);if(!spec)throw Error("Pass a protocol run id");
  fs.mkdirSync(path.join(__dirname,"results"),{recursive:true});
  const result=execute(spec,protocol,row=>console.log(JSON.stringify({id,time:row.time,Ebar:row.Ebar,Qbar:row.Qbar,maxF:row.max_F,minS:row.min_S,centers:row.centers})));
  fs.writeFileSync(path.join(__dirname,"results",id+".json"),JSON.stringify(result)+"\n");
  console.log(JSON.stringify({completed:id,seconds:result.wall_seconds,energy_drift:result.maximum_relative_energy_drift,charge_error:result.maximum_charge_error}));
}
module.exports={execute};
