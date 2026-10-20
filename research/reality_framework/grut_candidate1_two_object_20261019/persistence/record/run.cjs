"use strict";
const fs=require("node:fs"),path=require("node:path"),K=require("../kernel_modes.cjs"),D=require("../diagnostics.cjs"),base=require("../../protocol.json").base,profile=require("../../profile.json");
function execute(spec,protocol,onProgress){
  const objects=K.caseObjects("pair_phase_halfpi",base);objects[1].phase=spec.sign*Math.PI/2;
  objects[0].velocity=spec.direction*base.speed;objects[1].velocity=-spec.direction*base.speed;
  const cfg={...base,...spec,objects,modes:[]},r=K.makeRun(profile,cfg),trace=[],profiles=[];
  const steps=Math.round(spec.T/spec.dt),period=Math.round(protocol.diagnostic_period/spec.dt);let first=null,maxE=0,maxQ=0;
  const begin=Date.now();
  for(let k=0;k<=steps;k++){
    if(k%period===0||k===steps){
      const f=r.diagnostics(),s=f.summary,local=D.ballDiagnostics(r,f,cfg);if(!first)first=s;
      maxE=Math.max(maxE,Math.abs(s.Ebar-first.Ebar)/first.Ebar);maxQ=Math.max(maxQ,Math.abs(s.Qbar-first.Qbar)/first.absQbar);
      const balls=local.balls.filter(x=>x.radius===protocol.receiver_radius&&!x.overlaps_other_ball&&x.center_z>0).sort((a,b)=>b.center_z-a.center_z);
      const receiver=balls[0]||null,q=receiver?.Qbar??null,threshold=protocol.reference_charge,margin=protocol.undecided_fraction*threshold;
      const decode=q===null?null:(q>threshold+margin?1:(q<threshold-margin?-1:0));
      trace.push({...s,...local,receiver_Qbar:q,decoded:decode});
      if(protocol.profile_times.some(t=>Math.abs(t-s.time)<1e-7))profiles.push({time:s.time,lineDensity:f.lineDensity,lineCharge:f.lineCharge,axisF2:f.axisF2,axisS:f.axisS});
      if(!Object.values(s).every(Number.isFinite))throw Error("Nonfinite background");
      if(onProgress&&Math.abs(s.time/20-Math.round(s.time/20))<1e-8)onProgress(trace.at(-1));
    }
    if(k<steps)r.step();
  }
  const window=trace.filter(x=>x.time>=spec.readout_start-1e-7),target=spec.direction===1?spec.sign:0;
  return {id:spec.id,configuration:cfg,wall_seconds:(Date.now()-begin)/1000,trace,profiles,maximum_relative_energy_drift:maxE,maximum_charge_error:maxQ,
    readout:{start:spec.readout_start,end:spec.T,target,samples:window.length,all_samples_match:window.every(x=>x.decoded===target),receiver_Q_range:[Math.min(...window.map(x=>x.receiver_Qbar??Infinity)),Math.max(...window.map(x=>x.receiver_Qbar??-Infinity))]},
    note:"An exploratory classical phase-to-charge record test. Ideal charge observable and decoder are supplied; no microscopic detector or conscious observer is implemented."};
}
if(require.main===module){
  const protocol=require("./protocol.json"),id=process.argv[2],spec=protocol.runs.find(x=>x.id===id);if(!spec)throw Error("Pass a prescribed id");
  fs.mkdirSync(path.join(__dirname,"results"),{recursive:true});
  const result=execute(spec,protocol,row=>console.log(JSON.stringify({id,time:row.time,receiver:row.receiver_Qbar,decode:row.decoded})));
  fs.writeFileSync(path.join(__dirname,"results",id+".json"),JSON.stringify(result)+"\n");console.log(JSON.stringify({completed:id,readout:result.readout,seconds:result.wall_seconds}));
}
module.exports={execute};
