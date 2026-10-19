"use strict";
const fs=require("node:fs"),path=require("node:path");
const K=require("./kernel.cjs"),protocol=require("./protocol.json"),profile=require("./profile.json");
function execute(spec,onProgress){
  const cfg={...protocol.base,...spec,objects:K.caseObjects(spec.name,protocol.base)};
  const r=K.makeRun(profile,cfg),rows=[],profiles=[];
  const diagEvery=Math.round(protocol.measurements.diagnostic_period/spec.dt);
  const count=Math.round(spec.T/spec.dt);
  if(Math.abs(count*spec.dt-spec.T)>1e-9||Math.abs(diagEvery*spec.dt-1)>1e-9)throw Error("Time grid mismatch");
  for(let step=0;step<=count;step++){
    if(step%diagEvery===0){
      const d=r.diagnostics();rows.push(d.summary);
      if(protocol.measurements.profile_times.some(t=>Math.abs(t-d.summary.time)<1e-9))
        profiles.push({time:d.summary.time,lineDensity:d.lineDensity,lineCharge:d.lineCharge,axisF2:d.axisF2,axisS:d.axisS});
      if(onProgress)onProgress(d.summary);
    }
    if(step<count)r.step();
  }
  return {id:spec.id,configuration:cfg,grid:r.grid,analysis:K.analyzeRun(rows,cfg),trace:rows,profiles};
}
function controls(){
  const cfg={...protocol.base,R:24,Z:32,h:.5,dt:.1};
  const vacuum=K.makeRun(profile,{...cfg,objects:[]});vacuum.step(20);
  const vd=vacuum.diagnostics().summary;
  const r=K.makeRun(profile,{...cfg,objects:K.caseObjects("pair_phase_halfpi",cfg)});
  const snapshot=r.snapshot();r.step(100);r.reverseMomenta();r.step(100);
  const rev=r.compareSnapshot(snapshot,true);
  return {vacuum:vd,time_reversal:rev,
    pass:vd.Ebar===0&&vd.Qbar===0&&vd.max_F===0&&vd.min_S===1
      &&rev.maximum_field_error<1e-9&&rev.maximum_momentum_error<1e-9};
}
if(require.main===module){
  const out=path.join(__dirname,"results");fs.mkdirSync(out,{recursive:true});
  const selected=process.argv[2];
  const specs=selected?protocol.runs.filter(x=>x.id===selected):protocol.runs;
  if(!specs.length)throw Error("Unknown run id");
  const c=controls();fs.writeFileSync(path.join(out,"controls.json"),JSON.stringify(c,null,2)+"\n");
  if(!c.pass)throw Error("Control failure");
  for(const spec of specs){
    const result=execute(spec);fs.writeFileSync(path.join(out,spec.id+".json"),JSON.stringify(result)+"\n");
    console.log(JSON.stringify({id:spec.id,...result.analysis}));
  }
}
module.exports={execute,controls};
