"use strict";
const fs=require("node:fs"),path=require("node:path"),{execute}=require("./run.cjs"),p=require("./protocol.json");
function controls(){
  const initial=[];
  for(const h of [.5,1/3,.25])for(const sigma of [4,6,8]){
    const x=execute({id:"initial",sigma,h,dt:.1*h,R:48,Z:64,T:0,control:"full"},p);
    const s=x.trace[0],I=p.Qbar/p.preparation_frequency,w=p.preparation_frequency;
    const analytic=.5*(w*w+1)*I+3*I/(4*sigma*sigma)+p.physical_b*I*I/(4*Math.pow(2*Math.PI,1.5)*Math.pow(sigma,3));
    initial.push({h,sigma,initial_S_min:s.min_S,Qbar:s.Qbar,relative_charge_error:s.Qbar/p.Qbar-1,Ebar:s.Ebar,analytic_Ebar:analytic,relative_energy_error:s.Ebar/analytic-1});
  }
  const full=execute({id:"response",sigma:6,h:.5,dt:.05,R:24,Z:32,T:2,control:"full"},p);
  const rigid=execute({id:"rigid_response",sigma:6,h:.5,dt:.05,R:24,Z:32,T:2,control:"rigid_phi"},p);
  return {initial_comparisons:initial,full_final_min_S:full.trace.at(-1).min_S,rigid_final_min_S:rigid.trace.at(-1).min_S,
    pass:initial.every(x=>x.initial_S_min===1)&&initial.filter(x=>x.h===.25).every(x=>Math.abs(x.relative_charge_error)<.001&&Math.abs(x.relative_energy_error)<.001)&&full.trace.at(-1).min_S<.99&&rigid.trace.at(-1).min_S===1};
}
if(require.main===module){const c=controls();fs.writeFileSync(path.join(__dirname,"controls.json"),JSON.stringify(c,null,2)+"\n");console.log(JSON.stringify(c,null,2));if(!c.pass)process.exitCode=1;}
module.exports={controls};
