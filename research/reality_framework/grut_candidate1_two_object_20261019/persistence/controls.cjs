"use strict";
const fs=require("node:fs"),path=require("node:path");
const old=require("../kernel.cjs"),K=require("./kernel_modes.cjs"),profile=require("../profile.json"),base=require("../protocol.json").base;
function maxDifference(a,b){let d=0;for(let k=0;k<a.length;k++)d=Math.max(d,Math.abs(a[k]-b[k]));return d;}
function runControls(){
  const cfg={...base,R:12,Z:16,h:.5,dt:.05,objects:old.caseObjects("pair_phase_halfpi",base)};
  const ref=old.makeRun(profile,cfg),ext=K.makeRun(profile,{...cfg,modes:[{m:1,seed:"relative_shift",epsilon:1e-4},{m:2,seed:"elliptic",epsilon:1e-4}]});
  ref.step(40);ext.step(40);const a=ref.snapshot(),b=ext.snapshot();
  let backgroundError=0;for(const key of ["U","W","S","P","Q","T"])backgroundError=Math.max(backgroundError,maxDifference(a[key],b[key]));
  const phase=K.makeRun(profile,{...cfg,modes:[{m:0,seed:"phase",epsilon:1e-4},{m:2,seed:"zero",epsilon:0}]});
  phase.step(80);let phaseError=0,zeroError=0;const s=phase.state,d=phase.modes[0].d,c=phase.modes[0].seedScale;
  for(let k=0;k<s[0].length;k++)phaseError=Math.max(phaseError,Math.abs(d[0][k]+c*s[1][k]),Math.abs(d[1][k]-c*s[0][k]),Math.abs(d[3][k]+c*s[4][k]),Math.abs(d[4][k]-c*s[3][k]),Math.abs(d[2][k]),Math.abs(d[5][k]));
  for(const arr of phase.modes[1].d)for(const x of arr)zeroError=Math.max(zeroError,Math.abs(x));
  const tangent=K.makeRun(profile,{...cfg,objects:old.caseObjects("single_static",base),modes:[{m:0,seed:"amplitude",epsilon:1}]});
  const initial=tangent.state.map(x=>Float64Array.from(x)),direction=tangent.modes[0].d.slice(0,6).map(x=>Float64Array.from(x));
  tangent.step(40);const finiteDifference=[];
  for(const eps of [1e-3,5e-4]){
    const plus=K.makeRun(profile,{...cfg,objects:old.caseObjects("single_static",base)}),minus=K.makeRun(profile,{...cfg,objects:old.caseObjects("single_static",base)});
    for(let a=0;a<6;a++)for(let k=0;k<initial[a].length;k++){plus.state[a][k]=initial[a][k]+eps*direction[a][k];minus.state[a][k]=initial[a][k]-eps*direction[a][k];}
    plus.recompute();minus.recompute();plus.step(40);minus.step(40);let error=0,norm=0;
    for(let a=0;a<6;a++)for(let k=0;k<initial[a].length;k++){
      const predicted=tangent.modes[0].d[a][k],difference=(plus.state[a][k]-minus.state[a][k])/(2*eps)-predicted;
      error+=difference*difference;norm+=predicted*predicted;
    }
    finiteDifference.push({epsilon:eps,relative_L2_error:Math.sqrt(error/norm)});
  }
  const angular=[];
  for(const m of [1,2]){
    const r=K.makeRun(profile,{...cfg,R:6,Z:6,objects:[],modes:[{m,seed:"zero",epsilon:0}]}),d=r.modes[0].d,{Nr,Nz,h}=r.grid;
    for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){const rho=(i+.5)*h,z=-6+(j+.5)*h;d[0][i*Nz+j]=Math.pow(rho,m)*z*z;}
    r.recompute();let error=0;
    for(let i=0;i<Nr-2;i++)for(let j=2;j<Nz-2;j++){
      const rho=(i+.5)*h,z=-6+(j+.5)*h;if(Math.hypot(rho,z)>4)continue;
      const target=Math.pow(rho,m)*(2-z*z);error=Math.max(error,Math.abs(d[6][i*Nz+j]-target));
    }
    angular.push({m,maximum_error:error,target:"L_m(rho^m z^2)-rho^m z^2 = rho^m(2-z^2)"});
  }
  return {background_bitwise_error:backgroundError,global_phase_tangent_error:phaseError,zero_mode_error:zeroError,
    tangent_finite_difference:finiteDifference,error_reduction:finiteDifference[0].relative_L2_error/finiteDifference[1].relative_L2_error,
    angular_manufactured:angular,
    pass:backgroundError===0&&phaseError<1e-10&&zeroError===0&&finiteDifference[1].relative_L2_error<.4*finiteDifference[0].relative_L2_error&&angular.every(x=>x.maximum_error<1e-10)};
}
if(require.main===module){const result=runControls();fs.writeFileSync(path.join(__dirname,"controls_preflight.json"),JSON.stringify(result,null,2)+"\n");console.log(JSON.stringify(result,null,2));if(!result.pass)process.exitCode=1;}
module.exports={runControls};
