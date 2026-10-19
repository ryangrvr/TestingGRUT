"use strict";
const K=require("./kernel.cjs"),protocol=require("./protocol.json"),profile=require("./profile.json");
function calculate(){
  const base=protocol.base,rows=[];
  for(const h of [.5,1/3,.25]){
    const cfg={...base,R:24,Z:32,h,dt:1e-6,objects:K.caseObjects("single_static",base)};
    const r=K.makeRun(profile,cfg),before=r.snapshot();r.step();const after=r.snapshot();
    let aNorm=0,sNorm=0,scale=0,mass=0;
    for(let i=0;i<r.grid.Nr;i++)for(let j=0;j<r.grid.Nz;j++){
      const rho=(i+.5)*h,z=-cfg.Z+(j+.5)*h,rr=Math.hypot(rho,z);if(rr>12)continue;
      const idx=i*r.grid.Nz+j,weight=2*Math.PI*rho*h*h;
      const uacc=(after.P[idx]-before.P[idx])/cfg.dt,sacc=(after.T[idx]-before.T[idx])/cfg.dt;
      const target=-base.Omega*base.Omega*before.U[idx];
      aNorm+=weight*(uacc-target)**2;sNorm+=weight*sacc*sacc;scale+=weight*target*target;mass+=weight;
    }
    rows.push({h,core_radius:12,relative_L2_X_stationary_acceleration_defect:Math.sqrt(aNorm/scale),phi_acceleration_RMS:Math.sqrt(sNorm/mass)});
  }
  const manufacturedProfile=Array.from({length:501},(_,i)=>{const r=i*.1;return [r,r*r,2*r,1,0];});
  const cfg={a:0,b:0,Omega:0,R:6,Z:6,separation:4,h:.5,dt:1e-8,objects:[{z:0,velocity:0,phase:0,charge:1}]};
  const r=K.makeRun(manufacturedProfile,cfg),before=r.snapshot();r.step();const after=r.snapshot();
  let err=0,points=0;
  for(let i=0;i<r.grid.Nr-2;i++)for(let j=2;j<r.grid.Nz-2;j++){
    const idx=i*r.grid.Nz+j,rho=(i+.5)*cfg.h,z=-cfg.Z+(j+.5)*cfg.h;if(Math.hypot(rho,z)>4)continue;
    err=Math.max(err,Math.abs(after.P[idx]/cfg.dt-(6-rho*rho-z*z)));points++;
  }
  const gamma=1/Math.sqrt(1-base.speed**2),d=base.separation/2,F=K.interpolateProfile(profile,gamma*d)[0];
  const coefficient=2*Math.PI*gamma*gamma*d*F*F,flux=[];
  for(const h of [.5,1/3,.25])for(const name of ["pair_phase_0","pair_phase_halfpi","pair_phase_pi_minus","pair_phase_pi","pair_phase_pi_plus"]){
    const c={...base,h,dt:.1*h,objects:K.caseObjects(name,base)};
    const rr=K.makeRun(profile,c),s=rr.snapshot(),jl=rr.grid.Nz/2-1;let current=0;
    for(let i=0;i<rr.grid.Nr;i++){
      const kl=i*rr.grid.Nz+jl,kr=kl+1,rho=(i+.5)*h;
      current+=2*Math.PI*rho*(s.U[kl]*s.W[kr]-s.W[kl]*s.U[kr]);
    }
    flux.push({h,name,phase:c.objects[1].phase,initial_discrete_rightward_Qbar_flux:current,continuum_prediction:coefficient*Math.sin(c.objects[1].phase)});
  }
  const planeQuadratures=[];
  for(const h of [.5,1/3,.25,.125,.0625]){
    let ch=0;
    for(let i=0;i<Math.round(base.R/h);i++){
      const rho=(i+.5)*h;
      const A=K.interpolateProfile(profile,Math.hypot(rho,gamma*(d-h/2)))[0];
      const B=K.interpolateProfile(profile,Math.hypot(rho,gamma*(d+h/2)))[0];
      ch+=2*Math.PI*rho*(A*A-B*B);
    }
    planeQuadratures.push({h,discrete_coefficient:ch,relative_difference_from_continuum:ch/coefficient-1});
  }
  return {kind:"Post-freeze supplementary controls; no changes to physical evolutions",
    stationary_profile_defects:rows,
    manufactured_operator_control:{note:"Uses a=b=0 only to test the operator, not as a Candidate 1 physical run",profile:"F=r^2, S=1, Omega=0",target:"u_tt=6-r^2",dt:cfg.dt,points,maximum_absolute_error:err},
    initial_current:{assumptions:"Equal charges, symmetric centers and opposite velocities, equal radial profile, supplied initial superposition",gamma,half_separation:d,F_at_gamma_d:F,continuum_coefficient:coefficient,law:"Qbar_left_dot=-2*pi*gamma^2*d*F(gamma*d)^2*sin(theta)",flux,plane_quadratures:planeQuadratures},
    pass:rows[2].relative_L2_X_stationary_acceleration_defect<rows[0].relative_L2_X_stationary_acceleration_defect&&rows[2].phi_acceleration_RMS<rows[0].phi_acceleration_RMS&&err<1e-9&&Math.abs(planeQuadratures[4].relative_difference_from_continuum)<Math.abs(planeQuadratures[0].relative_difference_from_continuum)};
}
if(require.main===module){
  const fs=require("node:fs"),path=require("node:path"),result=calculate();
  fs.mkdirSync(path.join(__dirname,"results"),{recursive:true});
  fs.writeFileSync(path.join(__dirname,"results","additional_controls.json"),JSON.stringify(result,null,2)+"\n");
  console.log(JSON.stringify(result,null,2));if(!result.pass)process.exitCode=1;
}
module.exports={calculate};
