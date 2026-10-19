"use strict";
/* Axisymmetric finite-volume Hamiltonian evolution of Candidate 1.
   Coordinates x=g*v*r, tau=g*v*t. X=v*(u+i*w)/sqrt(2), S=phi/v.
   D=0, a=lambda_phi/g^2=1, b=lambda_X/g^2=0.1.
   Radial cell centers enforce regularity by zero axis flux.
   Vacuum ghost values at outer edges. Matched weighted bond energy.
   Kick-drift-kick Verlet; U(1) charge is conserved to roundoff. */
function interpolateProfile(rows,r) {
  if(r>=rows[rows.length-1][0]) return [0,0,1,0];
  const h=rows[1][0]-rows[0][0];
  let k=Math.min(rows.length-2,Math.max(0,Math.floor(r/h)));
  const p=rows[k],q=rows[k+1],d=q[0]-p[0],t=(r-p[0])/d;
  const t2=t*t,t3=t2*t;
  const h00=2*t3-3*t2+1,h10=t3-2*t2+t,h01=-2*t3+3*t2,h11=t3-t2;
  const d00=(6*t2-6*t)/d,d10=3*t2-4*t+1,d01=(-6*t2+6*t)/d,d11=3*t2-2*t;
  return [h00*p[1]+h10*d*p[2]+h01*q[1]+h11*d*q[2],
    d00*p[1]+d10*p[2]+d01*q[1]+d11*q[2],
    h00*p[3]+h10*d*p[4]+h01*q[3]+h11*d*q[4],
    d00*p[3]+d10*p[4]+d01*q[3]+d11*q[4]];
}
function makeRun(profile,cfg) {
  const h=cfg.h,dt=cfg.dt,Nr=Math.round(cfg.R/h),Nz=Math.round(2*cfg.Z/h),n=Nr*Nz;
  if(Math.abs(Nr*h-cfg.R)>1e-9||Math.abs(Nz*h-2*cfg.Z)>1e-9) throw Error("Nonintegral domain");
  const U=new Float64Array(n),W=new Float64Array(n),S=new Float64Array(n);
  const P=new Float64Array(n),Q=new Float64Array(n),T=new Float64Array(n);
  const AU=new Float64Array(n),AW=new Float64Array(n),AS=new Float64Array(n);
  const weights=new Float64Array(Nr),plus=new Float64Array(Nr),minus=new Float64Array(Nr);
  for(let i=0;i<Nr;i++){weights[i]=2*Math.PI*(i+.5)*h*h*h;plus[i]=(i+1)/(i+.5);minus[i]=i/(i+.5);}
  S.fill(1);
  const omega=cfg.Omega,a=cfg.a,b=cfg.b;
  const items=cfg.objects||[];
  for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
    const idx=i*Nz+j,rho=(i+.5)*h,z=-cfg.Z+(j+.5)*h;
    for(const o of items){
      const gamma=1/Math.sqrt(1-o.velocity*o.velocity),zz=gamma*(z-o.z);
      const r=Math.hypot(rho,zz),f=interpolateProfile(profile,r);
      const phase=o.phase+o.charge*omega*gamma*o.velocity*(z-o.z);
      const cs=Math.cos(phase),sn=Math.sin(phase);
      const radialTime= r>0 ? -gamma*o.velocity*zz/r:0;
      const fdot=radialTime*f[1],sdot=radialTime*f[3],freq=o.charge*omega*gamma;
      U[idx]+=f[0]*cs;W[idx]+=f[0]*sn;S[idx]+=f[2]-1;
      P[idx]+=fdot*cs+freq*f[0]*sn;Q[idx]+=fdot*sn-freq*f[0]*cs;T[idx]+=sdot;
    }
  }
  const ih2=1/(h*h);
  function acceleration(){
    for(let i=0;i<Nr;i++){
      const rp=plus[i],rm=minus[i],base=i*Nz;
      for(let j=0;j<Nz;j++){
        const k=base+j,u=U[k],w=W[k],s=S[k];
        const kp=k+Nz,km=k-Nz;
        const up=i+1<Nr?U[kp]:0,wp=i+1<Nr?W[kp]:0,sp=i+1<Nr?S[kp]:1;
        const um=i>0?U[km]:u,wm=i>0?W[km]:w,sm=i>0?S[km]:s;
        const ul=j>0?U[k-1]:0,wl=j>0?W[k-1]:0,sl=j>0?S[k-1]:1;
        const ur=j+1<Nz?U[k+1]:0,wr=j+1<Nz?W[k+1]:0,sr=j+1<Nz?S[k+1]:1;
        const f2=u*u+w*w,c=s*s+b*f2;
        AU[k]=ih2*(rp*(up-u)+rm*(um-u)+ul+ur-2*u)-c*u;
        AW[k]=ih2*(rp*(wp-w)+rm*(wm-w)+wl+wr-2*w)-c*w;
        AS[k]=ih2*(rp*(sp-s)+rm*(sm-s)+sl+sr-2*s)-(a*(s*s-1)+f2)*s;
      }
    }
  }
  let steps=0;acceleration();
  function step(count=1){
    const half=.5*dt;
    for(let m=0;m<count;m++){
      for(let k=0;k<n;k++){P[k]+=half*AU[k];Q[k]+=half*AW[k];T[k]+=half*AS[k];U[k]+=dt*P[k];W[k]+=dt*Q[k];S[k]+=dt*T[k];}
      acceleration();
      for(let k=0;k<n;k++){P[k]+=half*AU[k];Q[k]+=half*AW[k];T[k]+=half*AS[k];}
      steps++;
    }
  }
  function diagnostics(){
    let energy=0,charge=0,absoluteCharge=0,number=0,zNum=0,z2Num=0,r2Num=0,midNumber=0;
    let leftCharge=0,rightCharge=0,leftNumber=0,rightNumber=0;
    let maxF=0,minS=Infinity;
    const lineDensity=Array(Nz).fill(0),lineCharge=Array(Nz).fill(0),axisF2=Array(Nz).fill(0),axisS=Array(Nz).fill(0);
    for(let i=0;i<Nr;i++){
      const mass=weights[i],base=i*Nz,rho=(i+.5)*h;
      for(let j=0;j<Nz;j++){
        const k=base+j,u=U[k],w=W[k],s=S[k],f2=u*u+w*w,z=-cfg.Z+(j+.5)*h;
        const density=w*P[k]-u*Q[k],kin=.5*(P[k]*P[k]+Q[k]*Q[k]+T[k]*T[k]);
        const pot=.25*a*(s*s-1)*(s*s-1)+.5*s*s*f2+.25*b*f2*f2;
        energy+=mass*(kin+pot);charge+=mass*density;absoluteCharge+=mass*Math.abs(density);
        number+=mass*f2;zNum+=mass*f2*z;z2Num+=mass*f2*z*z;r2Num+=mass*f2*rho*rho;
        if(Math.abs(z)<3)midNumber+=mass*f2;
        if(z<0){leftCharge+=mass*density;leftNumber+=mass*f2;}else{rightCharge+=mass*density;rightNumber+=mass*f2;}
        lineDensity[j]+=mass*f2/h;lineCharge[j]+=mass*density/h;
        if(i===0){axisF2[j]=f2;axisS[j]=s;}
        maxF=Math.max(maxF,Math.sqrt(f2));minS=Math.min(minS,s);
        // Positive-direction interior bonds count each edge once.
        // Outer vacuum ghost bonds and lower-z ghost are included once.
        if(i+1<Nr){
          const du=U[k+Nz]-u,dw=W[k+Nz]-w,ds=S[k+Nz]-s;
          energy+=.5*mass*ih2*plus[i]*(du*du+dw*dw+ds*ds);
        }else energy+=.5*mass*ih2*plus[i]*(f2+(s-1)*(s-1));
        if(j+1<Nz){
          const du=U[k+1]-u,dw=W[k+1]-w,ds=S[k+1]-s;
          energy+=.5*mass*ih2*(du*du+dw*dw+ds*ds);
        }else energy+=.5*mass*ih2*(f2+(s-1)*(s-1));
        if(j===0)energy+=.5*mass*ih2*(f2+(s-1)*(s-1));
      }
    }
    const result={time:steps*dt,Ebar:energy,Qbar:charge,absQbar:absoluteCharge,
      normF2:number,z_mean:zNum/number,z_rms:Math.sqrt(z2Num/number),rho_rms:Math.sqrt(r2Num/number),
      middle_norm_fraction:midNumber/number,left_Qbar:leftCharge,right_Qbar:rightCharge,
      max_F:maxF,min_S:minS};
    if(number===0){result.z_mean=0;result.z_rms=0;result.rho_rms=0;result.middle_norm_fraction=0;}
    return {summary:result,lineDensity,lineCharge,axisF2,axisS};
  }
  function snapshot(){return {U:Array.from(U),W:Array.from(W),S:Array.from(S),P:Array.from(P),Q:Array.from(Q),T:Array.from(T),steps};}
  function reverseMomenta(){for(let k=0;k<n;k++){P[k]=-P[k];Q[k]=-Q[k];T[k]=-T[k];}}
  function compareSnapshot(snap,reverse=false){
    let maxField=0,maxMomentum=0;
    for(let k=0;k<n;k++){
      maxField=Math.max(maxField,Math.abs(U[k]-snap.U[k]),Math.abs(W[k]-snap.W[k]),Math.abs(S[k]-snap.S[k]));
      const sign=reverse?-1:1;
      maxMomentum=Math.max(maxMomentum,Math.abs(P[k]-sign*snap.P[k]),Math.abs(Q[k]-sign*snap.Q[k]),Math.abs(T[k]-sign*snap.T[k]));
    }
    return {maximum_field_error:maxField,maximum_momentum_error:maxMomentum};
  }
  return {step,diagnostics,snapshot,reverseMomenta,compareSnapshot,grid:{Nr,Nz,points:n,h,dt},configuration:cfg};
}
function caseObjects(name,c){
  const speed=c.speed,half=c.separation/2;
  if(name==="vacuum")return [];
  if(name==="single_static")return [{z:0,velocity:0,phase:0,charge:1}];
  if(name==="single_boost")return [{z:-half,velocity:speed,phase:0,charge:1}];
  const phase=name==="pair_phase_pi"?Math.PI:(name==="pair_phase_halfpi"?Math.PI/2:(name==="pair_phase_pi_minus"?Math.PI-.2:(name==="pair_phase_pi_plus"?Math.PI+.2:0)));
  return [{z:-half,velocity:speed,phase:0,charge:1},
    {z:half,velocity:-speed,phase,charge:name==="pair_opposite"?-1:1}];
}
function analyzeRun(rows,configuration){
  const first=rows[0],qscale=Math.max(first.absQbar,1e-300),e=first.Ebar;
  let edrift=0,qdrift=0,finite=true;
  for(const r of rows){edrift=Math.max(edrift,Math.abs(r.Ebar-e)/Math.max(Math.abs(e),1e-300));
    qdrift=Math.max(qdrift,Math.abs(r.Qbar-first.Qbar)/qscale);
    finite=finite&&Object.values(r).every(Number.isFinite);}
  return {finite,maximum_relative_energy_drift:edrift,
    maximum_charge_error_over_initial_absolute_charge:qdrift,initial:first,final:rows[rows.length-1]};
}
if(typeof module!=="undefined")module.exports={interpolateProfile,makeRun,caseObjects,analyzeRun};
