"use strict";
function regions(frame,grid,cfg){
  const a=frame.lineDensity,threshold=.05*Math.max(...a),groups=[];let start=-1;
  for(let j=0;j<=a.length;j++){
    const yes=j<a.length&&a[j]>=threshold;
    if(yes&&start<0)start=j;
    if(!yes&&start>=0){
      let norm=0,moment=0;for(let k=start;k<j;k++){const z=-cfg.Z+(k+.5)*grid.h;norm+=a[k]*grid.h;moment+=z*a[k]*grid.h;}
      if(norm>.05*frame.summary.normF2)groups.push({z:moment/norm,norm_fraction:norm/frame.summary.normF2});start=-1;
    }
  }
  return groups.sort((a,b)=>b.norm_fraction-a.norm_fraction).slice(0,3).sort((a,b)=>a.z-b.z);
}
function ballDiagnostics(run,frame,cfg){
  const centers=regions(frame,run.grid,cfg),{Nr,Nz,h}=run.grid,[U,W,S,P,Q,T]=run.state,ih2=1/(h*h),result=[];
  for(let c=0;c<centers.length;c++)for(const radius of [8,12,16]){
    const z0=centers[c].z,r2=radius*radius;
    let energy=0,charge=0,norm=0,canonical=0;const modeNorm=run.modes.map(()=>0);
    for(let i=0;i<Math.min(Nr,Math.ceil(radius/h));i++){
      const rho=(i+.5)*h,mass=2*Math.PI*rho*h*h,rp=(i+1)/(i+.5),rm=i/(i+.5);
      const j0=Math.max(0,Math.floor((z0-radius+cfg.Z)/h)),j1=Math.min(Nz-1,Math.ceil((z0+radius+cfg.Z)/h));
      for(let j=j0;j<=j1;j++){
        const z=-cfg.Z+(j+.5)*h;if(rho*rho+(z-z0)*(z-z0)>r2)continue;
        const k=i*Nz+j,u=U[k],w=W[k],s=S[k],f2=u*u+w*w;
        const pot=.25*cfg.a*(s*s-1)*(s*s-1)+.5*s*s*f2+.25*cfg.b*f2*f2;
        let e=mass*(.5*(P[k]*P[k]+Q[k]*Q[k]+T[k]*T[k])+pot);
        for(const [other,coupling,boundary]of [[i+1<Nr?k+Nz:-1,rp,i+1===Nr],[i>0?k-Nz:k,rm,false],[j+1<Nz?k+1:-1,1,j+1===Nz],[j>0?k-1:-1,1,j===0]]){
          const du=(other>=0?U[other]:0)-u,dw=(other>=0?W[other]:0)-w,ds=(other>=0?S[other]:1)-s;
          e+=(boundary?.5:.25)*mass*ih2*coupling*(du*du+dw*dw+ds*ds);
        }
        energy+=e;charge+=mass*(w*P[k]-u*Q[k]);norm+=mass*f2;
        canonical+=mass*(f2+(s-1)*(s-1)+P[k]*P[k]+Q[k]*Q[k]+T[k]*T[k]);
        for(let m=0;m<run.modes.length;m++){
          const d=run.modes[m].d;let x=0;for(let a=0;a<6;a++)x+=d[a][k]*d[a][k];modeNorm[m]+=(run.modes[m].m===0?1:.5)*mass*x;
        }
      }
    }
    result.push({center_z:z0,radius,overlaps_other_ball:centers.some((x,i)=>i!==c&&Math.abs(x.z-z0)<2*radius),
      Ebar:energy,Qbar:charge,energy_per_absolute_charge:Math.abs(charge)>1e-12?energy/Math.abs(charge):null,
      norm_fraction:norm/frame.summary.normF2,mode_relative_norms:modeNorm.map(x=>canonical>0?Math.sqrt(x/canonical):0)});
  }
  return {centers,balls:result};
}
module.exports={regions,ballDiagnostics};
