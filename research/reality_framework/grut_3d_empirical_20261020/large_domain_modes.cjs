// Additional symmetry audit of the ORIGINAL linear m=1 solver, not a second
// claim of unrestricted nonlinear 3D. Large domain removes returned waves.
"use strict";
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const here=__dirname,oldpath=path.join(here,'../grut_candidate1_two_object_20261019/persistence/kernel_modes.cjs');
const K=require(oldpath),profile=JSON.parse(fs.readFileSync(path.join(here,'../grut_candidate1_two_object_20261019/profile.json'),'utf8'));
const hash=p=>crypto.createHash('sha256').update(fs.readFileSync(p)).digest('hex');
const protocolPath=path.join(here,'LARGE_DOMAIN_PROTOCOL.json');
if(!fs.existsSync(protocolPath))fs.writeFileSync(protocolPath,JSON.stringify({local_freeze_unix:Date.now()/1000,public_preregistration:false,source_sha256:hash(__filename),parent_sha256:hash(oldpath),grids:[.5,1/3,.25],R:180,Z:192,T:160,dt:.05,sample:5,scope:'Original axisymmetric-background linear m=1 sector with full translation, full rotation/tilt and relative-shift seeds. Independent of unrestricted nonlinear Cartesian outcomes.',posthoc:'Added after smaller Cartesian late symmetry residuals showed a resolution-independent floor; tests returned-boundary-wave interpretation rather than silently accepting it.',criteria:'Report symmetry residual at every resolution/time and actual convergence. Positive norms need not remain constant. If a global symmetry fails to converge, do not label relative-shift amplification a resolved physical instability.'},null,2)+'\n');
const spec=JSON.parse(fs.readFileSync(protocolPath));if(spec.source_sha256!==hash(__filename)||spec.parent_sha256!==hash(oldpath))throw Error('frozen source mismatch');
function derivative(f,i,j,Nr,Nz,h,axis){
 const k=i*Nz+j;
 if(axis==='r'){const lo=i>0?f[k-Nz]:f[k],hi=i+1<Nr?f[k+Nz]:(f===null?1:0);return (hi-lo)/(2*h);}
 return ((j+1<Nz?f[k+1]:0)-(j>0?f[k-1]:0))/(2*h);
}
function measure(run){
 const {Nr,Nz,h}=run.grid,s=run.state,m=run.modes,scale=m.map(x=>x.seedScale),norm=[0,0,0],gn=[0,0],res=[0,0];
 let centers=[{peak:-1,z:0},{peak:-1,z:0}];
 for(let j=0;j<Nz;j++){const z=-spec.Z+(j+.5)*h,k=j,F=s[0][k]**2+s[1][k]**2,a=z>0?1:0;if(F>centers[a].peak)centers[a]={peak:F,z};}
 // Late core geometric shape: project real amplitude + S onto local x-shift.
 const core=centers.map(c=>({z:c.z,raw:0,translation_norm:0,dot:0}));
 for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
  const k=i*Nz+j,rho=(i+.5)*h,z=-spec.Z+(j+.5)*h,weight=Math.PI*rho*h*h;
  const dr=[],dz=[];
  for(let c=0;c<6;c++){
   const f=s[c],vac=c===2?1:0,lo=i>0?f[k-Nz]:f[k],hi=i+1<Nr?f[k+Nz]:vac;
   dr[c]=(hi-lo)/(2*h);dz[c]=((j+1<Nz?f[k+1]:vac)-(j>0?f[k-1]:vac))/(2*h);
   const G=[-dr[c],rho*dz[c]-z*dr[c]];
   for(let a=0;a<3;a++){const d=m[a].d[c][k]/scale[a];norm[a]+=weight*d*d;if(a<2){gn[a]+=weight*G[a]*G[a];res[a]+=weight*(d-G[a])**2;}}
  }
  const F=Math.hypot(s[0][k],s[1][k]);if(F<1e-8)continue;
  const df=(s[0][k]*m[2].d[0][k]+s[1][k]*m[2].d[1][k])/(F*scale[2]),ds=m[2].d[2][k]/scale[2];
  const tx=-(s[0][k]*dr[0]+s[1][k]*dr[1])/F,ts=-dr[2];
  for(const c of core)if(rho*rho+(z-c.z)**2<64){c.raw+=weight*(df*df+ds*ds);c.translation_norm+=weight*(tx*tx+ts*ts);c.dot+=weight*(df*tx+ds*ts);}
 }
 return {time:run.diagnostics().summary.time,norm:norm.map(Math.sqrt),generator_norm:gn.map(Math.sqrt),relative_symmetry_residual:res.map((x,a)=>Math.sqrt(x/gn[a])),late_core_shape:core.map(c=>({...c,fitted_displacement:c.dot/c.translation_norm,remaining_shape_norm:Math.sqrt(Math.max(0,c.raw-c.dot*c.dot/c.translation_norm)),fraction_explained_by_translation:c.dot*c.dot/(c.translation_norm*c.raw)}))};
}
const requested=process.argv.slice(2).map(Number);
for(const h of spec.grids){
 if(requested.length&&!requested.some(x=>Math.abs(x-h)<1e-12))continue;
 const id='large_modes_h'+h.toFixed(6),outpath=path.join(here,'raw',id+'.json');if(fs.existsSync(outpath))throw Error('do not overwrite '+id);
 const cfg={a:1,b:.1,Omega:.8000000019763136,R:spec.R,Z:spec.Z,h,dt:spec.dt,modes:[{m:1,seed:'translation',epsilon:1e-4},{m:1,seed:'zero',epsilon:1e-4},{m:1,seed:'relative_shift',epsilon:1e-4}],objects:[{z:-12,velocity:.2,phase:0,charge:1},{z:12,velocity:-.2,phase:Math.PI/2,charge:1}]};
 console.log('START',id);const begin=Date.now(),r=K.makeRun(profile,cfg),{Nr,Nz}=r.grid,d=r.modes[1].d;let norm=0,bgnorm=0;
 for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
  const k=i*Nz+j,rho=(i+.5)*h,z=-spec.Z+(j+.5)*h,w=Math.PI*rho*h*h;
  for(let c=0;c<6;c++){const f=r.state[c],vac=c===2?1:0,dr=((i+1<Nr?f[k+Nz]:vac)-(i>0?f[k-Nz]:f[k]))/(2*h),dz=((j+1<Nz?f[k+1]:vac)-(j>0?f[k-1]:vac))/(2*h);d[c][k]=rho*dz-z*dr;norm+=w*d[c][k]**2;bgnorm+=2*w*(f[k]-vac)**2;}
 }
 const scale=1e-4*Math.sqrt(bgnorm/norm);for(let c=0;c<6;c++)for(let k=0;k<Nr*Nz;k++)d[c][k]*=scale;r.modes[1].seedScale=scale;r.modes[1].initialNorm=norm*scale*scale;r.recompute();
 const trace=[measure(r)];for(let t=spec.sample;t<=spec.T+1e-8;t+=spec.sample){r.step(Math.round(spec.sample/spec.dt));trace.push(measure(r));fs.writeFileSync(outpath,JSON.stringify({id,cfg,source_sha256:hash(__filename),parent_sha256:hash(oldpath),wall_seconds:(Date.now()-begin)/1000,trace},null,2)+'\n');}
 console.log('END',id,'seconds',(Date.now()-begin)/1000,JSON.stringify(trace.at(-1)));
}
