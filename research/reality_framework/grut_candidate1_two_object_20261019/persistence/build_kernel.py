"""Generate a tangent extension without changing the archived physical kernel.

The background field/Verlet statements are copied verbatim. New modes are linear
azimuthal Fourier perturbations and do not backreact on the background.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "kernel.cjs"
source = PARENT.read_text()


def replace_once(old, new):
    global source
    if source.count(old) != 1:
        raise RuntimeError(f"Nonunique insertion anchor: {old[:70]!r}")
    source = source.replace(old, new)


helpers = r'''
  // A real cos(m*azimuth) perturbation. Each slot has a complete canonical
  // phase-space perturbation and its current acceleration. No mode backreaction.
  const modes=(cfg.modes||[]).map(spec=>({
    ...spec, d:Array.from({length:9},()=>new Float64Array(n)),
    initialNorm:0,seedScale:0
  }));
  const backgroundState=[U,W,S,P,Q,T];
  function sourceState(o,rho,z){
    const gamma=1/Math.sqrt(1-o.velocity*o.velocity),zz=gamma*(z-o.z);
    const radius=Math.hypot(rho,zz);
    if(radius>=profile[profile.length-1][0])return [0,0,0,0,0,0];
    const f=interpolateProfile(profile,radius),phase=o.phase+o.charge*omega*gamma*o.velocity*(z-o.z);
    const cs=Math.cos(phase),sn=Math.sin(phase),radialTime=radius>0?-gamma*o.velocity*zz/radius:0;
    const fdot=radialTime*f[1],freq=o.charge*omega*gamma;
    return [f[0]*cs,f[0]*sn,f[2]-1,fdot*cs+freq*f[0]*sn,
      fdot*sn-freq*f[0]*cs,radialTime*f[3]];
  }
  function initializeModes(){
    let bgNorm=0;
    for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
      const k=i*Nz+j,mass=weights[i];
      bgNorm+=mass*(U[k]*U[k]+W[k]*W[k]+(S[k]-1)*(S[k]-1)+P[k]*P[k]+Q[k]*Q[k]+T[k]*T[k]);
    }
    for(const mode of modes){
      const d=mode.d;
      if(mode.seed==='phase'){
        for(let k=0;k<n;k++){d[0][k]=-W[k];d[1][k]=U[k];d[3][k]=-Q[k];d[4][k]=P[k];}
      }else if(mode.seed==='amplitude'){
        for(let k=0;k<n;k++){d[0][k]=U[k];d[1][k]=.5*W[k];d[2][k]=.1*(S[k]-1);d[3][k]=.3*P[k];d[4][k]=.2*Q[k];d[5][k]=.1*T[k];}
      }else if(mode.seed!=='zero'){
        // Derivatives of the supplied boosted source, with fixed 1e-4 radius
        // increment. m=1 shifts the two source centers oppositely in x.
        // m=2 is infinitesimal transverse ellipticity: -rho*d/drho.
        const eps=1e-4;
        for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
          const k=i*Nz+j,rho=(i+.5)*h,z=-cfg.Z+(j+.5)*h;
          for(let o=0;o<items.length;o++){
            const x=sourceState(items[o],rho+eps,z),y=sourceState(items[o],rho-eps,z);
            const coefficient=mode.seed==='elliptic'?-rho:(mode.seed==='translation'?-1:(o===0?-1:1));
            for(let a=0;a<6;a++)d[a][k]+=coefficient*(x[a]-y[a])/(2*eps);
          }
        }
      }
      const angularFactor=mode.m===0?1:.5;
      let norm=0;
      for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
        const k=i*Nz+j;let value=0;for(let a=0;a<6;a++)value+=d[a][k]*d[a][k];
        norm+=angularFactor*weights[i]*value;
      }
      mode.seedScale=norm>0?(mode.epsilon??1e-4)*Math.sqrt(bgNorm/norm):0;
      for(let a=0;a<6;a++)for(let k=0;k<n;k++)d[a][k]*=mode.seedScale;
      mode.initialNorm=norm*mode.seedScale*mode.seedScale;
    }
  }
  function modeAcceleration(mode){
    const d=mode.d,A=d[0],B=d[1],C=d[2],AA=d[6],AB=d[7],AC=d[8],m2=mode.m*mode.m;
    for(let i=0;i<Nr;i++){
      const rp=plus[i],rm=minus[i],base=i*Nz,centrifugal=m2/((i+.5)*(i+.5)*h*h);
      for(let j=0;j<Nz;j++){
        const k=base+j,u=U[k],w=W[k],s=S[k],aa=A[k],bb=B[k],cc=C[k];
        const kp=k+Nz,km=k-Nz;
        const ap=i+1<Nr?A[kp]:0,bp=i+1<Nr?B[kp]:0,cp=i+1<Nr?C[kp]:0;
        const am=i>0?A[km]:aa,bm=i>0?B[km]:bb,cm=i>0?C[km]:cc;
        const al=j>0?A[k-1]:0,bl=j>0?B[k-1]:0,cl=j>0?C[k-1]:0;
        const ar=j+1<Nz?A[k+1]:0,br=j+1<Nz?B[k+1]:0,cr=j+1<Nz?C[k+1]:0;
        const f2=u*u+w*w,uw=2*b*u*w,su=2*s*u,sw=2*s*w;
        AA[k]=ih2*(rp*(ap-aa)+rm*(am-aa)+al+ar-2*aa)-centrifugal*aa
          -(s*s+b*(3*u*u+w*w))*aa-uw*bb-su*cc;
        AB[k]=ih2*(rp*(bp-bb)+rm*(bm-bb)+bl+br-2*bb)-centrifugal*bb
          -uw*aa-(s*s+b*(u*u+3*w*w))*bb-sw*cc;
        AC[k]=ih2*(rp*(cp-cc)+rm*(cm-cc)+cl+cr-2*cc)-centrifugal*cc
          -su*aa-sw*bb-(a*(3*s*s-1)+f2)*cc;
      }
    }
  }
  function modeDiagnostics(){
    let bgNorm=0;
    for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
      const k=i*Nz+j;bgNorm+=weights[i]*(U[k]*U[k]+W[k]*W[k]+(S[k]-1)*(S[k]-1)+P[k]*P[k]+Q[k]*Q[k]+T[k]*T[k]);
    }
    return modes.map(mode=>{
      let norm=0;const d=mode.d,angularFactor=mode.m===0?1:.5;
      for(let i=0;i<Nr;i++)for(let j=0;j<Nz;j++){
        const k=i*Nz+j;let val=0;for(let a=0;a<6;a++)val+=d[a][k]*d[a][k];norm+=angularFactor*weights[i]*val;
      }
      return {m:mode.m,seed:mode.seed,norm_squared:norm,amplification:mode.initialNorm>0?Math.sqrt(norm/mode.initialNorm):0,
        relative_to_background:bgNorm>0?Math.sqrt(norm/bgNorm):0};
    });
  }
  initializeModes();
'''

replace_once('  let steps=0;acceleration();', helpers + '\n  let steps=0;acceleration();for(const mode of modes)modeAcceleration(mode);')
old = '''      acceleration();
      for(let k=0;k<n;k++){P[k]+=half*AU[k];Q[k]+=half*AW[k];T[k]+=half*AS[k];}'''
new = '''      for(const mode of modes){const d=mode.d;for(let k=0;k<n;k++){d[3][k]+=half*d[6][k];d[4][k]+=half*d[7][k];d[5][k]+=half*d[8][k];d[0][k]+=dt*d[3][k];d[1][k]+=dt*d[4][k];d[2][k]+=dt*d[5][k];}}
      acceleration();for(const mode of modes)modeAcceleration(mode);
      for(let k=0;k<n;k++){P[k]+=half*AU[k];Q[k]+=half*AW[k];T[k]+=half*AS[k];}
      for(const mode of modes){const d=mode.d;for(let k=0;k<n;k++){d[3][k]+=half*d[6][k];d[4][k]+=half*d[7][k];d[5][k]+=half*d[8][k];}}'''
replace_once(old, new)
replace_once('return {step,diagnostics,snapshot,reverseMomenta,compareSnapshot,grid:',
    'return {step,diagnostics,snapshot,reverseMomenta,compareSnapshot,state:backgroundState,modes,modeDiagnostics,recompute:()=>{acceleration();for(const mode of modes)modeAcceleration(mode);},grid:')
source = '// Generated by build_kernel.py; background parent SHA recorded in kernel_origin.json.\n' + source
(HERE / 'kernel_modes.cjs').write_text(source)
(HERE / 'kernel_origin.json').write_text(json.dumps({
    'parent_path': '../kernel.cjs',
    'parent_sha256': hashlib.sha256(PARENT.read_bytes()).hexdigest(),
    'generated_sha256': hashlib.sha256(source.encode()).hexdigest(),
    'method': 'Insert synchronized linear Fourier tangent stages; background force and Verlet statements unchanged.'
}, indent=2) + '\n')
print('Generated kernel_modes.cjs')
