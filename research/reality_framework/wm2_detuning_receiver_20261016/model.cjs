"use strict";
// Dependency-free numerical model. The all-N proofs are in MATHEMATICS.md.
function createS3() {
  const pars={U:16,Delta:4,Lambda:64,epsilon:.25,J:.25};
  function assert(x,message="assertion failed"){if(!x)throw Error(message);}
  function bits(q){let n=0;for(;q;q>>>=1)n+=q&1;return n;}
  function mask(w){return w.reduce((a,v)=>a|(1<<v),0);}
  function words(n){const out=[];function go(w){out.push(w);for(let v=0;v<n;v++)if(!w.includes(v))go(w.concat(v));}go([0]);return out;}
  function occupations(k,L){const out=[];function go(a,left){if(a.length===L-1){out.push(a.concat(left));return;}for(let x=0;x<=left;x++)go(a.concat(x),left-x);}go([],k);return out;}
  function key(s){return s.w.join(",")+"|"+s.q+"|"+s.o.join(",")+"|"+s.p;}
  function mismatch(n,w,q){return bits(mask(w)^q);}
  function diag(n,s,delta){return pars.U*(n-s.w.length+bits(s.q)-n)+pars.Lambda*mismatch(n,s.w,s.q)+delta*(bits(s.q)-1);}
  function neighbors(n,L,s,{delta=0,g=1}={}) {
    const out=[],w=s.w,u=w[w.length-1],p=s.p;
    const add=(ww,q,o,pp,kind,a)=>out.push({s:{w:ww,q,o,p:pp},kind,a});
    if(w.length<n)for(let v=0;v<n;v++)if(!w.includes(v)){
      if(!(s.q>>v&1)){const o=s.o.slice();o[0]++;add(w.concat(v),s.q|(1<<v),o,p,"emission_append",-g/n*Math.sqrt(o[0]));}
      if(p!==null&&(p===u||p===v))add(w.concat(v),s.q,s.o,p,"bare_append",-pars.epsilon/Math.sqrt(n));
    }
    if(w.length>1){
      const v=u,parent=w[w.length-2];
      if((s.q>>v&1)&&s.o[0]){const o=s.o.slice();o[0]--;add(w.slice(0,-1),s.q^(1<<v),o,p,"emission_pop",-g/n*Math.sqrt(s.o[0]));}
      if(p!==null&&(p===parent||p===v))add(w.slice(0,-1),s.q,s.o,p,"bare_pop",-pars.epsilon/Math.sqrt(n));
    }
    for(let r=0;r+1<L;r++){
      if(s.o[r]){const o=s.o.slice();o[r]--;o[r+1]++;add(w,s.q,o,p,"receiver",-Math.sqrt(s.o[r]*(s.o[r+1]+1)));}
      if(s.o[r+1]){const o=s.o.slice();o[r+1]--;o[r]++;add(w,s.q,o,p,"receiver",-Math.sqrt(s.o[r+1]*(s.o[r]+1)));}
    }
    if(p!==null){const k=w.indexOf(p);if(k>0)add(w,s.q,s.o,w[k-1],"matter",-pars.J);if(k>=0&&k+1<w.length)add(w,s.q,s.o,w[k+1],"matter",-pars.J);}
    return out;
  }
  function make(n,L,{probe=true,delta=0,g=1}={}) {
    assert(n>=2&&n<=8&&L>=1,"enumeration budget N=2..8; analytic results are all N");
    const occ=Array.from({length:n},(_,k)=>occupations(k,L)),states=[];
    for(const w of words(n))for(const q of (probe?Array.from({length:1<<(n-1)},(_,i)=>2*i+1):[mask(w)]))
      for(const o of occ[bits(q)-1])for(const p of (probe?Array.from({length:n},(_,i)=>i):[null]))states.push({w,q,o,p});
    const index=new Map(states.map((s,i)=>[key(s),i])),cols=[],P=[],Q=[];
    let norm1=0;
    states.forEach((s,j)=>{
      const col=[[j,diag(n,s,delta)]];
      for(const v of neighbors(n,L,s,{delta,g})){const i=index.get(key(v.s));assert(i!==undefined,"unretained transition");col.push([i,v.a,v.kind]);}
      cols.push(col);norm1=Math.max(norm1,col.reduce((a,x)=>a+Math.abs(x[1]),0));
      (s.q===mask(s.w)?P:Q).push(j);
    });
    return {n,L,delta,g,states,index,cols,P,Q,norm1,size:states.length};
  }
  function oneParticle(n,L,delta=0,g=1){
    const size=L+1,cols=Array.from({length:size},(_,j)=>[[j,j?delta:0]]);
    for(let j=0;j<L;j++){const a=j===0?-g/n:-1;cols[j].push([j+1,a]);cols[j+1].push([j,a]);}
    let norm1=0;for(const c of cols)norm1=Math.max(norm1,c.reduce((a,x)=>a+Math.abs(x[1]),0));return {size,cols,norm1};
  }
  function multiply(h,x){const y=new Float64Array(h.size);for(let j=0;j<h.size;j++)if(x[j])for(const [i,a] of h.cols[j])y[i]+=a*x[j];return y;}
  function norm2(r,i){let s=0;for(let k=0;k<r.length;k++)s+=r[k]*r[k]+i[k]*i[k];return s;}
  function evolve(h,initial,t){
    let re=Float64Array.from(initial.re),im=Float64Array.from(initial.im||new Float64Array(h.size));
    const steps=Math.max(1,Math.ceil(Math.abs(t)*h.norm1/.8)),dt=t/steps;
    for(let j=0;j<steps;j++){
      let tr=re.slice(),ti=im.slice(),sr=re.slice(),si=im.slice();
      for(let k=1;k<=35;k++){
        const hr=multiply(h,tr),hi=multiply(h,ti),a=dt/k;
        for(let b=0;b<h.size;b++){tr[b]=a*hi[b];ti[b]=-a*hr[b];sr[b]+=tr[b];si[b]+=ti[b];}
        if(k>8&&norm2(tr,ti)<1e-33)break;
        assert(k<35,"Taylor stage did not converge");
      }
      re=sr;im=si;
    }
    return {re,im};
  }
  function spectral(n,delta,g,t,refine=1){
    const l=g/n,l2=l*l,c=1-l2,minimum=l2*l2*(1-delta*delta/(4*c));
    const steps=Math.ceil(512*n*n/(g*g))*refine,h=Math.PI/steps;
    let total=0,re=0,im=0;
    for(let j=1;j<steps;j++){
      const th=j*h,s=Math.sin(th),E=delta+2*Math.cos(th);
      const D=c*(E+l2*delta/(2*c))**2+minimum;
      const w=2*l2/Math.PI*s*s/D*h;
      total+=w;re+=w*Math.cos(E*t);im-=w*Math.sin(E*t);
    }
    return {re,im,total_density_mass:total,theta_intervals:steps,survival:re*re+im*im};
  }
  function binomialMacro(n,survival){
    const M=n-1,L=Math.max(2,Math.ceil(n/2)),p=1-survival;
    if(p<=0)return 0;if(p>=1)return 1;
    let probability=0;
    for(let k=L-1;k<=M;k++){let log=0;for(let j=1;j<=k;j++)log+=Math.log((M-j+1)/j);log+=k*Math.log(p)+(M-k)*Math.log(survival);probability+=Math.exp(log);}
    return Math.min(1,probability);
  }
  function recordEntropy(n,survival){
    const M=n-1,p=1-survival;let entropy=Math.log2(n);
    for(let k=0;k<=M;k++){
      let lc=0,lp=0;for(let j=1;j<=k;j++){lc+=Math.log((M-j+1)/j);lp+=Math.log(M-j+1);}
      const log=lc+(k? k*Math.log(p):0)+(M-k? (M-k)*Math.log(survival):0);
      const x=Math.exp(log);if(x>0)entropy+=x*(-Math.log2(x)+lp/Math.LN2);
    }return entropy;
  }
  function tail(d,z){if(d===0)return 1;if(z===0)return 0;let x=1;for(let j=1;j<=d;j++)x*=z/j;let sum=x;for(let j=d+1;j<10000;j++){x*=z/j;sum+=x;if(x<1e-17*Math.max(sum,1e-300))break;}return Math.min(1,sum);}
  function constants(){
    const A=1000/(729*Math.PI),near=.25*.005-Math.expm1(4.2*.005)+4.2*.005;
    return {U:16,Delta:4,Lambda:64,epsilon:.25,J:.25,delta_interval:[-.5,.5],g_interval:[.9,1.1],
      receiver_quanta_max:"N-1 (exact charge)",finite_receiver_modes:"64*N^2",observation_window:"[4*N^2,8*N^2]",
      density_variation_C:40/(9*Math.PI),continuum_amplitude_bound:A,
      conservative_formation_probability_lower:.5,sharper_infinite_receiver_lower:1-2*A*A,
      finite_return_amplitude_error_upper:"2*exp(-120*N^2), N>=2",
      all_geometry_probe_offdiagonal_norm_upper:2.1,near_full_H_response_lower_at_t_005:near,
      prefix_crossing_norm_upper:"2.7/sqrt(N)",
      complex_Schur_bound:"(1/(8*eta))*(1+5*eta/21)^(-d), eta>0",
      no_real_energy_Schur_claim:true};
  }
  function register(s,n){
    return {edges:s.w.slice(0,-1).map((u,k)=>[u,s.w[k+1]]),head:s.w.at(-1),
      free:Array.from({length:n},(_,i)=>i).filter(i=>!s.w.includes(i)),q:s.q,o:s.o,p:s.p};
  }
  function rkey(s){return s.edges.map(x=>x.join(",")).sort().join(";")+"|"+s.head+"|"+s.free.slice().sort((a,b)=>a-b).join(",")+"|"+s.q+"|"+s.o.join(",")+"|"+s.p;}
  function registerNeighbors(n,L,c,{g=1}={}){
    const out=[];
    const add=(edges,head,free,q,o,p,kind,a)=>out.push({s:{edges,head,free,q,o,p},kind,a});
    for(let u=0;u<n;u++)for(let v=0;v<n;v++)if(u!==v){
      if(c.head===u&&c.free.includes(v)&&!c.edges.some(b=>b[0]===u&&b[1]===v)){
        const edges=c.edges.concat([[u,v]]),free=c.free.filter(x=>x!==v);
        if(!(c.q>>v&1)){const o=c.o.slice();o[0]++;add(edges,v,free,c.q|(1<<v),o,c.p,"emission_append",-g/n*Math.sqrt(o[0]));}
        if(c.p!==null&&(c.p===u||c.p===v))add(edges,v,free,c.q,c.o,c.p,"bare_append",-pars.epsilon/Math.sqrt(n));
      }
      if(c.head===v&&c.edges.some(b=>b[0]===u&&b[1]===v)){
        const edges=c.edges.filter(b=>!(b[0]===u&&b[1]===v)),free=c.free.concat(v);
        if((c.q>>v&1)&&c.o[0]){const o=c.o.slice();o[0]--;add(edges,u,free,c.q^(1<<v),o,c.p,"emission_pop",-g/n*Math.sqrt(c.o[0]));}
        if(c.p!==null&&(c.p===u||c.p===v))add(edges,u,free,c.q,c.o,c.p,"bare_pop",-pars.epsilon/Math.sqrt(n));
      }
    }
    for(let r=0;r+1<L;r++)for(const direction of [1,-1]){
      const from=direction===1?r:r+1,to=direction===1?r+1:r;
      if(c.o[from]){const o=c.o.slice();o[from]--;o[to]++;add(c.edges,c.head,c.free,c.q,o,c.p,"receiver",-Math.sqrt(c.o[from]*(c.o[to]+1)));}
    }
    for(const [u,v] of c.edges)if(c.p===u)add(c.edges,c.head,c.free,c.q,c.o,v,"matter",-pars.J);
      else if(c.p===v)add(c.edges,c.head,c.free,c.q,c.o,u,"matter",-pars.J);
    return out;
  }
  return {assert,bits,mask,words,occupations,key,mismatch,diag,neighbors,make,oneParticle,multiply,norm2,evolve,spectral,binomialMacro,recordEntropy,tail,constants,register,rkey,registerNeighbors};
}
if(typeof module!=="undefined"&&module.exports)module.exports=createS3();
