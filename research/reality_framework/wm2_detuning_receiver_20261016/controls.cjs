"use strict";
function createControls(S) {
  const {assert,key,mask,bits,make,neighbors,register,registerNeighbors,rkey,evolve,norm2}=S;
  function unit(size,i){const re=new Float64Array(size);re[i]=1;return {re,im:new Float64Array(size)};}
  function sourceIndex(h,w,q,o,p){const i=h.index.get(key({w,q,o,p}));assert(i!==undefined);return i;}
  function factorial(n){let a=1;for(let k=2;k<=n;k++)a*=k;return a;}
  function complexPower(r,i,n){let a=1,b=0;for(let k=0;k<n;k++){const x=a*r-b*i;b=a*i+b*r;a=x;}return [a,b];}
  function fockAmplitude(z,o){
    const total=o.reduce((a,b)=>a+b,0);let factor=factorial(total),r=1,i=0;
    for(let j=0;j<o.length;j++){factor/=factorial(o[j]);const [a,b]=complexPower(z.re[j],z.im[j],o[j]);const x=r*a-i*b;i=r*b+i*a;r=x;}
    const a=Math.sqrt(factor);return [r*a,i*a];
  }
  function predicted(h,z){
    const M=h.n-1,re=new Float64Array(h.size),im=new Float64Array(h.size);
    for(let i=0;i<h.size;i++){
      const s=h.states[i],k=s.w.length-1,P=factorial(M)/factorial(M-k);
      const [a,b]=fockAmplitude(z,[M-k,...s.o]);re[i]=a/Math.sqrt(P);im[i]=b/Math.sqrt(P);
    }
    return {re,im};
  }
  function graph(n,L){
    const h=make(n,L),seed=sourceIndex(h,[0],1,Array(L).fill(0),0),seen=new Set([seed]),todo=[seed];
    let registerChecks=0,reverseChecks=0,covarianceChecks=0,emissionNorm=0,crossNorm=0;
    for(let a=0;a<todo.length;a++){
      const id=todo[a],s=h.states[id],ns=neighbors(n,L,s);
      const tuples=ns.map(v=>[rkey(register(v.s,n))+"|"+v.kind,v.a]).sort((a,b)=>a[0].localeCompare(b[0]));
      const raw=registerNeighbors(n,L,register(s,n)).map(v=>[rkey(v.s)+"|"+v.kind,v.a]).sort((a,b)=>a[0].localeCompare(b[0]));
      assert(tuples.length===raw.length);
      for(let j=0;j<raw.length;j++)assert(tuples[j][0]===raw[j][0]&&Math.abs(tuples[j][1]-raw[j][1])<1e-14);
      registerChecks++;
      const perm=s=>({w:s.w.map(v=>(v+1)%n),q:Array.from({length:n},(_,v)=>v).reduce((q,v)=>q|((s.q>>v&1)<<((v+1)%n)),0),o:s.o,p:s.p===null?null:(s.p+1)%n});
      const pp=ns.map(v=>[key(perm(v.s))+"|"+v.kind,v.a]).sort((a,b)=>a[0].localeCompare(b[0]));
      const qq=neighbors(n,L,perm(s)).map(v=>[key(v.s)+"|"+v.kind,v.a]).sort((a,b)=>a[0].localeCompare(b[0]));
      assert(pp.length===qq.length);for(let j=0;j<pp.length;j++)assert(pp[j][0]===qq[j][0]&&Math.abs(pp[j][1]-qq[j][1])<1e-14);
      covarianceChecks++;
      let ep=0,bp=0;
      for(const v of ns){
        const target=h.index.get(key(v.s));assert(target!==undefined);
        assert(v.s.o.reduce((a,b)=>a+b,0)===bits(v.s.q)-1);
        const inverse=neighbors(n,L,v.s).find(x=>key(x.s)===key(s)&&Math.abs(x.a-v.a)<1e-14);assert(inverse);reverseChecks++;
        if(v.kind==="emission_append")ep+=v.a*v.a;
        if(s.q===mask(s.w)&&v.s.q!==mask(v.s.w))bp+=v.a*v.a;
        if(!seen.has(target)){seen.add(target);todo.push(target);}
      }
      emissionNorm=Math.max(emissionNorm,ep);crossNorm=Math.max(crossNorm,bp);
    }
    assert(seen.size===h.size);assert(emissionNorm<=.25+1e-14&&crossNorm<=.0625+1e-14);
    const inP=new Set(h.P);let exitChecks=0,sourceCount=0,minimum=Infinity;
    for(const start of h.P){
      const s=h.states[start];if(!s.w.includes(s.p))continue;sourceCount++;
      const dist=Array(h.size).fill(Infinity),queue=[];
      for(const v of neighbors(n,L,s))if(v.s.q!==mask(v.s.w)){const i=h.index.get(key(v.s));dist[i]=1;queue.push(i);}
      for(let b=0;b<queue.length;b++){
        const id=queue[b],c=h.states[id],old=dist[id];
        for(const v of neighbors(n,L,c)){
          const j=h.index.get(key(v.s)),d=old+(v.kind==="receiver"?0:1);
          if(inP.has(j)){
            if(s.w.includes(v.s.p)){const intrinsic=Math.abs(s.w.indexOf(s.p)-s.w.indexOf(v.s.p));assert(d>=intrinsic+2);minimum=Math.min(minimum,d-intrinsic);exitChecks++;}
          }else if(d<dist[j]){dist[j]=d;queue.push(j);}
        }
      }
      const pd=Array(h.size).fill(Infinity);pd[start]=0;const pq=[start],event=Array(n).fill(Infinity);event[s.p]=0;
      for(let b=0;b<pq.length;b++){
        const id=pq[b];for(const v of neighbors(n,L,h.states[id])){
          if(v.s.q!==mask(v.s.w))continue;
          const j=h.index.get(key(v.s)),d=pd[id]+(v.kind==="receiver"?0:1);
          if(d<pd[j]){pd[j]=d;pq.push(j);event[v.s.p]=Math.min(event[v.s.p],d);}
        }
      }
      for(const j of s.w)assert(event[j]===Math.abs(s.w.indexOf(s.p)-s.w.indexOf(j)));
    }
    assert(minimum===2);
    return {N:n,receiver_modes:L,full_one_probe_dimension:h.size,P_dimension:h.P.length,Q_dimension:h.Q.length,
      reachable_basis_states:seen.size,independent_register_states:registerChecks,reverse_primitive_checks:reverseChecks,
      root_moving_covariance_checks:covarianceChecks,emission_append_norm_squared_upper_checked:emissionNorm,
      PQ_norm_squared_upper_checked:crossNorm,source_states:sourceCount,virtual_exit_distance_checks:exitChecks,
      minimum_virtual_word_length_minus_distance:minimum,receiver_moves_resummed_at_zero_cost:true,
      metric_matches_path_including_all_protected_geometry_moves:true};
  }
  function reduction(n,L,delta,g){
    const h=make(n,L,{probe:false,delta,g}),one=S.oneParticle(n,L,delta,g),start=unit(h.size,sourceIndex(h,[0],1,Array(L).fill(0),null));
    const rows=[];
    for(const t of [.25,1,4]){
      const x=evolve(h,start,t),z=evolve(one,unit(one.size,0),t),pred=predicted(h,z);
      let error=0;const pop=Array(n).fill(0);
      for(let j=0;j<h.size;j++){error+=(x.re[j]-pred.re[j])**2+(x.im[j]-pred.im[j])**2;pop[h.states[j].w.length-1]+=x.re[j]**2+x.im[j]**2;}
      error=Math.sqrt(error);assert(error<2e-10);
      const s=z.re[0]**2+z.im[0]**2;
      const success=pop.slice(Math.max(2,Math.ceil(n/2))-1).reduce((a,b)=>a+b,0),bin=S.binomialMacro(n,s);
      assert(Math.abs(success-bin)<2e-10);assert(Math.abs(norm2(x.re,x.im)-1)<2e-10);
      rows.push({t,full_wavefunction_error:error,full_norm:norm2(x.re,x.im),macroscopic_probability:success,binomial_probability:bin});
    }
    return {N:n,receiver_modes:L,delta,g,full_vacuum_configuration_dimension:h.size,rows};
  }
  function continuum(n,delta,g){
    const rows=[];for(const a of [4,6,8]){
      const t=a*n*n,z=S.spectral(n,delta,g,t);assert(Math.abs(z.total_density_mass-1)<2e-8);
      assert(Math.sqrt(z.survival)<=S.constants().continuum_amplitude_bound+1e-9);
      const p=S.binomialMacro(n,z.survival);assert(p>=.5);
      rows.push({t,pump_survival:z.survival,macroscopic_probability:p,density_mass:z.total_density_mass,theta_intervals:z.theta_intervals,complete_record_entropy_bits:S.recordEntropy(n,z.survival)});
    }
    return {N:n,delta,g,rows};
  }
  function finiteReceiver(n,delta,g){
    const L=64*n*n,one=S.oneParticle(n,L,delta,g),rows=[];let t0=0,z=unit(one.size,0);
    for(const a of [4,6,8]){
      const t=a*n*n;z=evolve(one,z,t-t0);t0=t;
      const sp=S.spectral(n,delta,g,t);const error=Math.hypot(z.re[0]-sp.re,z.im[0]-sp.im);
      assert(error<2e-8);assert(Math.abs(norm2(z.re,z.im)-1)<1e-9);
      const survival=z.re[0]**2+z.im[0]**2,p=S.binomialMacro(n,survival);assert(p>=.5);
      const M=n-1,k=M*(1-survival);let hr=0;
      for(let j=1;j<L;j++)hr-=2*(z.re[j]*z.re[j+1]+z.im[j]*z.im[j+1]);
      const hint=-2*(g/n)*(z.re[0]*z.re[1]+z.im[0]*z.im[1])*M,relativeReceiver=M*hr;
      const storedGain=(12+delta)*k,remainingBinding=16*(M-k),receiverPhysical=4*k+relativeReceiver;
      const ledger=16*M-remainingBinding-storedGain-receiverPhysical-hint;
      assert(Math.abs(ledger)<2e-7);
      rows.push({t,pump_survival:survival,macroscopic_probability:p,
        finite_vs_continuum_amplitude_error:error,unitarity_error:Math.abs(norm2(z.re,z.im)-1),
        mean_attachments:k,remaining_binding_energy:remainingBinding,storage_gain:storedGain,
        physical_receiver_energy:receiverPhysical,interaction_energy:hint,complete_energy_ledger_error:ledger,
        mean_geometry_record_work:-hint,complete_geometry_record_entropy_bits:S.recordEntropy(n,survival)});
    }
    return {N:n,receiver_modes:L,delta,g,one_particle_dimension:L+1,initial_receiver:"vacuum",rows};
  }
  function weyl(n,delta){
    const L=64*n*n,R=L/2,M=n-1,e=-delta-80/M;assert(Math.abs(e)<2);
    const c=-e/2,u=Math.sin(Math.PI/(2*(R+1)))**2,a2=8*u*(1-u)/(R+1);
    const r1=16*u*((1-c*c)+(2*c*c-1)*u)+a2*(4*c*c-3),bias=4*c*u;
    const residual=Math.sqrt(M*r1+M*(M-1)*bias*bias),bound=2*Math.PI*M/(R+1);assert(residual<=bound);
    return {N:n,delta,receiver_modes:L,box_modes:R,one_particle_target_energy:e,
      Q_residual_norm:residual,analytic_Q_residual_upper:bound,
      inverse_norm_lower_if_invertible:1/bound,coupled_real_pole_window_radius:2*Math.sqrt(n)*bound,contact_occupation_exactly_zero:true,
      real_gap_uniform_in_N:false};
  }

  function complexSchur(h,E,eta){
    const qi=h.Q,pi=h.P,n=qi.length,m=pi.length,Qmap=new Map(qi.map((x,i)=>[x,i])),Pmap=new Map(pi.map((x,i)=>[x,i]));
    const ar=new Float64Array(n*n),ai=new Float64Array(n*n),br=new Float64Array(n*m),bi=new Float64Array(n*m);
    for(let j=0;j<n;j++)for(const [i,a] of h.cols[qi[j]]){const row=Qmap.get(i);if(row!==undefined)ar[row*n+j]+=a;}
    for(let j=0;j<m;j++)for(const [i,a] of h.cols[pi[j]]){const row=Qmap.get(i);if(row!==undefined)br[row*m+j]+=a;}
    for(let j=0;j<n;j++){ar[j*n+j]-=E;ai[j*n+j]-=eta;}
    const rhs=br.slice(),origR=ar.slice(),origI=ai.slice();
    for(let k=0;k<n;k++){
      let pivot=k,max=0;for(let i=k;i<n;i++){const v=ar[i*n+k]**2+ai[i*n+k]**2;if(v>max){max=v;pivot=i;}}
      assert(max>1e-24);
      if(pivot!==k){
        for(let j=0;j<n;j++){let a=ar[k*n+j];ar[k*n+j]=ar[pivot*n+j];ar[pivot*n+j]=a;a=ai[k*n+j];ai[k*n+j]=ai[pivot*n+j];ai[pivot*n+j]=a;}
        for(let j=0;j<m;j++){let a=br[k*m+j];br[k*m+j]=br[pivot*m+j];br[pivot*m+j]=a;a=bi[k*m+j];bi[k*m+j]=bi[pivot*m+j];bi[pivot*m+j]=a;}
      }
      const pr=ar[k*n+k],pp=ai[k*n+k],den=pr*pr+pp*pp;
      for(let i=k+1;i<n;i++){
        const xr=ar[i*n+k],xi=ai[i*n+k],fr=(xr*pr+xi*pp)/den,fi=(xi*pr-xr*pp)/den;
        ar[i*n+k]=ai[i*n+k]=0;
        if(fr===0&&fi===0)continue;
        for(let j=k+1;j<n;j++){const ur=ar[k*n+j],ui=ai[k*n+j];ar[i*n+j]-=fr*ur-fi*ui;ai[i*n+j]-=fr*ui+fi*ur;}
        for(let j=0;j<m;j++){const ur=br[k*m+j],ui=bi[k*m+j];br[i*m+j]-=fr*ur-fi*ui;bi[i*m+j]-=fr*ui+fi*ur;}
      }
    }
    for(let i=n-1;i>=0;i--)for(let c=0;c<m;c++){
      let sr=br[i*m+c],si=bi[i*m+c];
      for(let j=i+1;j<n;j++){const ur=ar[i*n+j],ui=ai[i*n+j],xr=br[j*m+c],xi=bi[j*m+c];sr-=ur*xr-ui*xi;si-=ur*xi+ui*xr;}
      const a=ar[i*n+i],b=ai[i*n+i],den=a*a+b*b;br[i*m+c]=(sr*a+si*b)/den;bi[i*m+c]=(si*a-sr*b)/den;
    }
    let residual=0;
    for(let i=0;i<n;i++)for(let c=0;c<m;c++){
      let sr=-rhs[i*m+c],si=0;
      for(let j=0;j<n;j++){const a=origR[i*n+j],b=origI[i*n+j],xr=br[j*m+c],xi=bi[j*m+c];sr+=a*xr-b*xi;si+=a*xi+b*xr;}
      residual=Math.max(residual,Math.hypot(sr,si));
    }
    assert(residual<1e-9);
    const sr=new Float64Array(m*m),si=new Float64Array(m*m);
    for(let i=0;i<m;i++)for(let j=0;j<m;j++)for(let k=0;k<n;k++){sr[i*m+j]-=rhs[k*m+i]*br[k*m+j];si[i*m+j]-=rhs[k*m+i]*bi[k*m+j];}
    return {re:sr,im:si,size:m,residual};
  }
  function schur(n,L,delta,g){
    const h=make(n,L,{delta,g}),ps=h.P.map(i=>h.states[i]),rows=[];
    for(const E of [-1,0,1])for(const eta of [.5,1]){
      const z=complexSchur(h,E,eta);
      for(let i=0;i<n;i++)for(let j=0;j<n;j++)if(i!==j)for(let d=1;d<n;d++){
        const source=[],target=[];
        ps.forEach((s,k)=>{if(s.p===j)target.push(k);if(s.p===i&&s.w.includes(i)&&s.w.includes(j)&&Math.abs(s.w.indexOf(i)-s.w.indexOf(j))>=d)source.push(k);});
        if(!source.length)continue;
        let sum=0;for(const a of target)for(const b of source)sum+=z.re[a*z.size+b]**2+z.im[a*z.size+b]**2;
        const frob=Math.sqrt(sum),bound=1/(8*eta)*(1+5*eta/21)**(-d);assert(frob<=bound+1e-10);
        rows.push({E,eta,source:i,receiver:j,min_distance:d,source_columns:source.length,
          all_consistent_receiver_rows:target.length,Sigma_Frobenius_upper:frob,operator_bound:bound,solve_residual:z.residual});
      }
    }
    return {N:n,receiver_modes:L,delta,g,full_dimension:h.size,Q_dimension:h.Q.length,rows};
  }
  function source(n,L,length,delta,g){
    const h=make(n,L,{delta,g}),one=S.oneParticle(n,L,delta,g),fz=evolve(one,unit(one.size,0),4*n*n);
    const emitted=1-fz.re[0]**2-fz.im[0]**2;assert(emitted>1e-8);
    const phi={re:fz.re.slice(1),im:fz.im.slice(1)};for(let j=0;j<L;j++){phi.re[j]/=Math.sqrt(emitted);phi.im[j]/=Math.sqrt(emitted);}
    const w=Array.from({length},(_,i)=>i),q=mask(w),initial={re:new Float64Array(h.size),im:new Float64Array(h.size)};
    for(const o of S.occupations(length-1,L)){
      const [a,b]=fockAmplitude(phi,o);
      for(const p of [0,1]){const k=sourceIndex(h,w,q,o,p);initial.re[k]=a/Math.sqrt(2);initial.im[k]=b/Math.sqrt(2);}
    }
    assert(Math.abs(norm2(initial.re,initial.im)-1)<1e-10);
    const y0={re:initial.re.slice(),im:initial.im.slice()},plus={re:initial.re.slice(),im:initial.im.slice()},minus={re:initial.re.slice(),im:initial.im.slice()};
    for(let i=0;i<h.size;i++)if(h.states[i].p===0){
      plus.re[i]=initial.im[i];plus.im[i]=-initial.re[i];minus.re[i]=-initial.im[i];minus.im[i]=initial.re[i];
    }else y0.re[i]=y0.im[i]=0;
    const rows=[];
    for(const t of [.005,.05,.25]){
      const x=evolve(h,initial,t),y=evolve(h,y0,t),xp=evolve(h,plus,t),xm=evolve(h,minus,t),chi=Array(n).fill(0),fd=Array(n).fill(0);
      let bad=0,prefixBad=0;const keep=Math.ceil(length/2),need=length-keep+1;
      for(let i=0;i<h.size;i++){
        const p=h.states[i].p;chi[p]+=2*(x.re[i]*y.im[i]-x.im[i]*y.re[i]);fd[p]+=(xp.re[i]**2+xp.im[i]**2-xm.re[i]**2-xm.im[i]**2)/2;
        if(h.states[i].q!==mask(h.states[i].w))bad+=x.re[i]**2+x.im[i]**2;
        if(h.states[i].w.slice(0,keep).join(",")!==w.slice(0,keep).join(","))prefixBad+=x.re[i]**2+x.im[i]**2;
      }
      const error=Math.max(...chi.map((x,j)=>Math.abs(x-fd[j])));assert(error<1e-9);
      const near=.25*t-Math.expm1(4.2*t)+4.2*t;if(near>0)assert(chi[1]>=near-1e-10);
      const far={};for(const k of w){const bound=S.tail(k,2.1*t)*S.tail(Math.abs(k-1),2.1*t);far[k]=bound;assert(Math.abs(chi[k])<=bound+1e-10);}
      const pref=S.tail(need,2.7/Math.sqrt(n)*t)**2;assert(prefixBad<=pref+1e-10);assert(bad<=.0625*t*t+1e-10);
      rows.push({t,chi_all_storage_and_receiver_modes:chi,pulse_identity_error:error,near_response_lower:near,
        receiver_response_bounds:far,mismatch_probability:bad,mismatch_Duhamel_upper:.0625*t*t,
        bulk_prefix_failure_probability:prefixBad,bulk_prefix_bound:pref});
    }
    return {N:n,receiver_modes:L,delta,g,recorded_word:w,variance:.25,
      conditional_field_preparation_time:4*n*n,
      finite_small_receiver_macroscopic_probability:S.binomialMacro(n,1-emitted),
      small_receiver_notice:"This exact matrix check uses the stated small L, not L_N. The all-N formation theorem uses L_N=64*N^2.",
      field_state:"Exactly conditioned from the neutral seed via quadratic evolution; no reservoir reset",rows};
  }
  function propagation(n,L,delta,g){
    const h=make(n,L,{delta,g}),rows=[];
    for(const t of [.05,.25])for(let d=1;d<n;d++){
      const source=h.P.filter(a=>{const s=h.states[a];return s.p===0&&s.w.includes(1)&&s.w.indexOf(1)>=d;});
      let square=0;
      for(const a of source){const z=evolve(h,unit(h.size,a),t);for(let j=0;j<h.size;j++)if(h.states[j].p===1)square+=z.re[j]**2+z.im[j]**2;}
      const frob=Math.sqrt(square),bound=S.tail(d,2.1*t);assert(frob<=bound+1e-10);
      rows.push({t,min_distance:d,source_columns:source.length,all_receiver_Frobenius_upper:frob,operator_Dyson_bound:bound});
    }
    return {N:n,receiver_modes:L,delta,g,rows};
  }


  function weylDirect(n,delta){
    const a=weyl(n,delta),L=a.receiver_modes,R=a.box_modes,theta=Math.PI/(R+1),k=Math.acos(-a.one_particle_target_energy/2);
    const re=new Float64Array(L),im=new Float64Array(L),begin=L-R;
    for(let r=1;r<=R;r++){const x=Math.sqrt(2/(R+1))*Math.sin(theta*r);re[begin+r-1]=x*Math.cos(k*r);im[begin+r-1]=x*Math.sin(k*r);}
    let residual=0,mean=0,total=0;
    for(let j=0;j<L;j++){
      const hr=-(j?re[j-1]:0)-(j+1<L?re[j+1]:0),hi=-(j?im[j-1]:0)-(j+1<L?im[j+1]:0);
      residual+=(hr-a.one_particle_target_energy*re[j])**2+(hi-a.one_particle_target_energy*im[j])**2;
      mean+=re[j]*hr+im[j]*hi;total+=re[j]**2+im[j]**2;
    }
    const bias=mean-a.one_particle_target_energy,predicted=a.Q_residual_norm**2/(n-1)-(n-2)*(4*(-a.one_particle_target_energy/2)*Math.sin(theta/2)**2)**2;
    assert(Math.abs(total-1)<1e-10);assert(Math.abs(residual-predicted)<1e-12);
    const qResidual=Math.sqrt((n-1)*residual+(n-1)*(n-2)*bias*bias);
    assert(Math.abs(qResidual-a.Q_residual_norm)<2e-8);
    return {...a,direct_one_particle_norm:total,direct_Q_residual:qResidual,
      direct_vs_formula_error:Math.abs(qResidual-a.Q_residual_norm),direct_contact_amplitude:Math.hypot(re[0],im[0])};
  }

  return {unit,graph,reduction,continuum,finiteReceiver,weyl,weylDirect,schur,source,propagation};
}
if(typeof module!=="undefined"&&module.exports)module.exports=createControls;
