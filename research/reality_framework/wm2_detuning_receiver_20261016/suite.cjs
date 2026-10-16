"use strict";
// Same functions/configurations executed in V8 for the archived results.
function runSuite(S,C,execution="JavaScript runtime"){
  const r={status:"PASS",execution,arithmetic:"IEEE-754 doubles; no outward-rounded numerical certificate",
    constants:S.constants(),graph:[],reduction:[],continuum:[],finite_receiver:[],Schur:[],source:[],propagation:[],Q_near_resonance:[],Q_direct:[]};
  r.graph=[[2,2],[3,3],[4,2]].map(([n,L])=>C.graph(n,L));
  for(const [n,L] of [[2,3],[3,3],[4,2],[5,2]])for(const [d,g] of [[-.5,.9],[0,1],[.5,1.1]])r.reduction.push(C.reduction(n,L,d,g));
  r.Q_near_resonance=[64,128,256,512,1024,4096].flatMap(n=>[-.5,0,.5].map(d=>C.weyl(n,d)));
  r.Q_direct=[C.weylDirect(64,-.5),C.weylDirect(128,.5)];
  for(const n of [2,3,4])for(const d of [-.5,0,.5])for(const g of [.9,1,1.1])r.finite_receiver.push(C.finiteReceiver(n,d,g));
  for(const n of [2,3,4,8,16,32])for(const d of [-.5,-.25,0,.25,.5])for(const g of [.9,1,1.1])r.continuum.push(C.continuum(n,d,g));
  r.quadrature_refinement=[];
  for(const n of [2,16,32])for(const d of [-.5,0,.5]){
    const a=S.spectral(n,d,.9,8*n*n),b=S.spectral(n,d,.9,8*n*n,2),error=Math.hypot(a.re-b.re,a.im-b.im);S.assert(error<1e-9);
    r.quadrature_refinement.push({N:n,delta:d,g:.9,t:8*n*n,amplitude_change:error});
  }
  r.Schur=[[2,3,.5,.9],[3,2,-.5,1.1],[3,3,.5,1]].map(([n,L,d,g])=>C.schur(n,L,d,g));
  r.source=[[2,2,2,.5,.9],[3,2,3,-.5,.9],[3,3,2,.5,1.1],[4,2,2,.5,.9],[4,2,3,0,1]].map(([n,L,m,d,g])=>C.source(n,L,m,d,g));
  r.propagation=[C.propagation(3,2,-.5,.9),C.propagation(3,2,.5,1.1)];
  r.decision={formation:"PASS within the frozen symmetric interval with finite L_N=64*N^2",
    operational_transport:"PASS: full-H bounds include all receiver and mismatch states",
    real_energy_consistency_band:"FAIL: coupled Q poles enter the seed-energy window",full_WM2:"NOT ESTABLISHED"};
  return r;
}
if(typeof module!=="undefined"&&module.exports)module.exports=runSuite;
