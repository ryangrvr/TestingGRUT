// Floating-point evaluation of the analytic probability theorem.
// Underflow is recorded as null with a logarithmic bound, never as zero.
export function calculate() {

function logFactorial(n){let x=0;for(let k=2;k<=n;k++)x+=Math.log(k);return x;}
function bound(N,Q,L,k,T,t=1){
 const rank=Math.min(N,Q),d=2*L+k,b=Math.abs(t)*Math.sqrt(2*Q*rank);
 if(T===0||t===0)return {N,Q,L,k,T,d,hop_norm_bound:b,log_probability_upper:null,probability_upper:0};
 const logp=-Math.log(rank)+2*(d*Math.log(b*Math.abs(T))-logFactorial(d));
 return {N,Q,L,k,T,d,hop_norm_bound:b,log_probability_upper:logp,probability_upper:logp < -700 ? null : Math.min(1,Math.exp(logp)),underflow_note:logp < -700 ? "Use the logarithmic bound; probability is not zero" : null,log10_probability_upper:Math.min(0,logp/Math.LN10)};
}
const rows=[];
for(const N of [5,16,64,256,1024])for(const T of [.25,.5,.75,1.])rows.push(bound(N,2*N,2*N,Math.ceil(N/4),T));
return {frozen_parameters:{t:1,m:1,epsilon_E:.25,Gamma:.1,kappa:.5},
 target:"Every active edge lies on a directed electric circuit; n_v >= ceil(N/4); Q=2N",
 amplitude_bound:"sqrt(2 Q)|t| b^(d-1) |T|^d/d!, b=|t|sqrt(2 Q min(N,Q)), d=2L+k",
 probability_bound:"min(1, [1/min(N,Q)] (b |T|)^ (2d)/(d!)^2)",
 extensive_threshold_limit:4.25/(2*Math.E),
 threshold_scope:"Exponential suppression for each fixed T < 4.25/(2e); a configuration-distance bound, not a spatial light cone or long-time no-go",
 rows};

}
if (process.argv[1] && import.meta.url === new URL("file://" + process.argv[1]).href) console.log(JSON.stringify(calculate(), null, 2));
