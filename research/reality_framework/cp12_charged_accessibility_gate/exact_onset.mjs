// Exact BigInt shortest-path certificate; no dependencies.
export function calculate() {

const N=5, es=[]; for(let i=0;i<N;i++)for(let j=i+1;j<N;j++)es.push([i,j]);
const pop=x=>{let c=0;while(x){x&=x-1;c++;}return c;};
const factorial=n=>{let a=1n;for(let j=2;j<=n;j++)a*=BigInt(j);return a;};
const gcd=(a,b)=>{a=a<0n?-a:a;while(b){let c=a%b;a=b;b=c;}return a;};
const rational=(a,b)=>{const g=gcd(a,b);return [a/g,b/g];};
const rs=(a,b)=>{const [x,y]=rational(a,b);return y===1n?String(x):x+"/"+y;};
function state(E){let bits=0,v=[0,0,0,0,0];for(let e=0;e<10;e++){const [i,j]=es[e],a=E[2*e],b=E[2*e+1];v[i]+=a;v[j]+=b;if(a+b===0)bits|=1<<(5+e);else if(a+b!==1)throw Error("edge Gauss");}for(let i=0;i<5;i++){if(v[i]!==0&&v[i]!==1)throw Error("vertex Gauss");bits|=v[i]<<i;}if(pop(bits)!==10)throw Error("fermion number");return bits;}
function amp(E,reverse=false){let finalbits=state(E),twice=0;E.forEach((v,j)=>{if(Math.abs(v)===2)twice|=1<<j;});
let v=new Map([[reverse? ((1048575)+","+twice+","+finalbits):("0,0,"+(32736)),1n]]),counts=[1];
const D=E.reduce((s,x)=>s+Math.abs(x),0);
for(let step=0;step<D;step++){let w=new Map();for(let [key,a] of v){let [once,two,bits]=key.split(",").map(Number);
for(let j=0;j<20;j++){const bit=1<<j;let no,nt;
if(reverse){if(!(once&bit))continue;if(two&bit){no=once;nt=two^bit;}else{no=once^bit;nt=two;}}
else if(once&bit){if(!(twice&bit)||(two&bit))continue;no=once;nt=two|bit;}else {no=once|bit;nt=two;}
const e=j>>1,leg=j&1,vertex=es[e][leg];let dst=E[j]>0?vertex:5+e,source=E[j]>0?5+e:vertex;
if(reverse)[dst,source]=[source,dst];if(!((bits>>source)&1)||((bits>>dst)&1))continue;
let sign;
if(reverse){let low=Math.min(dst,source),hi=Math.max(dst,source),mask=((1<<hi)-1)^((1<<(low+1))-1);sign=pop(bits&mask)&1?-1n:1n;}
else{const after=bits^(1<<source);sign=(pop(bits&((1<<source)-1))+pop(after&((1<<dst)-1)))&1?-1n:1n;}
const nb=bits^(1<<source)^(1<<dst), k=no+","+nt+","+nb;
w.set(k,(w.get(k)||0n)+a*sign*(leg?-1n:1n));
}}v=new Map([...w].filter(x=>x[1]!==0n));counts.push(v.size);}
const key=reverse?"0,0,32736":("1048575,"+twice+","+finalbits);
if(v.size!==1||!v.has(key))throw Error("unexpected terminal");return {amplitude:v.get(key),counts};}
function perms(a){if(a.length===0)return [[]];return a.flatMap((v,j)=>perms(a.filter((_,k)=>k!==j)).map(p=>[v,...p]));}
const maps=perms([0,1,2,3,4]).map(p=>es.flatMap(([i,j])=>{let e=es.findIndex(([a,b])=>a===Math.min(p[i],p[j])&&b===Math.max(p[i],p[j]));return [2*e+ +(p[i]>p[j]),2*e+ +(p[j]>p[i])];}));
let targets=[], tournaments=0;
for(let o=0;o<1024;o++){let tails=es.map(([i,j],e)=>(o>>e)&1?i:j);if([0,1,2,3,4].some(i=>tails.filter(x=>x===i).length!==2))continue;tournaments++;
for(let a=0;a<10;a++)for(let b=a+1;b<10;b++){if(tails[a]===tails[b])continue;let E=es.flatMap(([i,j],e)=>[tails[e]===i?1+ +(e===a||e===b):-1,tails[e]===j?1+ +(e===a||e===b):-1]);state(E);targets.push(E);}}
let groups=new Map();
for(let E of targets){let orb=maps.map(m=>{let v=Array(20);for(let j=0;j<20;j++)v[m[j]]=E[j];return v.join(",");});orb.sort();const key=orb[0];if(!groups.has(key))groups.set(key,[]);groups.get(key).push(E);}
let rows=[],total=0n;
for(let [key,grp] of [...groups].sort()){let E=key.split(",").map(Number),f=amp(E),r=amp(E,true);if(f.amplitude!==r.amplitude)throw Error("reverse mismatch");
let amount=0n;
for(const fields of grp){
 const direct=amp(fields).amplitude;
 if(direct*direct!==f.amplitude*f.amplitude)throw Error("direct 960-state sign check");
 amount+=direct*direct;
}
total+=amount;
rows.push({electric_fields:E,orbit_size:grp.length,exact_signed_amplitude:String(f.amplitude),independent_reverse_amplitude:String(r.amplitude),maximum_forward_layer:Math.max(...f.counts),maximum_reverse_layer:Math.max(...r.counts)});
}
let den=factorial(22)**2n;const coeff=rational(total,den);
let cnum=20n,cden=1n;for(let j=1;j<22;j++){cnum*=BigInt(80+16*Math.min(5,j)+j*j)**2n;cden*=64n;}const overlap=rational(total*cden,cnum);
const pythonTargetFields=[2,-1,1,-1,-1,1,-1,1,2,-1,1,-1,-1,1,1,-1,1,-1,1,-1];
const crossLanguage=amp(pythonTargetFields);
if(crossLanguage.amplitude!==-49583226880n)throw Error("Python/JavaScript single-target mismatch");
return {python_cross_language_target_fields:pythonTargetFields,python_cross_language_exact_amplitude:String(crossLanguage.amplitude),N:5,Q:10,charged_target:"all ten edges circulate, n_v>=2",minimum_permitted_hops:22,tournaments,minimum_configurations:targets.length,orbits:rows,exact_norm_square_PG_H22_psi:String(total),probability_leading_order:44,leading_coefficient:rs(total,den),leading_coefficient_float:Number(coeff[0])/Number(coeff[1]),actual_initial_cyclic_target_operator_norm_squared_lower:rs(total*cden,cnum),actual_initial_cyclic_target_operator_norm_squared_lower_float:Number(overlap[0])/Number(overlap[1]),independent_forward_reverse_all_orbits_agree:true,all_960_forward_amplitudes_directly_checked:true};

}
if (process.argv[1] && import.meta.url === new URL("file://" + process.argv[1]).href) console.log(JSON.stringify(calculate(), null, 2));
