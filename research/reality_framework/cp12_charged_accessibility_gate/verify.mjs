import {calculate as calculateOnset} from "./exact_onset.mjs";
import {calculate as calculateBounds} from "./analytic_bounds.mjs";
import {readFileSync} from "node:fs";
import {createHash} from "node:crypto";
function checkCore(){

const exact=calculateOnset(), bounds=calculateBounds();
if(exact.minimum_permitted_hops!==22||exact.probability_leading_order!==44||
   exact.minimum_configurations!==960||exact.tournaments!==24||
   exact.orbits.length!==8||!exact.independent_forward_reverse_all_orbits_agree||
   exact.exact_norm_square_PG_H22_psi!=="1974321263356092214149120"||
   exact.leading_coefficient!=="1326415312/848779513441154588027109375")
 throw Error("exact certificate mismatch");
if(!exact.all_960_forward_amplitudes_directly_checked)throw Error("exhaustive direct check missing");
if(exact.orbits.some(x=>x.orbit_size!==120||x.exact_signed_amplitude!==x.independent_reverse_amplitude))
 throw Error("orbit mismatch");
for(const x of bounds.rows){
 if(x.probability_upper===0)throw Error("underflow misreported as zero");
 if(x.d!==4*x.N+Math.ceil(x.N/4))throw Error("distance");
 if(x.probability_upper!==null&&(x.probability_upper<0||x.probability_upper>1))throw Error("probability");
}
function choose(n,k){let p=1;for(let j=1;j<=k;j++)p=p*(n-j+1)/j;return p;}
function fac(n){let p=1;for(let j=2;j<=n;j++)p*=j;return p;}
let toy=0,maxViolation=0;
for(let n=1;n<=16;n++)for(let d=1;d<=n;d++)for(let T of [.01,.05,.1,.2,.4]){
 const p=Math.sin(T)**2;let actual=0;
 for(let m=d;m<=n;m++)actual+=choose(n,m)*p**m*(1-p)**(n-m);
 const upper=Math.min(1,(n*T)**(2*d)/(n*fac(d)**2));
 maxViolation=Math.max(maxViolation,actual-upper);
 if(actual>upper+1e-12)throw Error("unitary-remainder independent-spin check");
 toy++;
}
return {status:"CP12_CHARGED_ACCESS_EXACT_CORE_PASS",exact_integer_arithmetic:"BigInt",
 all_960_minimum_configurations_in_8_orbits_accounted_for:true,
 all_960_forward_amplitudes_directly_checked:true,
 forward_and_reverse_independent_parity_implementations_agree:true,
 single_configuration_cross_language_check:"Python before the container loss and independent JavaScript afterward both returned -49583226880 for identical fields",
 analytic_bound_evaluations:bounds.rows.length,finite_unitary_remainder_sanity_cases:toy,
 maximum_sanity_probability_violation:maxViolation,
 numerical_session_commands_completed:true,
 packaged_python_numerical_reproduction_executed:false,
 full_packet_reproduction_pass_claimed:false,
 external_mathematical_certification:false,
 scope:"Exact core and analytic bound evaluation executed independently in the orchestration JavaScript runtime after container loss. The sanity cases are independent-spin systems, not CP12 formation simulations."};
}
const manifest=JSON.parse(readFileSync(new URL("./MANIFEST.json",import.meta.url),"utf8"));
for(const [name,expected] of Object.entries(manifest.packet_git_blob_sha1)){
 const data=readFileSync(new URL("./"+name,import.meta.url));
 const actual=createHash("sha1").update("blob "+data.length+"\0").update(data).digest("hex");
 if(actual!==expected)throw Error("File hash mismatch: "+name);
}
for(const entry of manifest.parent_dependency_git_blob_sha1){
 const data=readFileSync(new URL("../"+entry.path,import.meta.url));
 const actual=createHash("sha1").update("blob "+data.length+"\0").update(data).digest("hex");
 if(actual!==entry.git_blob_sha1)throw Error("Dependency hash mismatch: "+entry.path);
}
const actual=checkCore();
const exactExpected=JSON.parse(readFileSync(new URL("./EXACT_ONSET.json",import.meta.url),"utf8"));
if(JSON.stringify(calculateOnset())!==JSON.stringify(exactExpected))throw Error("Exact saved-result mismatch");
console.log(JSON.stringify(actual,null,2));
