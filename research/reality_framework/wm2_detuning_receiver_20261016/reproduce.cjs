"use strict";
// Run with: node research/reality_framework/wm2_detuning_receiver_20261016/reproduce.cjs
// Add --recompute to rerun every bounded numerical control before refreshing hashes.
const fs=require("node:fs"),path=require("node:path"),sha=require("./sha256.cjs"),S=require("./model.cjs");
const verify=require("./verify.cjs"),C=require("./controls.cjs")(S),runSuite=require("./suite.cjs");
const base="9cbfef8b70f90c875df7b87b59b945c7bc99d02e",dir=__dirname;
const read=name=>fs.readFileSync(path.join(dir,name),"utf8"),write=(name,value)=>fs.writeFileSync(path.join(dir,name),value);
const refresh=process.argv.includes("--recompute");
if(refresh){
  const r=runSuite(S,C,"Node.js "+process.version+"; recomputed from published source");
  write("RESULTS.json",JSON.stringify(r,null,2)+"\n");
}
const r=JSON.parse(read("RESULTS.json")),a=JSON.parse(read("ASSESSMENT.json")),stats=verify(S,r,a);
const names=fs.readdirSync(dir).filter(x=>/\.(md|cjs|json)$/.test(x)&&!["manifest.json","VERIFICATION.json"].includes(x)).sort();
let manifest=JSON.parse(read("manifest.json"));
const entries=Object.fromEntries(names.map(name=>[name,{sha256:sha(read(name)),bytes:Buffer.byteLength(read(name),"utf8")}]));
if(refresh){manifest={...manifest,files:entries};write("manifest.json",JSON.stringify(manifest,null,2)+"\n");}
else S.assert(JSON.stringify(entries)===JSON.stringify(manifest.files),"Content hash mismatch; inspect changes before --recompute");
S.assert(manifest.research_base_commit===base);
const old=path.join(dir,"../wm2_fuel_lock_20261016/manifest.json");
if(fs.existsSync(old))S.assert(sha(fs.readFileSync(old,"utf8"))===manifest.frozen_S2_manifest_sha256,"Frozen S2 manifest changed");
const verification={...stats,packet_files_hashed:names.length,manifest_sha256:sha(read("manifest.json"))};
if(refresh)write("VERIFICATION.json",JSON.stringify(verification,null,2)+"\n");
else S.assert(JSON.stringify(verification)===JSON.stringify(JSON.parse(read("VERIFICATION.json"))),"Saved verification record mismatch");
console.log(JSON.stringify(verification,null,2));
