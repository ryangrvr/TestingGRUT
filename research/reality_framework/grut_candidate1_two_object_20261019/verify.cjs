"use strict";
const fs=require("node:fs"),path=require("node:path");
const protocol=require("./protocol.json"),reference=require("./source_observables.json");
const output=path.join(__dirname,"results");
function linear(a,x0,h,x){
  const f=(x-x0)/h;
  if(f<0||f>a.length-1)return 0;
  const j=Math.min(a.length-2,Math.max(0,Math.floor(f))),t=f-j;
  return (1-t)*a[j]+t*a[j+1];
}
function difference(a,b){
  const sa=a.analysis.final,sb=b.analysis.final;
  let maxRms=0,maxMiddle=0,maxRelativeNorm=0;
  for(let i=0;i<Math.min(a.trace.length,b.trace.length);i++){
    maxRms=Math.max(maxRms,Math.abs(a.trace[i].z_rms-b.trace[i].z_rms));
    maxMiddle=Math.max(maxMiddle,Math.abs(a.trace[i].middle_norm_fraction-b.trace[i].middle_norm_fraction));
    maxRelativeNorm=Math.max(maxRelativeNorm,Math.abs(a.trace[i].normF2-b.trace[i].normF2)/Math.max(b.trace[i].normF2,1e-100));
  }
  let profileL1=0;
  for(let q=0;q<b.profiles.length;q++){
    const pa=a.profiles.find(x=>Math.abs(x.time-b.profiles[q].time)<1e-8);
    const pb=b.profiles[q];if(!pa)continue;
    const h=b.grid.h,z0=-b.configuration.Z+h/2;
    let d=0,n=0;
    for(let j=0;j<pb.lineDensity.length;j++){
      const z=z0+j*h;
      if(Math.abs(z)>Math.min(a.configuration.Z,b.configuration.Z)-Math.max(a.grid.h,b.grid.h))continue;
      const x=linear(pa.lineDensity,-a.configuration.Z+a.grid.h/2,a.grid.h,z);
      d+=h*Math.abs(x-pb.lineDensity[j]);n+=h*Math.abs(pb.lineDensity[j]);
    }
    profileL1=Math.max(profileL1,d/Math.max(n,1e-100));
  }
  return {max_z_rms_difference:maxRms,max_middle_fraction_difference:maxMiddle,
    max_relative_norm_difference:maxRelativeNorm,max_normalized_line_profile_L1:profileL1,
    final_z_rms_difference:sa.z_rms-sb.z_rms,final_middle_fraction_difference:sa.middle_norm_fraction-sb.middle_norm_fraction};
}
function components(run,fraction){
  const p=run.profiles[run.profiles.length-1],arr=p.lineDensity,h=run.grid.h;
  const maximum=Math.max(...arr),groups=[];
  let start=-1;
  for(let j=0;j<=arr.length;j++){
    const inside=j<arr.length&&arr[j]>=fraction*maximum;
    if(inside&&start<0)start=j;
    if(!inside&&start>=0){
      let norm=0,q=0,zm=0,peak=start;
      for(let k=start;k<j;k++){
        const z=-run.configuration.Z+(k+.5)*h;
        norm+=arr[k]*h;q+=p.lineCharge[k]*h;zm+=z*arr[k]*h;
        if(arr[k]>arr[peak])peak=k;
      }
      groups.push({z_min:-run.configuration.Z+(start+.5)*h,z_max:-run.configuration.Z+(j-.5)*h,
        peak_z:-run.configuration.Z+(peak+.5)*h,norm_fraction:norm/run.analysis.final.normF2,
        Qbar:q,z_centroid:zm/norm});
      start=-1;
    }
  }
  return groups;
}
function summarize(results,controls){
  const map=Object.fromEntries(results.map(x=>[x.id,x])),checks=[],add=(name,pass,detail)=>checks.push({name,pass,detail});
  add("All prescribed runs present",results.length===protocol.runs.length,{expected:protocol.runs.length,actual:results.length});
  add("Vacuum and reversal controls",controls.pass,controls);
  let maxE=0,maxQ=0;
  for(const r of results){const a=r.analysis;maxE=Math.max(maxE,a.maximum_relative_energy_drift);maxQ=Math.max(maxQ,a.maximum_charge_error_over_initial_absolute_charge);
    add(r.id+" finite, conserved within declared tolerances",a.finite&&a.maximum_relative_energy_drift<protocol.acceptance.maximum_relative_energy_drift&&a.maximum_charge_error_over_initial_absolute_charge<protocol.acceptance.maximum_absolute_charge_error_over_initial_absolute_charge,
      {energy_drift:a.maximum_relative_energy_drift,charge_error:a.maximum_charge_error_over_initial_absolute_charge});
  }
  const errors={};
  for(const level of ["coarse","medium","fine"]){
    const r=map[level+"_single_static"];if(!r)continue;
    errors[level]={E_relative:Math.abs(r.analysis.initial.Ebar/reference.Ebar-1),
      Q_relative:Math.abs(r.analysis.initial.Qbar/reference.Qbar-1)};
  }
  if(errors.coarse&&errors.fine)add("Initial one-object continuum E and Q improve with refinement",
    errors.fine.E_relative<errors.coarse.E_relative&&errors.fine.Q_relative<errors.coarse.Q_relative,errors);
  const comparisons={};
  for(const name of protocol.base.cases){
    if(map["coarse_"+name]&&map["medium_"+name])comparisons["coarse_to_medium_"+name]=difference(map["coarse_"+name],map["medium_"+name]);
    if(map["medium_"+name]&&map["fine_"+name])comparisons["medium_to_fine_"+name]=difference(map["medium_"+name],map["fine_"+name]);
  }
  for(const id of ["medium_half_dt_pair_phase_0","medium_larger_box_pair_phase_0","medium_larger_box_pair_phase_pi"]){
    const name=id.endsWith("phase_pi")?"pair_phase_pi":"pair_phase_0";
    if(map[id]&&map["medium_"+name])comparisons[id]=difference(map["medium_"+name],map[id]);
  }
  const mirror={};
  for(const level of ["coarse","medium","fine"]){
    const a=map[level+"_pair_phase_pi_minus"],b=map[level+"_pair_phase_pi_plus"];
    if(!a||!b)continue;
    let err=0;
    for(let q=0;q<a.profiles.length;q++){
      const x=a.profiles[q].lineDensity,y=b.profiles[q].lineDensity;let d=0,n=0;
      for(let j=0;j<x.length;j++){d+=Math.abs(x[j]-y[y.length-1-j]);n+=Math.abs(x[j]);}
      err=Math.max(err,d/n);
    }
    mirror[level]=err;add(level+" phase-neighborhood reflection",err<1e-9,{normalized_L1:err});
  }
  const outcomes=results.map(r=>({id:r.id,grid:r.grid,initial:r.analysis.initial,final:r.analysis.final,
    energy_drift:r.analysis.maximum_relative_energy_drift,charge_error:r.analysis.maximum_charge_error_over_initial_absolute_charge,
    final_axial_components:{threshold_2_percent:components(r,.02),threshold_5_percent:components(r,.05),threshold_10_percent:components(r,.1)}}));
  return {source_freeze:"52f552802759392a2353af068f60efea2c5c95b2",runtime:{node:process.version,platform:process.platform,arch:process.arch,workflow_commit:process.env.GITHUB_SHA||null},
    all_consistency_checks_pass:checks.every(x=>x.pass),checks,max_energy_drift:maxE,max_charge_error:maxQ,
    one_object_initial_errors:errors,refinement_comparisons:comparisons,phase_reflection_errors:mirror,outcomes,
    interpretation:"Numerical consistency checks do not certify particle identity, an open stable phase, infinite lifetime, formation, quantum behavior or new physics. Threshold components are observable diagnostics."};
}
if(require.main===module){
  const results=protocol.runs.filter(x=>fs.existsSync(path.join(output,x.id+".json"))).map(x=>JSON.parse(fs.readFileSync(path.join(output,x.id+".json"),"utf8")));
  const controls=JSON.parse(fs.readFileSync(path.join(output,"controls.json"),"utf8"));
  const summary=summarize(results,controls);fs.writeFileSync(path.join(output,"summary.json"),JSON.stringify(summary,null,2)+"\n");
  console.log(JSON.stringify({runs:results.length,all_consistency_checks_pass:summary.all_consistency_checks_pass,max_energy_drift:summary.max_energy_drift,max_charge_error:summary.max_charge_error}));
  if(!summary.all_consistency_checks_pass)process.exitCode=1;
}
module.exports={difference,components,summarize};
