"use strict";
function figure(data){
 const colors={pair_phase_0:"#2563eb",pair_phase_halfpi:"#d97706",pair_phase_pi_minus:"#16a34a",pair_phase_pi:"#111827",pair_phase_pi_plus:"#9333ea"};
 const names={pair_phase_0:"0",pair_phase_halfpi:"pi/2",pair_phase_pi_minus:"pi - 0.2",pair_phase_pi:"pi",pair_phase_pi_plus:"pi + 0.2"};
 const keys=Object.keys(colors),W=1100,H=750,L=86,R=1045,plotW=R-L;
 const esc=s=>String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");
 const txt=(x,y,s,size=14,extra="")=>'<text x="'+x+'" y="'+y+'" font-size="'+size+'" '+extra+'>'+esc(s)+'</text>';
 const line=(x1,y1,x2,y2,color="#e5e7eb",extra="")=>'<line x1="'+x1+'" y1="'+y1+'" x2="'+x2+'" y2="'+y2+'" stroke="'+color+'" '+extra+'/>';
 let out=['<svg xmlns="http://www.w3.org/2000/svg" width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'">',
 '<title>Candidate 1 axisymmetric two-object results</title><desc>Finest-grid axial field-density profiles at time 80 and left-half charge fraction over time for five prepared relative phases.</desc>',
 '<rect width="1100" height="750" fill="white"/>','<g font-family="Arial, sans-serif" fill="#111827">',txt(42,29,"Candidate 1: one action, phase-dependent encounters",23)];
 keys.forEach((k,i)=>{const x=90+i*198;out.push(line(x,55,x+24,55,colors[k],'stroke-width="3"'),txt(x+33,60,"phase "+names[k],15));});
 const profiles=keys.map(k=>({key:k,run:data[k],p:data[k].profiles.at(-1)}));
 let ymax=0;
 for(const {run,p} of profiles)for(let j=0;j<p.lineDensity.length;j++){
  const z=-run.configuration.Z+(j+.5)*run.grid.h;if(Math.abs(z)<=24)ymax=Math.max(ymax,p.lineDensity[j]/run.analysis.final.normF2);
 }
 ymax*=1.10;
 const ya=93,ha=235,yb=437,hb=225;
 const ax=z=>L+(z+24)/48*plotW,ay=y=>ya+ha-y/ymax*ha;
 const bx=t=>L+t/80*plotW,by=q=>yb+hb-(q-.25)/.5*hb;
 out.push(txt(L,84,"A. Final axial F^2 density, normalized by total F^2",16));
 for(let i=0;i<=4;i++){
  const v=ymax*i/4,y=ay(v);out.push(line(L,y,R,y),txt(L-10,y+5,v.toFixed(3),12,'text-anchor="end"'));
 }
 for(const z of [-24,-12,0,12,24])out.push(line(ax(z),ya,ax(z),ya+ha),txt(ax(z),ya+ha+20,z,13,'text-anchor="middle"'));
 out.push(line(L,ya,L,ya+ha,"#6b7280"),line(L,ya+ha,R,ya+ha,"#6b7280"),txt((L+R)/2,ya+ha+42,"axial coordinate z (g v units)",14,'text-anchor="middle"'));
 for(const {key,run,p} of profiles){
  const points=[];
  for(let j=0;j<p.lineDensity.length;j++){
   const z=-run.configuration.Z+(j+.5)*run.grid.h;if(Math.abs(z)>24)continue;
   points.push(ax(z).toFixed(3)+","+ay(p.lineDensity[j]/run.analysis.final.normF2).toFixed(3));
  }
  out.push('<polyline points="'+points.join(" ")+'" fill="none" stroke="'+colors[key]+'" stroke-width="2.1" '+(key==="pair_phase_pi"?'stroke-dasharray="6 4"':"")+'/>');
 }
 out.push(txt(L,424,"B. Charge in the left half divided by total charge",16));
 for(const q of [.25,.375,.5,.625,.75]){
  const y=by(q);out.push(line(L,y,R,y),txt(L-10,y+5,q.toFixed(3),12,'text-anchor="end"'));
 }
 for(const t of [0,20,40,60,80])out.push(line(bx(t),yb,bx(t),yb+hb),txt(bx(t),yb+hb+20,t,13,'text-anchor="middle"'));
 out.push(line(L,yb,L,yb+hb,"#6b7280"),line(L,yb+hb,R,yb+hb,"#6b7280"),txt((L+R)/2,yb+hb+42,"dimensionless time tau = g v t",14,'text-anchor="middle"'));
 for(const key of keys){
  const run=data[key],pts=run.trace.map(r=>bx(r.time).toFixed(3)+","+by(r.left_Qbar/r.Qbar).toFixed(3)).join(" ");
  out.push('<polyline points="'+pts+'" fill="none" stroke="'+colors[key]+'" stroke-width="2.1" '+(key==="pair_phase_pi"?'stroke-dasharray="6 4"':"")+'/>');
 }
 out.push(txt(42,735,"Prepared equal-charge objects; a=1, b=0.1; speed 0.2; separation 24; h=0.25. Finite-time numerical results.",13),'</g></svg>');
 return out.join("\n");
}
if(require.main===module){
 const fs=require("node:fs"),path=require("node:path"),data={};
 for(const key of ["pair_phase_0","pair_phase_halfpi","pair_phase_pi_minus","pair_phase_pi","pair_phase_pi_plus"])
  data[key]=JSON.parse(fs.readFileSync(path.join(__dirname,"results","fine_"+key+".json"),"utf8"));
 fs.writeFileSync(path.join(__dirname,"RESULTS.svg"),figure(data));
}
module.exports={figure};
