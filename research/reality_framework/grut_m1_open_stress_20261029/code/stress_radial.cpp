// GRUT-M1 spherical open-boundary control. Reduced scalar Hamiltonian only.
// Derived from the reviewed daughter-core finite-volume engine; no quantum decay term.
// Symplectic drift/kick; a remote velocity absorber has an explicit receiver.
// r=gv*r_phys, t=gv*t_phys, X=v*(x+i y)/sqrt(2), phi=v*s.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using std::vector;
constexpr double PI=3.14159265358979323846;
struct State {vector<double>x,y,s,px,py,ps;};
struct Geometry {vector<double>r,vol,face,damp; double h,L,dt; int n;};
struct Budget {double E=0,Q=0;std::array<double,3>ec{},qc{};double xmax=0,smin=1;};
struct Stress {double trace_center=0,trace_face=0,virial=0,surface=0;};

Stress stress(const State& u,const Geometry& g,double radius){
 // Canonical scalar trace: 3K-G-3V, physical units divided by g^2.
 // Centered derivatives and the original face-energy convention are both
 // retained; neither turns a finite-radius stress into a full-sector theorem.
 Stress out;
 for(int i=0;i<g.n;i++){
  if(g.r[i]>=radius)break;
  int lo=std::max(0,i-1),hi=std::min(g.n-1,i+1);
  double denom=i==0?2*g.h:(hi-lo)*g.h;
  // At r=0, even regularity gives u[-1]=u[0], hence central derivative
  // (u[1]-u[0])/(2h).
  double dx=(u.x[hi]-u.x[lo])/denom,dy=(u.y[hi]-u.y[lo])/denom,ds=(u.s[hi]-u.s[lo])/denom;
  double n=u.x[i]*u.x[i]+u.y[i]*u.y[i];
  double K=.5*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i]);
  double V=.75*(u.s[i]*u.s[i]-1)*(u.s[i]*u.s[i]-1)-.2*std::pow(u.s[i]-1,3)+.5*u.s[i]*u.s[i]*n+.025*n*n;
  double G=.5*(dx*dx+dy*dy+ds*ds);
  double xn=i+1<g.n?u.x[i+1]:0,yn=i+1<g.n?u.y[i+1]:0,sn=i+1<g.n?u.s[i+1]:1;
  double Gface=.5*g.face[i]*(std::pow(xn-u.x[i],2)+std::pow(yn-u.y[i],2)+std::pow(sn-u.s[i],2));
  out.trace_center+=g.vol[i]*(3*K-G-3*V);
  out.trace_face+=g.vol[i]*(3*K-3*V)-Gface;
  out.virial-=g.vol[i]*g.r[i]*(u.px[i]*dx+u.py[i]*dy+u.ps[i]*ds);
 }
 int i=std::clamp((int)std::lround(radius/g.h)-1,0,g.n-2);
 double x=(u.x[i]+u.x[i+1])/2,y=(u.y[i]+u.y[i+1])/2,s=(u.s[i]+u.s[i+1])/2;
 double px=(u.px[i]+u.px[i+1])/2,py=(u.py[i]+u.py[i+1])/2,ps=(u.ps[i]+u.ps[i+1])/2;
 double dx=(u.x[i+1]-u.x[i])/g.h,dy=(u.y[i+1]-u.y[i])/g.h,ds=(u.s[i+1]-u.s[i])/g.h;
 double n=x*x+y*y,K=.5*(px*px+py*py+ps*ps),G=.5*(dx*dx+dy*dy+ds*ds);
 double V=.75*(s*s-1)*(s*s-1)-.2*std::pow(s-1,3)+.5*s*s*n+.025*n*n;
 double actual=(i+1)*g.h;
 out.surface=4*PI*actual*actual*actual*(K+G-V);
 out.trace_center/=.01;out.trace_face/=.01;out.virial/=.01;out.surface/=.01;
 return out;
}

Budget budget(const State& u,const Geometry& g){
 Budget out;const std::array<double,3> cuts={12.,20.,30.};
 for(int i=0;i<g.n;i++){
  double n=u.x[i]*u.x[i]+u.y[i]*u.y[i];
  double e=g.vol[i]*(.5*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i])
       +.75*std::pow(u.s[i]*u.s[i]-1,2)-.2*std::pow(u.s[i]-1,3)+.5*u.s[i]*u.s[i]*n+.025*n*n);
  double q=g.vol[i]*(u.y[i]*u.px[i]-u.x[i]*u.py[i]);
  double xn=i+1<g.n?u.x[i+1]:0,yn=i+1<g.n?u.y[i+1]:0,sn=i+1<g.n?u.s[i+1]:1;
  double eg=.5*g.face[i]*(std::pow(xn-u.x[i],2)+std::pow(yn-u.y[i],2)+std::pow(sn-u.s[i],2));
  out.E+=e+eg;out.Q+=q;out.xmax=std::max(out.xmax,std::sqrt(n));out.smin=std::min(out.smin,u.s[i]);
  for(int j=0;j<3;j++)if(g.r[i]<cuts[j]){out.ec[j]+=e+eg;out.qc[j]+=q;}
 }
 return out;
}

void force(const State& u,const Geometry& g,vector<double>& fx,vector<double>& fy,vector<double>& fs){
 for(int i=0;i<g.n;i++){
  double xn=i+1<g.n?u.x[i+1]:0,yn=i+1<g.n?u.y[i+1]:0,sn=i+1<g.n?u.s[i+1]:1;
  double lx=g.face[i]*(xn-u.x[i]),ly=g.face[i]*(yn-u.y[i]),ls=g.face[i]*(sn-u.s[i]);
  if(i){lx-=g.face[i-1]*(u.x[i]-u.x[i-1]);ly-=g.face[i-1]*(u.y[i]-u.y[i-1]);ls-=g.face[i-1]*(u.s[i]-u.s[i-1]);}
  double n=u.x[i]*u.x[i]+u.y[i]*u.y[i],mass=u.s[i]*u.s[i]+.1*n;
  fx[i]=lx/g.vol[i]-mass*u.x[i];fy[i]=ly/g.vol[i]-mass*u.y[i];
  fs[i]=ls/g.vol[i]-(3*(u.s[i]*u.s[i]-1)+n)*u.s[i]+.6*std::pow(u.s[i]-1,2);
 }
}

void damping(State& u,const Geometry& g,double& energy,double& charge){
 for(int i=0;i<g.n;i++)if(g.damp[i]<1){
  double f=g.damp[i],q=g.vol[i]*(u.y[i]*u.px[i]-u.x[i]*u.py[i]);
  energy+=.5*g.vol[i]*(1-f*f)*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i]);
  charge+=(1-f)*q;u.px[i]*=f;u.py[i]*=f;u.ps[i]*=f;
 }
}

int main(int argc,char**argv){try{
 if(argc!=10&&argc!=11)throw std::runtime_error("usage radial INITIAL PREFIX h dt L tmax probe absorber_start eta [stress_sample]");
 std::string input=argv[1],prefix=argv[2];Geometry g;
 g.h=std::stod(argv[3]);g.dt=std::stod(argv[4]);g.L=std::stod(argv[5]);
 double end=std::stod(argv[6]),probe=std::stod(argv[7]),abstart=std::stod(argv[8]),eta=std::stod(argv[9]);
 g.n=std::lround(g.L/g.h);g.r.resize(g.n);g.vol.resize(g.n);g.face.resize(g.n);g.damp.resize(g.n);
 State u;for(auto*p:{&u.x,&u.y,&u.s,&u.px,&u.py,&u.ps})p->resize(g.n);
 std::ifstream init(input);std::string line;std::getline(init,line);
 for(int i=0;i<g.n;i++){
  if(!std::getline(init,line))throw std::runtime_error("initial grid too short");
  std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double r;
  row>>r>>u.x[i]>>u.y[i]>>u.s[i]>>u.px[i]>>u.py[i]>>u.ps[i];
  g.r[i]=(i+.5)*g.h;if(std::abs(r-g.r[i])>1e-8)throw std::runtime_error("initial grid mismatch");
  double lo=i*g.h,hi=(i+1)*g.h;
  g.vol[i]=4*PI/3*(hi*hi*hi-lo*lo*lo);g.face[i]=4*PI*hi*hi/g.h*(i+1==g.n?2:1);
  double z=std::max(0.,(g.r[i]-abstart)/(g.L-abstart));
  g.damp[i]=std::exp(-.5*g.dt*eta*z*z*z);
 }
 vector<double> fx(g.n),fy(g.n),fs(g.n);double ea=0,qa=0,energyFlux=0,chargeFlux=0;
 std::ofstream history(prefix+"_history.csv"),trace(prefix+"_trace.csv"),pressure(prefix+"_stress.csv");
 history<<std::setprecision(17);trace<<std::setprecision(17);pressure<<std::setprecision(17);
 history<<"tau,E,Q,E_abs,Q_abs,Eflux_probe,Qflux_probe,Fmax,Smin,E12,Q12,E20,Q20,E30,Q30\n";
 trace<<"tau,x,y,s_minus1,dx,dy,ds,px,py,ps\n";
 pressure<<"tau,trace20_center,trace20_face,virial20,surface20,trace60_center,trace60_face,virial60,surface60\n";
 int ip=std::lround(probe/g.h-.5);ip=std::clamp(ip,1,g.n-2);
 int steps=std::lround(end/g.dt),histstride=std::max(1,(int)std::lround(10/g.dt));
 int tracestride=std::max(1,(int)std::lround(.25/g.dt));
 double stresssample=argc==11?std::stod(argv[10]):1.;
 int stressstride=std::max(1,(int)std::lround(stresssample/g.dt));
 auto probevalues=[&](){return std::array<double,3>{(u.x[ip+1]-u.x[ip-1])/(2*g.h),(u.y[ip+1]-u.y[ip-1])/(2*g.h),(u.s[ip+1]-u.s[ip-1])/(2*g.h)};};
 double prevEF=0,prevQF=0;
 for(int step=0;step<=steps;step++){
  double t=step*g.dt;auto d=probevalues();double area=4*PI*g.r[ip]*g.r[ip];
  double ef=-area*(u.px[ip]*d[0]+u.py[ip]*d[1]+u.ps[ip]*d[2]);
  double qf=area*(u.x[ip]*d[1]-u.y[ip]*d[0]);
  if(step){energyFlux+=.5*g.dt*(ef+prevEF);chargeFlux+=.5*g.dt*(qf+prevQF);}prevEF=ef;prevQF=qf;
  if(step%histstride==0||step==steps){
   auto b=budget(u,g);if(!std::isfinite(b.E)||b.xmax>100)throw std::runtime_error("unresolved runaway/nonfinite state");
   history<<t<<','<<b.E<<','<<b.Q<<','<<ea<<','<<qa<<','<<energyFlux<<','<<chargeFlux<<','<<b.xmax<<','<<b.smin;
   for(int j=0;j<3;j++){history<<','<<b.ec[j]<<','<<b.qc[j];}
   history<<'\n';
  }
  if(step%tracestride==0||step==steps)trace<<t<<','<<u.x[ip]<<','<<u.y[ip]<<','<<u.s[ip]-1<<','<<d[0]<<','<<d[1]<<','<<d[2]<<','<<u.px[ip]<<','<<u.py[ip]<<','<<u.ps[ip]<<'\n';
  if(step%stressstride==0||step==steps){
   auto a=stress(u,g,20.);auto b=stress(u,g,60.);
   pressure<<t<<','<<a.trace_center<<','<<a.trace_face<<','<<a.virial<<','<<a.surface<<','
           <<b.trace_center<<','<<b.trace_face<<','<<b.virial<<','<<b.surface<<'\n';
  }
  if(step==steps)break;
  damping(u,g,ea,qa);force(u,g,fx,fy,fs);
  for(int i=0;i<g.n;i++){u.px[i]+=.5*g.dt*fx[i];u.py[i]+=.5*g.dt*fy[i];u.ps[i]+=.5*g.dt*fs[i];}
  for(int i=0;i<g.n;i++){u.x[i]+=g.dt*u.px[i];u.y[i]+=g.dt*u.py[i];u.s[i]+=g.dt*u.ps[i];}
  force(u,g,fx,fy,fs);
  for(int i=0;i<g.n;i++){u.px[i]+=.5*g.dt*fx[i];u.py[i]+=.5*g.dt*fy[i];u.ps[i]+=.5*g.dt*fs[i];}
  damping(u,g,ea,qa);
 }
 std::ofstream final(prefix+"_final.csv");final<<std::setprecision(17)<<"rho,x,y,s,px,py,ps\n";
 for(int i=0;i<g.n;i++)final<<g.r[i]<<','<<u.x[i]<<','<<u.y[i]<<','<<u.s[i]<<','<<u.px[i]<<','<<u.py[i]<<','<<u.ps[i]<<'\n';
 history.flush();trace.flush();pressure.flush();final.flush();
 if(!history||!trace||!pressure||!final)throw std::runtime_error("output stream failed before completion");
 auto b=budget(u,g);std::cout<<std::setprecision(17)<<"Final E "<<b.E<<" Q "<<b.Q<<" Eabs "<<ea<<" Qabs "<<qa<<" Fmax "<<b.xmax<<" Smin "<<b.smin<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
