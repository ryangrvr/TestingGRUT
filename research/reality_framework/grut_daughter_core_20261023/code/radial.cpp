// Spherical finite-volume X/phi Hamiltonian prototype, no quantum decay term.
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

Budget budget(const State& u,const Geometry& g){
 Budget out;const std::array<double,3> cuts={40.,60.,80.};
 for(int i=0;i<g.n;i++){
  double n=u.x[i]*u.x[i]+u.y[i]*u.y[i];
  double e=g.vol[i]*(.5*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i])
       +.75*std::pow(u.s[i]*u.s[i]-1,2)+.5*u.s[i]*u.s[i]*n+.025*n*n);
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
  fs[i]=ls/g.vol[i]-(3*(u.s[i]*u.s[i]-1)+n)*u.s[i];
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
 if(argc!=10)throw std::runtime_error("usage radial INITIAL PREFIX h dt L tmax probe absorber_start eta");
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
 std::ofstream history(prefix+"_history.csv"),trace(prefix+"_trace.csv");
 history<<std::setprecision(17);trace<<std::setprecision(17);
 history<<"tau,E,Q,E_abs,Q_abs,Eflux_probe,Qflux_probe,Fmax,Smin,E40,Q40,E60,Q60,E80,Q80\n";
 trace<<"tau,x,y,s_minus1,dx,dy,ds,px,py,ps\n";
 int ip=std::lround(probe/g.h-.5);ip=std::clamp(ip,1,g.n-2);
 int steps=std::lround(end/g.dt),histstride=std::max(1,(int)std::lround(10/g.dt));
 int tracestride=std::max(1,(int)std::lround(.25/g.dt));
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
   for(int j=0;j<3;j++)history<<','<<b.ec[j]<<','<<b.qc[j];history<<'\n';
  }
  if(step%tracestride==0||step==steps)trace<<t<<','<<u.x[ip]<<','<<u.y[ip]<<','<<u.s[ip]-1<<','<<d[0]<<','<<d[1]<<','<<d[2]<<','<<u.px[ip]<<','<<u.py[ip]<<','<<u.ps[ip]<<'\n';
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
 auto b=budget(u,g);std::cout<<std::setprecision(17)<<"Final E "<<b.E<<" Q "<<b.Q<<" Eabs "<<ea<<" Qabs "<<qa<<" Fmax "<<b.xmax<<" Smin "<<b.smin<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
