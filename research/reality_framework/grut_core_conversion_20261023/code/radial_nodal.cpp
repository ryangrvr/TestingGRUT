// Fourth-order radial Hamiltonian in U=r*psi, Z=r*(S-1).
// Odd endpoint continuation implements regularity. No physical damping term.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
using std::vector;
constexpr double PI=3.14159265358979323846;
struct State{vector<double>x,y,s,px,py,ps;};
double sample(const vector<double>&u,int j){
 int n=u.size();if(j==-1||j==n)return 0.;
 if(j==-2)return -u[0];if(j==n+1)return -u[n-1];return u[j];
}
double lap(const vector<double>&u,int i,double h){
 return (-sample(u,i-2)+16*sample(u,i-1)-30*u[i]+16*sample(u,i+1)-sample(u,i+2))/(12*h*h);
}
void force(const State&u,double h,vector<double>&fx,vector<double>&fy,vector<double>&fs){
 for(int i=0;i<(int)u.x.size();i++){
  double r=(i+1)*h,n=(u.x[i]*u.x[i]+u.y[i]*u.y[i])/(r*r),s=1+u.s[i]/r;
  fx[i]=lap(u.x,i,h)-(s*s+.1*n)*u.x[i];
  fy[i]=lap(u.y,i,h)-(s*s+.1*n)*u.y[i];
  fs[i]=lap(u.s,i,h)-r*(3*(s*s-1)+n)*s;
 }
}
struct Budget{double E=0,Q=0,Fmax=0,Smin=1;std::array<double,3>ec{},qc{};};
Budget budget(const State&u,double h){
 Budget b;std::array<double,3>cuts={40,60,80};double w=4*PI*h;
 for(int i=0;i<(int)u.x.size();i++){
  double r=(i+1)*h,n=(u.x[i]*u.x[i]+u.y[i]*u.y[i])/(r*r),s=1+u.s[i]/r;
  double e=w*(.5*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i])
       -.5*(u.x[i]*lap(u.x,i,h)+u.y[i]*lap(u.y,i,h)+u.s[i]*lap(u.s,i,h))
       +r*r*(.75*(s*s-1)*(s*s-1)+.5*s*s*n+.025*n*n));
  double q=w*(u.y[i]*u.px[i]-u.x[i]*u.py[i]);b.E+=e;b.Q+=q;
  b.Fmax=std::max(b.Fmax,std::sqrt(n));b.Smin=std::min(b.Smin,s);
  for(int j=0;j<3;j++){
   double cutw=r<cuts[j]-h*.1?1.:(std::abs(r-cuts[j])<h*.1?.5:0.);
   b.ec[j]+=cutw*e;b.qc[j]+=cutw*q;
  }
 }return b;
}
int main(int argc,char**argv){try{
 if(argc!=10)throw std::runtime_error("usage INITIAL PREFIX h dt L tmax probe absorber eta");
 std::string prefix=argv[2];double h=std::stod(argv[3]),dt=std::stod(argv[4]),L=std::stod(argv[5]);
 double end=std::stod(argv[6]),probe=std::stod(argv[7]),ab=std::stod(argv[8]),eta=std::stod(argv[9]);
 int n=std::lround(L/h)-1;State u;for(auto*p:{&u.x,&u.y,&u.s,&u.px,&u.py,&u.ps})p->resize(n);
 vector<double>fx(n),fy(n),fs(n),damp(n);std::ifstream init(argv[1]);std::string line;std::getline(init,line);
 for(int i=0;i<n;i++){
  if(!std::getline(init,line))throw std::runtime_error("short initial grid");
  std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double r,x,y,s,px,py,ps;
  row>>r>>x>>y>>s>>px>>py>>ps;if(std::abs(r-(i+1)*h)>1e-8)throw std::runtime_error("grid mismatch");
  u.x[i]=r*x;u.y[i]=r*y;u.s[i]=r*(s-1);u.px[i]=r*px;u.py[i]=r*py;u.ps[i]=r*ps;
  double z=std::max(0.,(r-ab)/(L-ab));damp[i]=std::exp(-.5*dt*eta*z*z*z);
 }
 double ea=0,qa=0,eflux=0,qflux=0,prevE=0,prevQ=0;
 auto damping=[&](){for(int i=0;i<n;i++){
  double f=damp[i],w=4*PI*h;ea+=.5*w*(1-f*f)*(u.px[i]*u.px[i]+u.py[i]*u.py[i]+u.ps[i]*u.ps[i]);
  qa+=w*(1-f)*(u.y[i]*u.px[i]-u.x[i]*u.py[i]);u.px[i]*=f;u.py[i]*=f;u.ps[i]*=f;
 }};
 auto vv=[&](double d){force(u,h,fx,fy,fs);
  for(int i=0;i<n;i++){u.px[i]+=.5*d*fx[i];u.py[i]+=.5*d*fy[i];u.ps[i]+=.5*d*fs[i];}
  for(int i=0;i<n;i++){u.x[i]+=d*u.px[i];u.y[i]+=d*u.py[i];u.s[i]+=d*u.ps[i];}
  force(u,h,fx,fy,fs);
  for(int i=0;i<n;i++){u.px[i]+=.5*d*fx[i];u.py[i]+=.5*d*fy[i];u.ps[i]+=.5*d*fs[i];}
 };
 std::ofstream hist(prefix+"_history.csv"),trace(prefix+"_trace.csv");hist<<std::setprecision(17);trace<<std::setprecision(17);
 hist<<"tau,E,Q,E_abs,Q_abs,Eflux_probe,Qflux_probe,Fmax,Smin,E40,Q40,E60,Q60,E80,Q80\n";
 trace<<"tau,x,y,s_minus1,dx,dy,ds,px,py,ps\n";
 int ip=std::clamp((int)std::lround(probe/h)-1,1,n-2),steps=std::lround(end/dt);
 int hs=std::max(1,(int)std::lround(10/dt)),ts=std::max(1,(int)std::lround(.25/dt));
 for(int step=0;step<=steps;step++){
  double r=(ip+1)*h,t=step*dt,x=u.x[ip]/r,y=u.y[ip]/r,s=u.s[ip]/r;
  double dx=(u.x[ip+1]/(r+h)-u.x[ip-1]/(r-h))/(2*h),dy=(u.y[ip+1]/(r+h)-u.y[ip-1]/(r-h))/(2*h);
  double ds=(u.s[ip+1]/(r+h)-u.s[ip-1]/(r-h))/(2*h),px=u.px[ip]/r,py=u.py[ip]/r,ps=u.ps[ip]/r;
  double ef=-4*PI*r*r*(px*dx+py*dy+ps*ds),qf=4*PI*r*r*(x*dy-y*dx);
  if(step){eflux+=.5*dt*(ef+prevE);qflux+=.5*dt*(qf+prevQ);}prevE=ef;prevQ=qf;
  if(step%hs==0||step==steps){auto b=budget(u,h);if(!std::isfinite(b.E)||b.Fmax>100)throw std::runtime_error("runaway");
   hist<<t<<','<<b.E<<','<<b.Q<<','<<ea<<','<<qa<<','<<eflux<<','<<qflux<<','<<b.Fmax<<','<<b.Smin;
   for(int j=0;j<3;j++)hist<<','<<b.ec[j]<<','<<b.qc[j];hist<<'\n';}
  if(step%ts==0||step==steps)trace<<t<<','<<x<<','<<y<<','<<s<<','<<dx<<','<<dy<<','<<ds<<','<<px<<','<<py<<','<<ps<<'\n';
  if(step==steps)break;damping();double a=1/(2-std::cbrt(2.)),b=-std::cbrt(2.)*a;vv(a*dt);vv(b*dt);vv(a*dt);damping();
 }
 std::ofstream out(prefix+"_final.csv");out<<std::setprecision(17)<<"rho,x,y,s,px,py,ps\n";
 for(int i=0;i<n;i++){double r=(i+1)*h;out<<r<<','<<u.x[i]/r<<','<<u.y[i]/r<<','<<1+u.s[i]/r<<','<<u.px[i]/r<<','<<u.py[i]/r<<','<<u.ps[i]/r<<'\n';}
 auto b=budget(u,h);std::cout<<std::setprecision(17)<<"Final E "<<b.E<<" Q "<<b.Q<<" Eabs "<<ea<<" Qabs "<<qa<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
