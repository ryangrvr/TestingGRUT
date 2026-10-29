// Even-Legendre Galerkin, axisymmetric M1 scalar evolution.
// Same radial finite volumes, splitting and damping as parent radial.cpp.
// Angular quadrature is required to integrate quartic modal products.
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
struct Grid {int n,m,nq;double h,dt,L;vector<int>l;vector<double>r,vol,face,damp,basis,proj,w;};
struct Budget {double E=0,Q=0,E20=0,Q20=0,Eabs=0,Qabs=0;};

double potential(double x,double y,double s){double z=x*x+y*y;return .75*(s*s-1)*(s*s-1)-.2*std::pow(s-1,3)+.5*s*s*z+.025*z*z;}

void force(const State& u,const Grid& g,vector<double>& fx,vector<double>& fy,vector<double>& fs){
 const int m=g.m,nq=g.nq;
 for(int i=0;i<g.n;i++){
  std::array<double,64> a{},b{},c{};
  for(int j=0;j<nq;j++)for(int k=0;k<m;k++){
   double P=g.basis[j*m+k];size_t z=(size_t)i*m+k;
   a[j]+=P*u.x[z];b[j]+=P*u.y[z];c[j]+=P*u.s[z];
  }
  for(int j=0;j<nq;j++){
   double z=a[j]*a[j]+b[j]*b[j],f=c[j]*c[j]+.1*z;
   double vx=-f*a[j],vy=-f*b[j],vs=-(3*(c[j]*c[j]-1)+z)*c[j]+.6*(c[j]-1)*(c[j]-1);
   a[j]=vx;b[j]=vy;c[j]=vs;
  }
  for(int k=0;k<m;k++){
   size_t z=(size_t)i*m+k;double xnext=i+1<g.n?u.x[z+m]:0,ynext=i+1<g.n?u.y[z+m]:0,snext=i+1<g.n?u.s[z+m]:(k==0?1:0);
   double lx=g.face[i]*(xnext-u.x[z]),ly=g.face[i]*(ynext-u.y[z]),ls=g.face[i]*(snext-u.s[z]);
   if(i){lx-=g.face[i-1]*(u.x[z]-u.x[z-m]);ly-=g.face[i-1]*(u.y[z]-u.y[z-m]);ls-=g.face[i-1]*(u.s[z]-u.s[z-m]);}
   double angular=4*PI*g.h*g.l[k]*(g.l[k]+1)/g.vol[i];
   fx[z]=lx/g.vol[i]-angular*u.x[z];fy[z]=ly/g.vol[i]-angular*u.y[z];fs[z]=ls/g.vol[i]-angular*u.s[z];
   for(int j=0;j<nq;j++){double p=g.proj[j*m+k];fx[z]+=p*a[j];fy[z]+=p*b[j];fs[z]+=p*c[j];}
  }
 }
}

void damp(State& u,const Grid& g,double& energy,double& charge){
 for(int i=0;i<g.n;i++)if(g.damp[i]<1){
  double f=g.damp[i],E=0,Q=0;
  for(int k=0;k<g.m;k++){size_t z=(size_t)i*g.m+k;double w=1./(2*g.l[k]+1);
   E+=.5*w*(u.px[z]*u.px[z]+u.py[z]*u.py[z]+u.ps[z]*u.ps[z]);
   Q+=w*(u.y[z]*u.px[z]-u.x[z]*u.py[z]);
   u.px[z]*=f;u.py[z]*=f;u.ps[z]*=f;
  }
  energy+=g.vol[i]*(1-f*f)*E;charge+=g.vol[i]*(1-f)*Q;
 }
}

Budget budget(const State& u,const Grid& g){
 Budget B;
 for(int i=0;i<g.n;i++){
  double kin=0,grad=0,ang=0,Q=0,pot=0;
  for(int k=0;k<g.m;k++){size_t z=(size_t)i*g.m+k;double w=1./(2*g.l[k]+1);
   kin+=.5*w*(u.px[z]*u.px[z]+u.py[z]*u.py[z]+u.ps[z]*u.ps[z]);Q+=w*(u.y[z]*u.px[z]-u.x[z]*u.py[z]);
   double xn=i+1<g.n?u.x[z+g.m]:0,yn=i+1<g.n?u.y[z+g.m]:0,sn=i+1<g.n?u.s[z+g.m]:(k==0?1:0);
   grad+=.5*g.face[i]*w*((xn-u.x[z])*(xn-u.x[z])+(yn-u.y[z])*(yn-u.y[z])+(sn-u.s[z])*(sn-u.s[z]));
   ang+=.5*4*PI*g.h*w*g.l[k]*(g.l[k]+1)*(u.x[z]*u.x[z]+u.y[z]*u.y[z]+u.s[z]*u.s[z]);
  }
  for(int j=0;j<g.nq;j++){
   double x=0,y=0,s=0;for(int k=0;k<g.m;k++){size_t z=(size_t)i*g.m+k;double p=g.basis[j*g.m+k];x+=p*u.x[z];y+=p*u.y[z];s+=p*u.s[z];}
   pot+=.5*g.w[j]*potential(x,y,s);
  }
  double E=g.vol[i]*(kin+pot)+grad+ang,q=g.vol[i]*Q;B.E+=E;B.Q+=q;
  if(g.r[i]<20){B.E20+=E;B.Q20+=q;}
 }
 return B;
}

int main(int argc,char**argv){try{
 if(argc!=13)throw std::runtime_error("usage axis INITIAL ANGLES PREFIX h dt L tmax probe absorber_start eta lmax nq");
 std::string initial=argv[1],angles=argv[2],prefix=argv[3];Grid g;
 g.h=std::stod(argv[4]);g.dt=std::stod(argv[5]);g.L=std::stod(argv[6]);double T=std::stod(argv[7]),probe=std::stod(argv[8]),ab=std::stod(argv[9]),eta=std::stod(argv[10]);
 int lmax=std::stoi(argv[11]);g.nq=std::stoi(argv[12]);g.m=lmax/2+1;g.n=std::lround(g.L/g.h);
 if(lmax%2||g.nq>64||g.nq<2*lmax+1||g.n<4)throw std::runtime_error("invalid angular or radial grid");
 g.l.resize(g.m);for(int k=0;k<g.m;k++)g.l[k]=2*k;
 g.basis.resize(g.nq*g.m);g.proj.resize(g.nq*g.m);g.w.resize(g.nq);
 std::ifstream ain(angles);std::string line;std::getline(ain,line);
 for(int j=0;j<g.nq;j++){if(!std::getline(ain,line))throw std::runtime_error("angle file short");std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);double mu;s>>mu>>g.w[j];
  for(int k=0;k<g.m;k++){double p;s>>p;g.basis[j*g.m+k]=p;g.proj[j*g.m+k]=(2*g.l[k]+1)*.5*g.w[j]*p;}}
 size_t ns=(size_t)g.n*g.m;State u;for(auto*p:{&u.x,&u.y,&u.s,&u.px,&u.py,&u.ps})p->resize(ns);
 g.r.resize(g.n);g.vol.resize(g.n);g.face.resize(g.n);g.damp.resize(g.n);
 std::ifstream in(initial);std::getline(in,line);
 for(int i=0;i<g.n;i++){
  if(!std::getline(in,line))throw std::runtime_error("initial file short");
  std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);double r;s>>r;
  for(auto*p:{&u.x,&u.y,&u.s,&u.px,&u.py,&u.ps})for(int k=0;k<g.m;k++)s>>(*p)[(size_t)i*g.m+k];
  g.r[i]=(i+.5)*g.h;if(std::abs(r-g.r[i])>1e-8)throw std::runtime_error("grid mismatch");
  double lo=i*g.h,hi=(i+1)*g.h;g.vol[i]=4*PI/3*(hi*hi*hi-lo*lo*lo);g.face[i]=4*PI*hi*hi/g.h*(i+1==g.n?2:1);
  double z=std::max(0.,(g.r[i]-ab)/(g.L-ab));g.damp[i]=std::exp(-.5*g.dt*eta*z*z*z);
 }
 vector<double>fx(ns),fy(ns),fs(ns);double ea=0,qa=0;int ip=std::clamp((int)std::lround(probe/g.h-.5),1,g.n-2);
 long steps=std::lround(T/g.dt),histstride=std::max<long>(1,std::lround(10/g.dt)),tracestride=std::max<long>(1,std::lround(.25/g.dt));
 std::ofstream hist(prefix+"_history.csv"),trace(prefix+"_trace.csv");hist<<std::setprecision(17);trace<<std::setprecision(17);
 hist<<"tau,E,Q,Eabs,Qabs,E20,Q20\n";trace<<"tau";
 for(int l:g.l)for(auto p:{"x","y","s_minus1","dx","dy","ds","px","py","ps"}){trace<<','<<p<<"_l"<<l;}
 trace<<'\n';
 for(long step=0;step<=steps;step++){
  double t=step*g.dt;
  if(step%histstride==0||step==steps){auto b=budget(u,g);hist<<t<<','<<b.E<<','<<b.Q<<','<<ea<<','<<qa<<','<<b.E20<<','<<b.Q20<<'\n';if(!std::isfinite(b.E))throw std::runtime_error("nonfinite energy");}
  if(step%tracestride==0||step==steps){trace<<t;for(int k=0;k<g.m;k++){size_t z=(size_t)ip*g.m+k;double dx=(u.x[z+g.m]-u.x[z-g.m])/(2*g.h),dy=(u.y[z+g.m]-u.y[z-g.m])/(2*g.h),ds=(u.s[z+g.m]-u.s[z-g.m])/(2*g.h);
   trace<<','<<u.x[z]<<','<<u.y[z]<<','<<u.s[z]-(k==0?1:0)<<','<<dx<<','<<dy<<','<<ds<<','<<u.px[z]<<','<<u.py[z]<<','<<u.ps[z];}trace<<'\n';}
  if(step==steps)break;
  damp(u,g,ea,qa);force(u,g,fx,fy,fs);
  for(size_t z=0;z<ns;z++){u.px[z]+=.5*g.dt*fx[z];u.py[z]+=.5*g.dt*fy[z];u.ps[z]+=.5*g.dt*fs[z];}
  for(size_t z=0;z<ns;z++){u.x[z]+=g.dt*u.px[z];u.y[z]+=g.dt*u.py[z];u.s[z]+=g.dt*u.ps[z];}
  force(u,g,fx,fy,fs);
  for(size_t z=0;z<ns;z++){u.px[z]+=.5*g.dt*fx[z];u.py[z]+=.5*g.dt*fy[z];u.ps[z]+=.5*g.dt*fs[z];}
  damp(u,g,ea,qa);
 }
 hist.flush();trace.flush();if(!hist||!trace)throw std::runtime_error("output stream failed before completion");
 auto b=budget(u,g);std::cout<<std::setprecision(17)<<"E "<<b.E<<" Q "<<b.Q<<" Eabs "<<ea<<" Qabs "<<qa<<'\n';
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
