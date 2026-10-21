// Independent extension of Codex's 3-real-field Cartesian discretization.
// D=0, externally prescribed radiation expansion. No PRS modification.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
using Vec=std::vector<double>;
struct State {std::array<Vec,3> q,p; explicit State(size_t v){for(auto* b:{&q,&p})for(auto& f:*b)f.assign(v,0.);}};
struct Grid {
 int n; size_t n2,v; double h,L;
 Grid(int N,double dx):n(N),n2(size_t(N)*N),v(n2*N),h(dx),L(N*dx){}
 std::array<size_t,6> nb(size_t z)const {
  int i=z/n2,j=(z%n2)/n,k=z%n;
  return {i?z-n2:z+(n-1)*n2,i+1<n?z+n2:z-(n-1)*n2,
   j?z-n:z+size_t(n-1)*n,j+1<n?z+n:z-size_t(n-1)*n,
   k?z-1:z+n-1,k+1<n?z+1:z-n+1};
 }
 double kick(State& b,double dt,double a,double rate)const {
  double source=0,ih2=1/(h*h);
  #pragma omp parallel for collapse(2) reduction(+:source) schedule(static)
  for(int i=0;i<n;i++)for(int j=0;j<n;j++){
   size_t base=(size_t(i)*n+j)*n,ip=size_t(i+1<n?i+1:0)*n2+j*n,
    im=size_t(i?i-1:n-1)*n2+j*n,jp=size_t(i)*n2+(j+1<n?j+1:0)*n,
    jm=size_t(i)*n2+(j?j-1:n-1)*n;
   for(int k=0;k<n;k++){
    size_t z=base+k,kp=base+(k+1<n?k+1:0),km=base+(k?k-1:n-1);
    double u=b.q[0][z],w=b.q[1][z],s=b.q[2][z],f=u*u+w*w;
    std::array<double,3> pot={(s*s+.1*f)*u,(s*s+.1*f)*w,(s*s-a*a+f)*s};
    for(int c=0;c<3;c++){
     const auto& q=b.q[c];double lap=(q[ip+k]+q[im+k]+q[jp+k]+q[jm+k]+q[kp]+q[km]-6*q[z])*ih2;
     b.p[c][z]+=dt*(lap-pot[c]);
    }
    source-=a*rate*(s*s-a*a);
   }
  }
  return source*h*h*h;
 }
 void drift(State& b,double dt)const {
  #pragma omp parallel for schedule(static)
  for(size_t z=0;z<v;z++)for(int c=0;c<3;c++)b.q[c][z]+=dt*b.p[c][z];
 }
};
struct Modes {
 std::vector<std::array<double,9>> data;
 explicit Modes(const std::string& path){if(path.empty())return;std::ifstream f(path);std::array<double,9> r;while(f>>r[0]>>r[1]>>r[2]>>r[3]>>r[4]>>r[5]>>r[6]>>r[7]>>r[8])data.push_back(r);if(data.empty())throw std::runtime_error("empty mode table");}
};
void initial(const Grid& g,State& b,const Modes& modes,std::string kind,double kap,double rate,double noise,int mirror,double speed){
 double x=std::pow(kap,2./3.),S0=0,omega=std::sqrt(.1*x),F=std::sqrt(x);
 bool cold=kind=="cold";
 if(cold){double lo=0,hi=2/2.7,target=kap*std::sqrt(.1);for(int it=0;it<100;it++){double m=.5*(lo+hi);if(m*std::sqrt(1-.9*m)>target)hi=m;else lo=m;}x=.5*(lo+hi);S0=std::sqrt(1-x);omega=std::sqrt(1-.9*x);F=std::sqrt(x);}
 if(kind=="blob"){
  #pragma omp parallel for schedule(static)
  for(size_t z=0;z<g.v;z++){int ix=z/g.n2,iy=(z%g.n2)/g.n,iz=z%g.n;double rr=0;for(int i:{ix,iy,iz}){double x=(i+.5)*g.h;x=std::min(x,g.L-x);rr+=x*x;}b.q[2][z]=std::tanh((g.L/8-std::sqrt(rr))/std::sqrt(2.));}
  return;
 }
 if(kind=="kink"){
  double gam=1/std::sqrt(1-speed*speed),d=std::sqrt(2.);
  #pragma omp parallel for collapse(2) schedule(static)
  for(int i=0;i<g.n;i++)for(int j=0;j<g.n;j++)for(int k=0;k<g.n;k++){
   size_t z=(size_t(i)*g.n+j)*g.n+k;double xx=(i+.5)*g.h;
   double A=std::tanh(gam*(xx-.25*g.L)/d),B=std::tanh(gam*(xx-.75*g.L)/d);
   b.q[2][z]=A*B;b.p[2][z]=-speed*gam/d*((1-A*A)*B+(1-B*B)*A);
  }return;
 }
 for(auto row:modes.data){
  std::array<Vec,3> co,si;
  for(int ax=0;ax<3;ax++){co[ax].resize(g.n);si[ax].resize(g.n);for(int i=0;i<g.n;i++){double t=2*std::acos(-1.)*row[ax]*(i+.5)/g.n;co[ax][i]=std::cos(t);si[ax][i]=std::sin(t);}}
  #pragma omp parallel for collapse(2) schedule(static)
  for(int i=0;i<g.n;i++)for(int j=0;j<g.n;j++){
   double cxy=co[0][i]*co[1][j]-si[0][i]*si[1][j],sxy=si[0][i]*co[1][j]+co[0][i]*si[1][j];
   for(int k=0;k<g.n;k++){size_t z=(size_t(i)*g.n+j)*g.n+k;double c=cxy*co[2][k]-sxy*si[2][k],s=sxy*co[2][k]+cxy*si[2][k];for(int ax=0;ax<3;ax++)b.q[ax][z]+=row[3+2*ax]*c+row[4+2*ax]*s;}
  }
 }
 if(kind=="homogeneous"||kind=="linear")for(auto& f:b.q)std::fill(f.begin(),f.end(),0.);
 double meanf=0;
 #pragma omp parallel for schedule(static)
 for(size_t z=0;z<g.v;z++){
  b.q[0][z]=F*(1+.001*noise*b.q[0][z]);b.q[1][z]=F*.001*noise*b.q[1][z];
  if(kind=="linear"){int i=z/g.n2;b.q[2][z]=noise*std::cos(2*std::acos(-1.)*(i+.5)/g.n);}
  else b.q[2][z]=S0+.01*noise*b.q[2][z];
  b.q[2][z]*=mirror;
 }
 // Deterministic normalization: parallel reduction order otherwise changes the
 // initial X field by ulps between the two members of the exact mirror control.
 for(size_t z=0;z<g.v;z++)meanf+=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z];
 double scale=std::sqrt(x*g.v/meanf);
 #pragma omp parallel for schedule(static)
 for(size_t z=0;z<g.v;z++){
  b.q[0][z]*=scale;b.q[1][z]*=scale;
  b.p[0][z]=omega*b.q[1][z]+(cold?rate*b.q[0][z]:0);
  b.p[1][z]=-omega*b.q[0][z]+(cold?rate*b.q[1][z]:0);
  b.p[2][z]=rate*b.q[2][z];
 }
}
struct Topology {size_t occupied=0,largest=0,components=0;int any=0,large=0;};
Topology topology(const Grid& g,const std::vector<int8_t>& mask,int sign){
 Topology out;const int unset=std::numeric_limits<int>::min();
 std::vector<int> lx(g.v,unset),ly(g.v),lz(g.v);std::vector<uint32_t> queue;
 for(size_t start=0;start<g.v;start++){
  if(mask[start]!=sign||lx[start]!=unset)continue;
  queue.clear();queue.push_back(start);lx[start]=0;ly[start]=0;lz[start]=0;int winding=0;
  for(size_t a=0;a<queue.size();a++){
   size_t z=queue[a];auto ns=g.nb(z);
   for(int d=0;d<6;d++){
    size_t b=ns[d];if(mask[b]!=sign)continue;
    int xx=lx[z]+(d==0?-1:d==1?1:0),yy=ly[z]+(d==2?-1:d==3?1:0),zz=lz[z]+(d==4?-1:d==5?1:0);
    if(lx[b]==unset){lx[b]=xx;ly[b]=yy;lz[b]=zz;queue.push_back(b);}
    else {if(xx!=lx[b])winding|=1;if(yy!=ly[b])winding|=2;if(zz!=lz[b])winding|=4;}
   }
  }
  out.components++;out.occupied+=queue.size();out.any|=winding;
  if(queue.size()>out.largest){out.largest=queue.size();out.large=winding;}
 }
 return out;
}
void print_top(std::ostream& o,const Topology& t,size_t v){o<<"{\"components\":"<<t.components<<",\"occupied_fraction\":"<<double(t.occupied)/v<<",\"largest_fraction\":"<<double(t.largest)/v<<",\"any_winding\":"<<t.any<<",\"largest_winding\":"<<t.large<<"}";}
void topo_diag(const Grid& g,const State& b,double t,double a,std::ostream& o){
 std::vector<int8_t> labels(g.v);std::vector<uint32_t> queue;queue.reserve(g.v);
 for(size_t z=0;z<g.v;z++){
  double S=b.q[2][z]/a,F2=(b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z])/(a*a);
  if(F2<.5&&std::abs(S)>.5){labels[z]=S>0?1:-1;queue.push_back(z);}
 }
 auto p=topology(g,labels,1),m=topology(g,labels,-1);
 // Nearest seed in the graph distance. Seed order depends on position, not sign.
 for(size_t h=0;h<queue.size();h++){size_t z=queue[h];for(auto nb:g.nb(z))if(labels[nb]==0){labels[nb]=labels[z];queue.push_back(nb);}}
 size_t links=0;std::vector<uint8_t> shell(g.v,255);queue.clear();
 for(size_t z=0;z<g.v;z++){
  auto ns=g.nb(z);bool boundary=false;
  for(int d=0;d<6;d++)if(labels[z]&&labels[ns[d]]&&labels[z]!=labels[ns[d]]){boundary=true;if(d%2==1)links++;}
  if(boundary){shell[z]=0;queue.push_back(z);}
 }
 int radius=int(std::ceil(2/(a*g.h)));
 for(size_t h=0;h<queue.size();h++){size_t z=queue[h];if(shell[z]>=radius)continue;for(auto nb:g.nb(z))if(shell[nb]==255){shell[nb]=shell[z]+1;queue.push_back(nb);}}
 double Q=0,Qshell=0;for(size_t z=0;z<g.v;z++){double q=b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z];Q+=q;if(shell[z]!=255)Qshell+=q;}
 o<<std::setprecision(17)<<"{\"t\":"<<t<<",\"a\":"<<a<<",\"plus\":";print_top(o,p,g.v);o<<",\"minus\":";print_top(o,m,g.v);
 o<<",\"basin_link_area_comoving_over_V\":"<<(2./3)*links/(g.v*g.h)<<",\"interface_links\":"<<links<<",\"q_fraction_interface_shell\":"<<(Q?Qshell/Q:0)<<",\"shell_radius_cells\":"<<radius<<"}\n";o.flush();
}
struct Diagnostics {double H=0,Q=0,rho=0,P=0;};
Diagnostics scalar(const Grid& g,const State& b,double t,double a,double rate,double work,std::ostream& out){
 double K=0,G=0,U=0,HaK=0,Q=0,qneg=0,qcore=0,qF=0,source=0,Smean=0,Smin=1e99,Smax=-1e99;
 double fracp[3]={},fracm[3]={},genuinep=0,genuinem=0,corearea=0,coreEphi=0,corepi=0,coregrad=0,linearproj=0;
 double gradS=0,kinS=0,Uphi=0,Ucross=0,Ux=0,smootharea=0,smoothEphi=0,smoothpi=0,smoothgrad=0;
 #pragma omp parallel for reduction(+:K,G,U,HaK,Q,qneg,qcore,qF,source,Smean,genuinep,genuinem,corearea,coreEphi,corepi,coregrad,linearproj,gradS,kinS,Uphi,Ucross,Ux,smootharea,smoothEphi,smoothpi,smoothgrad,fracp[:3],fracm[:3]) reduction(min:Smin) reduction(max:Smax) schedule(static)
 for(size_t z=0;z<g.v;z++){
  auto ns=g.nb(z);double f=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z],s=b.q[2][z],sp=s/a;
  double us=.25*(s*s-a*a)*(s*s-a*a),uc=.5*s*s*f,ux=.025*f*f;
  double k=0,gg=0,ks=0,gs=0,gsc=0;
  for(int c=0;c<3;c++){
   double pp=b.p[c][z]-rate/a*b.q[c][z];k+=.5*pp*pp;HaK+=.5*b.p[c][z]*b.p[c][z];
   for(int ax=0;ax<3;ax++){double df=(b.q[c][ns[2*ax+1]]-b.q[c][z])/g.h;gg+=.5*df*df;if(c==2){gs+=.5*df*df;double d=(b.q[c][ns[2*ax+1]]-b.q[c][ns[2*ax]])/(2*g.h);gsc+=d*d;}}
   if(c==2)ks=.5*pp*pp;
  }
  K+=k;G+=gg;U+=us+uc+ux;kinS+=ks;gradS+=gs;Uphi+=us;Ucross+=uc;Ux+=ux;
  double q=b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z];Q+=q;qneg+=std::max(0.,-q);
  if(std::abs(sp)<.5){qcore+=q;corearea+=std::sqrt(gsc)/a;coreEphi+=ks+gs+us;corepi+=2*ks;coregrad+=gsc;}
  if(std::abs(sp)<1){double weight=1-sp*sp;smootharea+=.75*weight*std::sqrt(gsc)/a;smoothEphi+=weight*(ks+gs+us);smoothpi+=weight*2*ks;smoothgrad+=weight*gsc;}
  if(f/(a*a)>1)qF+=q;
  if(f/(a*a)<.5){if(sp>.5)genuinep++;if(sp<-.5)genuinem++;}
  for(int i=0;i<3;i++){double cut=.4+.1*i;if(sp>cut)fracp[i]++;if(sp<-cut)fracm[i]++;}
  Smean+=sp;Smin=std::min(Smin,sp);Smax=std::max(Smax,sp);
  source-=a*rate*(s*s-a*a);int i=z/g.n2;linearproj+=2*s*std::cos(2*std::acos(-1.)*(i+.5)/g.n);
 }
 double dv=g.h*g.h*g.h,V=g.v*dv,norm=1/(g.v*std::pow(a,4)),rho=(K+G+U)*norm,pressure=(K-G/3-U)*norm;
 double Hconf=(HaK+G+U)*dv;
 out<<std::setprecision(17)<<t<<','<<a<<','<<Hconf<<','<<work<<','<<source*dv<<','<<Q*dv<<','<<qneg*dv<<','<<rho<<','<<pressure<<','<<K*norm<<','<<G*norm<<','<<U*norm<<','<<kinS*norm<<','<<gradS*norm<<','<<Uphi*norm<<','<<Ucross*norm<<','<<Ux*norm<<','<<Smean/g.v<<','<<Smin<<','<<Smax;
 for(int i=0;i<3;i++)out<<','<<fracp[i]/g.v<<','<<fracm[i]/g.v;
 out<<','<<genuinep/g.v<<','<<genuinem/g.v<<','<<corearea/g.v<<','<<coreEphi*dv/(a*.6875)<<','<<(coregrad?corepi/coregrad:0)<<','<<(Q?qcore/Q:0)<<','<<(Q?qF/Q:0)<<','<<linearproj/g.v<<','<<smootharea/g.v<<','<<smoothEphi*dv/(a*.8)<<','<<(smoothgrad?smoothpi/smoothgrad:0)<<'\n';out.flush();
 return {Hconf,Q*dv,rho,pressure};
}
void slice(const Grid& g,const State& b,double a,const std::string& name){
 std::ofstream f(name,std::ios::binary);for(int i=0;i<g.n;i++)for(int j=0;j<g.n;j++){
  size_t z=(size_t(i)*g.n+j)*g.n+g.n/2;float S=b.q[2][z]/a,F2=(b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z])/(a*a),q=(b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z])/std::pow(a,3);
  float row[]={S,F2,q};f.write(reinterpret_cast<const char*>(row),sizeof(row));
 }
}
int main(int argc,char** argv){try{
 std::map<std::string,std::string> kv;for(int i=1;i<argc;i+=2){if(i+1==argc)throw std::runtime_error("key/value arguments required");kv[argv[i]]=argv[i+1];}
 auto val=[&](std::string k,std::string d){return kv.count(k)?kv[k]:d;};
 int n=std::stoi(val("--N","32")),threads=std::stoi(val("--threads","8")),mirror=std::stoi(val("--mirror","1"));
 double dx=std::stod(val("--dx","0.4")),dt=std::stod(val("--dt","0.08")),T=std::stod(val("--T","48")),taui=std::stod(val("--taui","40")),rate=kv.count("--rate")?std::stod(kv["--rate"]):1/taui,kap=std::stod(val("--kappa","2")),noise=std::stod(val("--noise","1")),speed=std::stod(val("--speed","0.456"));
 if(n<4||uint64_t(n)*n*n>std::numeric_limits<uint32_t>::max())throw std::runtime_error("grid out of range");
 int steps=std::llround(T/dt),stride=std::max(1,int(std::llround(std::stod(val("--sample","1"))/dt))),topstride=std::max(1,int(std::llround(std::stod(val("--top-sample","4"))/dt)));
 if(std::abs(steps*dt-T)>1e-9)throw std::runtime_error("T/dt must be integral");
 omp_set_num_threads(threads);Grid g(n,dx);State b(g.v);std::string kind=val("--kind","charge"),prefix=val("--out","run");
 Modes modes(val("--modes",""));auto start=std::chrono::steady_clock::now();initial(g,b,modes,kind,kap,rate,noise,mirror,speed);
 std::ofstream csv(prefix+".csv"),top(prefix+".topology.jsonl");
 csv<<"t,a,H_conf,work_conf,source_conf,Q_com,Q_negative,rho,p,K_phys,G_phys,U_phys,K_phi,G_phi,U_phi,U_cross,U_X,mean_S,min_S,max_S,fplus_04,fminus_04,fplus_05,fminus_05,fplus_06,fminus_06,genuine_plus,genuine_minus,area_core_com_over_V,phi_core_energy_proxy,v2_core_proxy,q_fraction_core,q_fraction_F2gt1,linear_chi_projection,area_smooth_com_over_V,phi_smooth_energy_proxy,v2_smooth_proxy\n";
 scalar(g,b,0,1,rate,0,csv);if(val("--topology","1")=="1")topo_diag(g,b,0,1,top);
 double work=0,source=g.kick(b,0,1,rate),abswork=0;
 for(int it=1;it<=steps;it++){
  double t0=(it-1)*dt,t=it*dt,a0=1+rate*t0,a=1+rate*t;
  g.kick(b,.5*dt,a0,rate);g.drift(b,dt);double src=g.kick(b,.5*dt,a,rate);
  work+=.5*dt*(source+src);abswork+=.5*dt*(std::abs(source)+std::abs(src));source=src;
  if(it%stride==0||it==steps)scalar(g,b,t,a,rate,work,csv);
  if((it%topstride==0||it==steps)&&val("--topology","1")=="1")topo_diag(g,b,t,a,top);
 }
 slice(g,b,1+rate*T,prefix+".slice.f32");
 if(kv.count("--state")){std::ofstream f(kv["--state"],std::ios::binary);for(auto* z:{&b.q,&b.p})for(auto& a:*z)f.write(reinterpret_cast<const char*>(a.data()),sizeof(double)*g.v);}
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::ofstream meta(prefix+".meta.json");meta<<std::setprecision(17)<<"{\"N\":"<<n<<",\"dx\":"<<dx<<",\"dt\":"<<dt<<",\"T\":"<<T<<",\"rate\":"<<rate<<",\"kappa\":"<<kap<<",\"kind\":\""<<kind<<"\",\"steps\":"<<steps<<",\"seconds\":"<<seconds<<",\"abs_work\":"<<abswork<<",\"slice_layout\":\"C-order [N,N,3] little-endian float32, central k=N/2, S,F2,q_phys\"}\n";
 std::cout<<prefix<<" N="<<n<<" steps="<<steps<<" seconds="<<seconds<<" site_steps_per_second="<<g.v*double(steps)/seconds<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
