// Candidate-1 D=0 Cartesian solver. No spatial or angular symmetry reduction.
// Second-order Hamiltonian differences and kick-drift-kick; natural units.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
using V=std::vector<double>;
struct State { std::array<V,3> q,p; State(size_t n){for(auto* z:{&q,&p})for(auto& f:*z) f.assign(n,0.);} };
struct Profile {
 std::vector<std::array<double,5>> d;
 Profile(std::string path){std::ifstream f(path);std::array<double,5> r;while(f>>r[0]>>r[1]>>r[2]>>r[3]>>r[4])d.push_back(r);if(d.size()<2)throw std::runtime_error("profile input");}
 std::array<double,4> at(double r)const{
  if(r>=d.back()[0])return {0,0,1,0};
  auto hi=std::upper_bound(d.begin(),d.end(),r,[](double x,const auto& y){return x<y[0];});
  int k=std::max(0,int(hi-d.begin())-1);auto l=d[k],u=d[k+1];double h=u[0]-l[0],s=(r-l[0])/h;
  std::array<double,4> v{};
  for(int j=0;j<2;j++){int c=1+2*j;v[2*j]=(2*s*s*s-3*s*s+1)*l[c]+(s*s*s-2*s*s+s)*h*l[c+1]+(-2*s*s*s+3*s*s)*u[c]+(s*s*s-s*s)*h*u[c+1];v[2*j+1]=((6*s*s-6*s)*l[c]+(3*s*s-4*s+1)*h*l[c+1]+(-6*s*s+6*s)*u[c]+(3*s*s-2*s)*h*u[c+1])/h;}
  return v;
 }
};
struct Grid {
 int n,nz,s;size_t s2,vol;double h,L,Lz;bool periodic;
 Grid(double hh,double ll,double lz,bool pp):n(int(std::llround(2*ll/hh))),nz(int(std::llround(2*lz/hh))),s(nz+2),s2(size_t(n+2)*s),vol(size_t(n+2)*s2),h(hh),L(ll),Lz(lz),periodic(pp){if(std::abs(n*h-2*L)>1e-8||std::abs(nz*h-2*Lz)>1e-8)throw std::runtime_error("box/h must be integer");if(periodic&&n!=nz)throw std::runtime_error("periodic requires cube");}
 size_t id(int i,int j,int k)const{return size_t(i)*s2+size_t(j)*s+k;}
 double x(int i)const{return -L+(i-.5)*h;}
 double zcoord(int k)const{return -Lz+(k-.5)*h;}
 void halo(State& b,bool tangent=false){
  if(!periodic){if(!tangent){
   for(int j=1;j<=n;j++)for(int k=1;k<=nz;k++){b.q[2][id(0,j,k)]=1;b.q[2][id(n+1,j,k)]=1;}
   for(int i=1;i<=n;i++)for(int k=1;k<=nz;k++){b.q[2][id(i,0,k)]=1;b.q[2][id(i,n+1,k)]=1;}
   for(int i=1;i<=n;i++)for(int j=1;j<=n;j++){b.q[2][id(i,j,0)]=1;b.q[2][id(i,j,nz+1)]=1;}
  }return;}
  for(int c=0;c<3;c++){
   #pragma omp parallel for collapse(2) schedule(static)
   for(int i=1;i<=n;i++)for(int j=1;j<=n;j++){
    b.q[c][id(0,i,j)]=b.q[c][id(n,i,j)];b.q[c][id(n+1,i,j)]=b.q[c][id(1,i,j)];
    b.q[c][id(i,0,j)]=b.q[c][id(i,n,j)];b.q[c][id(i,n+1,j)]=b.q[c][id(i,1,j)];
    b.q[c][id(i,j,0)]=b.q[c][id(i,j,n)];b.q[c][id(i,j,n+1)]=b.q[c][id(i,j,1)];
   }
  }
 }
 double lap(const V& f,size_t z)const{return (f[z-s2]+f[z+s2]+f[z-s]+f[z+s]+f[z-1]+f[z+1]-6*f[z])/(h*h);}
 void kick(State& b,std::vector<State>& tt,double dt){
  halo(b);for(auto& t:tt)halo(t,true);
  #pragma omp parallel for collapse(2) schedule(static)
  for(int i=1;i<=n;i++)for(int j=1;j<=n;j++)for(int k=1;k<=nz;k++){
   size_t z=id(i,j,k);double u=b.q[0][z],w=b.q[1][z],S=b.q[2][z],F2=u*u+w*w,C=S*S+.1*F2;
   b.p[0][z]+=dt*(lap(b.q[0],z)-C*u);b.p[1][z]+=dt*(lap(b.q[1],z)-C*w);b.p[2][z]+=dt*(lap(b.q[2],z)-(S*S-1+F2)*S);
   for(auto& t:tt){double du=t.q[0][z],dw=t.q[1][z],ds=t.q[2][z],df=u*du+w*dw;t.p[0][z]+=dt*(lap(t.q[0],z)-C*du-.2*u*df-2*S*u*ds);t.p[1][z]+=dt*(lap(t.q[1],z)-C*dw-.2*w*df-2*S*w*ds);t.p[2][z]+=dt*(lap(t.q[2],z)-(3*S*S-1+F2)*ds-2*S*df);}
  }
 }
 void drift(State& b,double dt){
  #pragma omp parallel for collapse(2) schedule(static)
  for(int i=1;i<=n;i++)for(int j=1;j<=n;j++)for(int k=1;k<=nz;k++){size_t z=id(i,j,k);for(int c=0;c<3;c++)b.q[c][z]+=dt*b.p[c][z];}
 }
 double gen(const V& f,int i,int j,int k,int type)const{
  size_t z=id(i,j,k);double dx=(f[z+s2]-f[z-s2])/(2*h);
  if(type==0)return -dx;
  return zcoord(k)*dx-x(i)*(f[z+1]-f[z-1])/(2*h);
 }
};
struct Obj {std::array<double,3> center,velocity;double phase,eps;};
void initialize(Grid& g,State& b,const Profile& pr,std::string kind,double eps,double angle,double impact=0){
 std::vector<Obj> objs;std::array<double,3> axis={std::sin(angle),0,std::cos(angle)};
 if(kind=="single")objs.push_back({{.23,-.31,.17},{0,0,0},0,eps});
 else if(kind=="pair")for(int s:{-1,1}){std::array<double,3> c{},v{};for(int j=0;j<3;j++){c[j]=s*12*axis[j];v[j]=-s*.2*axis[j];}c[0]+=s*impact; if(eps!=0){c[1]+=s*.5;v[0]+=s*.003;}objs.push_back({c,v,s==1?1.5707963267948966:0,eps*(s==1?-1:1)});}
 else throw std::runtime_error("case");
 #pragma omp parallel for collapse(2) schedule(static)
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){
  size_t z=g.id(i,j,k);std::array<double,3> r={g.x(i),g.x(j),g.zcoord(k)};b.q[2][z]=1;
  for(auto o:objs){std::array<double,3> rr{},v=o.velocity;double v2=0,dot=0;for(int a=0;a<3;a++){rr[a]=r[a]-o.center[a];v2+=v[a]*v[a];dot+=v[a]*rr[a];}double gamma=1/std::sqrt(1-v2);if(v2>0)for(int a=0;a<3;a++)rr[a]+=(gamma-1)*dot*v[a]/v2;
   const double R[3][3]={{1/std::sqrt(14.),2/std::sqrt(14.),3/std::sqrt(14.)},{2/std::sqrt(5.),-1/std::sqrt(5.),0},{3/std::sqrt(70.),6/std::sqrt(70.),-5/std::sqrt(70.)}};
   double rr2=0,rdot=0,stretch[3]={std::exp(o.eps),std::exp(-.6*o.eps),std::exp(-.4*o.eps)};
   for(int a=0;a<3;a++){double y=0,vy=0;for(int d=0;d<3;d++){y+=R[a][d]*rr[d];vy+=R[a][d]*v[d];}rr2+=stretch[a]*stretch[a]*y*y;rdot+=stretch[a]*stretch[a]*y*vy;}
   double rad=std::sqrt(rr2),rt=rad>0?-gamma*rdot/rad:0;auto f=pr.at(rad);double env=std::exp(-rr2/64.);
   double phase=o.phase+.8000000019763136*gamma*dot+o.eps*env*(std::sin(.7*rr[0]+.2*rr[1]-.3*rr[2])+.7*std::cos(.3*rr[0]-.5*rr[1]+.4*rr[2]));
   double kick=o.eps*env*(std::sin(.4*rr[0]+.7*rr[1])+.5*std::cos(.2*rr[1]+.6*rr[2]));
   double F=f[0],Ft=f[1]*rt+kick*F,co=std::cos(phase),si=std::sin(phase),om=.8000000019763136*gamma;
   b.q[0][z]+=F*co;b.q[1][z]+=F*si;b.q[2][z]+=f[2]-1;b.p[0][z]+=Ft*co+om*F*si;b.p[1][z]+=Ft*si-om*F*co;b.p[2][z]+=f[3]*rt+.5*kick*(1-f[2]);
  }
 }
 g.halo(b);
}
void condensate(Grid& g,State& b,std::string path,double eps){
 std::ifstream in(path);std::vector<std::array<double,6>> modes;std::array<double,6> a;while(in>>a[0]>>a[1]>>a[2]>>a[3]>>a[4]>>a[5])modes.push_back(a);if(modes.empty())throw std::runtime_error("no modes");
 const double pi=std::acos(-1.),om=std::sqrt(.775);
 #pragma omp parallel for collapse(2) schedule(static)
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){
  size_t z=g.id(i,j,k);double nu=0,ns=0,np=0;for(auto m:modes){double phase=pi/g.L*(m[0]*g.x(i)+m[1]*g.x(j)+m[2]*g.zcoord(k));nu+=m[3]*std::cos(phase);ns+=m[4]*std::sin(phase);np+=m[5]*std::cos(phase+.37);}
  double F=.5*(1+eps*nu),phase=eps*np;b.q[0][z]=F*std::cos(phase);b.q[1][z]=F*std::sin(phase);b.q[2][z]=std::sqrt(.75)+eps*ns;b.p[0][z]=om*b.q[1][z];b.p[1][z]=-om*b.q[0][z];
 }
 // Fix the spatial mean initial charge to the declared .25*sqrt(.775).
 double sum=0;
 #pragma omp parallel for collapse(2) reduction(+:sum) schedule(static)
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){size_t z=g.id(i,j,k);sum+=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z];}
 double scale=std::sqrt(.25*double(g.n)*g.n*g.nz/sum);
 #pragma omp parallel for collapse(2) schedule(static)
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){size_t z=g.id(i,j,k);for(int c=0;c<2;c++){b.q[c][z]*=scale;b.p[c][z]*=scale;}}
 g.halo(b);
}
void final_data(Grid& g,State& b,std::string path){
 // Exact six-neighbour 3D components. Threshold-selected charge is not an
 // asymptotic soliton mass function, especially in a periodic production box.
 double peak=0;for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){size_t z=g.id(i,j,k);peak=std::max(peak,b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z]);}
 std::ofstream out(path+".components.csv");out<<"threshold,rank,cells,Q,x,y,z\n";double dv=g.h*g.h*g.h;
 for(double frac:{.02,.05,.10}){
  std::vector<unsigned char> seen(g.vol);std::vector<size_t> queue;std::vector<std::array<double,5>> cs;
  for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){
   size_t start=g.id(i,j,k);if(seen[start])continue;seen[start]=1;if(b.q[0][start]*b.q[0][start]+b.q[1][start]*b.q[1][start]<frac*peak)continue;queue.clear();queue.push_back(start);double Q=0,m=0,mx=0,my=0,mz=0;
   for(size_t a=0;a<queue.size();a++){
    size_t z=queue[a];int x=int(z/g.s2),y=int((z%g.s2)/g.s),zz=int(z%g.s);double q=b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z],qp=std::max(0.,q);Q+=q;m+=qp;mx+=qp*g.x(x);my+=qp*g.x(y);mz+=qp*g.zcoord(zz);
    for(int d=0;d<6;d++){int xn=x,yn=y,zn=zz;if(d==0)xn--;if(d==1)xn++;if(d==2)yn--;if(d==3)yn++;if(d==4)zn--;if(d==5)zn++;if(g.periodic){if(xn<1)xn=g.n;if(xn>g.n)xn=1;if(yn<1)yn=g.n;if(yn>g.n)yn=1;if(zn<1)zn=g.nz;if(zn>g.nz)zn=1;}else if(xn<1||xn>g.n||yn<1||yn>g.n||zn<1||zn>g.nz)continue;size_t nb=g.id(xn,yn,zn);if(seen[nb])continue;seen[nb]=1;if(b.q[0][nb]*b.q[0][nb]+b.q[1][nb]*b.q[1][nb]>=frac*peak)queue.push_back(nb);}
   }cs.push_back({double(queue.size()),Q*dv,m>0?mx/m:0,m>0?my/m:0,m>0?mz/m:0});
  }
  std::sort(cs.begin(),cs.end(),[](const auto& x,const auto& y){return x[1]>y[1];});int rank=0;for(auto c:cs){out<<std::setprecision(16)<<frac<<','<<++rank;for(double x:c)out<<','<<x;out<<'\n';}
 }
 int stride=std::max(1,(std::max(g.n,g.nz)+63)/64),nn=(g.n+stride-1)/stride,nnz=(g.nz+stride-1)/stride;std::ofstream bin(path+".spatial.f32",std::ios::binary);std::ofstream meta(path+".spatial.json");meta<<"{\"N\":"<<g.n<<",\"Nz\":"<<g.nz<<",\"n\":"<<nn<<",\"nz\":"<<nnz<<",\"stride\":"<<stride<<",\"h\":"<<std::setprecision(17)<<g.h<<",\"L\":"<<g.L<<",\"Lz\":"<<g.Lz<<",\"layout\":\"C-order [i,j,k,channel], float32 native little endian; channels F2,S,q\"}\n";
 for(int i=1;i<=g.n;i+=stride)for(int j=1;j<=g.n;j+=stride)for(int k=1;k<=g.nz;k+=stride){size_t z=g.id(i,j,k);float f[3]={float(b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z]),float(b.q[2][z]),float(b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z])};bin.write(reinterpret_cast<const char*>(f),sizeof(f));}
}
void diag(Grid& g,State& b,std::vector<State>& tt,std::ofstream& out,double t,std::string sym){
 // Maxima in opposing half-spaces give seed centers; iterate signed-charge barycenters
 // inside R=8. Report R=8 and R=12 integrals, not density thresholds alone.
 struct Core{std::array<double,3> c{};double peak=-1;};std::array<Core,2> cores;
 double E=0,Q=0,neg=0,maxF=0,minS=1;
 for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){size_t z=g.id(i,j,k);double F2=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z];int side=g.zcoord(k)>0?1:0;if(F2>cores[side].peak){cores[side].peak=F2;cores[side].c={g.x(i),g.x(j),g.zcoord(k)};}}
 for(auto& core:cores)for(int it=0;it<4;it++){double mass=0;std::array<double,3> cm{};for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){std::array<double,3> r={g.x(i),g.x(j),g.zcoord(k)};double d2=0;for(int c=0;c<3;c++)d2+=(r[c]-core.c[c])*(r[c]-core.c[c]);if(d2>64)continue;size_t z=g.id(i,j,k);double q=std::max(0.,b.q[1][z]*b.p[0][z]-b.q[0][z]*b.p[1][z]);mass+=q;for(int c=0;c<3;c++)cm[c]+=q*r[c];}if(mass>0)for(int c=0;c<3;c++)core.c[c]=cm[c]/mass;}
 std::array<double,2> q8{},q12{},e8{},e12{},size2{};std::array<std::array<double,6>,2> cov{};
 std::vector<double> tn(tt.size()),gn(tt.size()),rn(tt.size());std::array<double,9> gram{};std::array<double,3> dot{};
 #pragma omp parallel
 {
  double le=0,lq=0,ln=0,lmax=0,lmin=1;std::array<double,2> lq8{},lq12{},le8{},le12{},lsize{};std::array<std::array<double,6>,2> lcov{};std::vector<double> ltn(tt.size()),lgn(tt.size()),lrn(tt.size());std::array<double,9> lg{};std::array<double,3> ld{};
  #pragma omp for collapse(2) schedule(static)
  for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){
   size_t z=g.id(i,j,k);double u=b.q[0][z],w=b.q[1][z],S=b.q[2][z],F2=u*u+w*w,q=w*b.p[0][z]-u*b.p[1][z],e=.25*(S*S-1)*(S*S-1)+.5*S*S*F2+.025*F2*F2;
   for(int c=0;c<3;c++){e+=.5*b.p[c][z]*b.p[c][z];for(size_t d:{g.s2,size_t(g.s),size_t(1)}){double dq=b.q[c][z+d]-b.q[c][z];e+=.5*dq*dq/(g.h*g.h);}if(!g.periodic){if(i==1){double dq=b.q[c][z]-b.q[c][z-g.s2];e+=.5*dq*dq/(g.h*g.h);}if(j==1){double dq=b.q[c][z]-b.q[c][z-g.s];e+=.5*dq*dq/(g.h*g.h);}if(k==1){double dq=b.q[c][z]-b.q[c][z-1];e+=.5*dq*dq/(g.h*g.h);}}}
   le+=e;lq+=q;ln+=std::max(0.,-q);lmax=std::max(lmax,F2);lmin=std::min(lmin,S);
   std::array<double,3> r={g.x(i),g.x(j),g.zcoord(k)};
   for(int a=0;a<2;a++){std::array<double,3> d{};double d2=0;for(int c=0;c<3;c++){d[c]=r[c]-cores[a].c[c];d2+=d[c]*d[c];}if(d2<144){lq12[a]+=q;le12[a]+=e;}if(d2<64){lq8[a]+=q;le8[a]+=e;lsize[a]+=std::max(0.,q)*d2;int c=0;for(int x=0;x<3;x++)for(int y=x;y<3;y++)lcov[a][c++]+=std::max(0.,q)*d[x]*d[y];}}
   if(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]<1296)for(int c=0;c<3;c++)for(int p=0;p<2;p++){
    const auto& f=p?b.p[c]:b.q[c];double phase=c==0?-(p?b.p[1][z]:b.q[1][z]):c==1?(p?b.p[0][z]:b.q[0][z]):0;std::array<double,3> G={g.gen(f,i,j,k,0),g.gen(f,i,j,k,1),phase};
    for(size_t a=0;a<tt.size();a++){double gen=sym=="u1"?phase:a<2?G[a]:0;double v=p?tt[a].p[c][z]:tt[a].q[c][z];ltn[a]+=v*v;lgn[a]+=gen*gen;lrn[a]+=(v-gen)*(v-gen);}
    if(tt.size()==3){double v=p?tt[2].p[c][z]:tt[2].q[c][z];for(int a=0;a<3;a++){ld[a]+=v*G[a];for(int b=0;b<3;b++)lg[3*a+b]+=G[a]*G[b];}}
   }
  }
  #pragma omp critical
  {E+=le;Q+=lq;neg+=ln;maxF=std::max(maxF,lmax);minS=std::min(minS,lmin);for(int a=0;a<2;a++){q8[a]+=lq8[a];q12[a]+=lq12[a];e8[a]+=le8[a];e12[a]+=le12[a];size2[a]+=lsize[a];for(int c=0;c<6;c++)cov[a][c]+=lcov[a][c];}for(size_t a=0;a<tt.size();a++){tn[a]+=ltn[a];gn[a]+=lgn[a];rn[a]+=lrn[a];}for(int a=0;a<9;a++)gram[a]+=lg[a];for(int a=0;a<3;a++)dot[a]+=ld[a];}
 }
 double dv=g.h*g.h*g.h;out<<std::setprecision(16)<<t<<','<<E*dv<<','<<Q*dv<<','<<neg*dv<<','<<maxF<<','<<minS;
 for(int a=0;a<2;a++){for(double x:cores[a].c)out<<','<<x;out<<','<<q8[a]*dv<<','<<e8[a]*dv<<','<<q12[a]*dv<<','<<e12[a]*dv<<','<<size2[a]*dv;for(double x:cov[a])out<<','<<x*dv;}
 for(size_t a=0;a<tt.size();a++)out<<','<<std::sqrt(tn[a]*dv)<<','<<std::sqrt(gn[a]*dv)<<','<<std::sqrt(rn[a]*dv);
 if(tt.size()==3){for(double x:gram)out<<','<<x*dv;for(double x:dot)out<<','<<x*dv;}
 out<<'\n';out.flush();
}
int main(int argc,char** argv){try{
 std::map<std::string,std::string> m;for(int i=1;i<argc;i+=2){if(i+1>=argc)throw std::runtime_error("key value args");m[argv[i]]=argv[i+1];}
 auto val=[&](std::string k,std::string d){return m.count(k)?m[k]:d;};double h=std::stod(val("--h","0.5")),L=std::stod(val("--L","36")),T=std::stod(val("--T","120")),dt=std::stod(val("--dt","0.05")),eps=std::stod(val("--eps","0")),angle=std::stod(val("--angle","0"));int threads=std::stoi(val("--threads","8"));omp_set_num_threads(threads);
 std::string kind=val("--case","single"),sym=val("--sym","none");Grid g(h,L,std::stod(val("--Lz",val("--L","36"))),kind=="condensate");State b(g.vol);Profile pr(val("--profile","profile.dat"));if(kind=="condensate")condensate(g,b,val("--modes","modes.dat"),eps);else initialize(g,b,pr,kind,eps,angle,std::stod(val("--impact","0")));
 std::vector<State> tt;int nt=sym=="all"?3:sym=="space"?2:(sym=="u1"||sym=="fd")?1:0;tt.reserve(nt);for(int a=0;a<nt;a++){
  tt.emplace_back(g.vol);auto& t=tt.back();
  if(a==2){if(kind!="pair")throw std::runtime_error("relative shift requires pair");State plus(g.vol),minus(g.vol);double delta=std::stod(val("--seed-delta","0.0001"));initialize(g,plus,pr,kind,eps,angle,delta);initialize(g,minus,pr,kind,eps,angle,-delta);for(int c=0;c<3;c++)for(size_t z=0;z<g.vol;z++){t.q[c][z]=(plus.q[c][z]-minus.q[c][z])/(2*delta);t.p[c][z]=(plus.p[c][z]-minus.p[c][z])/(2*delta);}continue;}
  #pragma omp parallel for collapse(2) schedule(static)
  for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){size_t z=g.id(i,j,k);for(int c=0;c<3;c++)for(int p=0;p<2;p++){auto& f=p?b.p[c]:b.q[c];double v=sym=="fd"?(p?(.1+.15*c)*f[z]:(.3-.1*c)*(f[z]-(c==2?1:0))):sym=="u1"?(c==0?-(p?b.p[1][z]:b.q[1][z]):c==1?(p?b.p[0][z]:b.q[0][z]):0):g.gen(f,i,j,k,a);(p?t.p[c]:t.q[c])[z]=v;}}
 }
 if(sym=="fd"){double offset=std::stod(val("--fd-offset","0"));for(int c=0;c<3;c++)for(size_t z=0;z<g.vol;z++){b.q[c][z]+=offset*tt[0].q[c][z];b.p[c][z]+=offset*tt[0].p[c][z];}}
 std::string path=val("--out","history.csv");std::ofstream out(path);out<<"t,E,Q,Qnegative,maxF2,minS";for(int a=0;a<2;a++){std::string p=",c"+std::to_string(a)+"_";out<<p<<"x"<<p<<"y"<<p<<"z"<<p<<"Q8"<<p<<"E8"<<p<<"Q12"<<p<<"E12"<<p<<"Qr2";for(auto c:{"xx","xy","xz","yy","yz","zz"})out<<p<<c;}for(int a=0;a<nt;a++)out<<",tan"<<a<<"_norm,gen"<<a<<"_norm,res"<<a<<"_norm";if(nt==3){for(int a=0;a<9;a++)out<<",gram"<<a;for(int a=0;a<3;a++)out<<",dot"<<a;}out<<'\n';
 auto begin=std::chrono::steady_clock::now();diag(g,b,tt,out,0,sym);int steps=int(std::llround(T/dt)),stride=std::max(1,int(std::llround(std::stod(val("--sample","10"))/dt)));if(std::abs(steps*dt-T)>1e-8)throw std::runtime_error("T/dt integer");
 g.kick(b,tt,.5*dt);
 for(int n=1;n<=steps;n++){
  g.drift(b,dt);for(auto& t:tt)g.drift(t,dt);
  if(n%stride==0||n==steps){g.kick(b,tt,.5*dt);diag(g,b,tt,out,n*dt,sym);if(n<steps)g.kick(b,tt,.5*dt);}
  else g.kick(b,tt,dt);
 }
 final_data(g,b,path);
 if(m.count("--state-out")){std::ofstream bin(m["--state-out"],std::ios::binary);for(int a=-1;a<nt;a++){const State& f=a<0?b:tt[a];for(int p=0;p<2;p++)for(int c=0;c<3;c++)for(int i=1;i<=g.n;i++)for(int j=1;j<=g.n;j++)for(int k=1;k<=g.nz;k++){double v=(p?f.p[c]:f.q[c])[g.id(i,j,k)];bin.write(reinterpret_cast<const char*>(&v),sizeof(v));}}}
 double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();std::cout<<"N="<<g.n<<" Nz="<<g.nz<<" steps="<<steps<<" seconds="<<sec<<" site_steps_per_second="<<double(g.n)*g.n*g.nz*steps/sec<<"\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
