// GRUT-M1 local scalar dynamics. Derived from reviewed conservative pilot.
// Includes the one-vacuum cubic deformation; no expansion/heat/quantum kinetics.
// Conservative cold-field pilot. Canonical fields q=(sqrt2 ReX,sqrt2 ImX,phi).
// Parameters and coordinates use m_X units. Periodic, fourth-order spatial
// Hamiltonian; velocity Verlet or fourth-order symmetric composition.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
using V=std::vector<double>;
int main(int argc,char**argv){
 if(argc!=9&&argc!=10){std::cerr<<"init prefix dt tend sample_dt snapshot_dt order threads [snapshot_bits]\n";return 2;}
 std::ifstream in(argv[1],std::ios::binary);char magic[8];int n;double L,par[5];
 in.read(magic,8);in.read((char*)&n,4);in.read((char*)&L,8);in.read((char*)par,40);
 if(std::string(magic,8)!="GRUTM1D\n"||!in||n<4)return 3;
 const size_t ns=(size_t)n*n*n;double dx=L/n,idx2=1/(dx*dx),nu=par[0],g2=par[1],lx=par[2],lp=par[3],eps=par[4];
 std::array<V,3> q,p;for(auto&f:q){f.resize(ns);in.read((char*)f.data(),8*ns);}for(auto&f:p){f.resize(ns);in.read((char*)f.data(),8*ns);}if(!in)return 4;
 double dt=std::atof(argv[3]),end=std::atof(argv[4]),sample=std::atof(argv[5]),snap=std::atof(argv[6]);int order=std::atoi(argv[7]);omp_set_num_threads(std::atoi(argv[8]));
 long steps=std::llround(end/dt);dt=end/steps;long every=std::max<long>(1,std::llround(sample/dt)),sv=std::max<long>(1,std::llround(snap/dt));
 int bits=argc==10?std::atoi(argv[9]):32;if(bits!=32&&bits!=64)return 5;
 std::vector<int> plus(n),minus(n),p2(n),m2(n);for(int i=0;i<n;i++){plus[i]=(i+1)%n;minus[i]=(i+n-1)%n;p2[i]=(i+2)%n;m2[i]=(i+n-2)%n;}
 auto lap=[&](const V&f,size_t at,int x,int y,int z){
  size_t xx=(size_t)x*n*n,yy=(size_t)y*n;
  double a=f[(size_t)plus[x]*n*n+yy+z]+f[(size_t)minus[x]*n*n+yy+z]+f[xx+(size_t)plus[y]*n+z]+f[xx+(size_t)minus[y]*n+z]+f[xx+yy+plus[z]]+f[xx+yy+minus[z]];
  double b=f[(size_t)p2[x]*n*n+yy+z]+f[(size_t)m2[x]*n*n+yy+z]+f[xx+(size_t)p2[y]*n+z]+f[xx+(size_t)m2[y]*n+z]+f[xx+yy+p2[z]]+f[xx+yy+m2[z]];
  return idx2*((4./3)*a-b/12.-7.5*f[at]);
 };
 auto kick=[&](double h){
  #pragma omp parallel for schedule(static)
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int z=0;z<n;z++){
   size_t i=((size_t)x*n+y)*n+z;double a=q[0][i],b=q[1][i],s=q[2][i],r2=a*a+b*b;double f=g2*s*s+lx*r2;
   p[0][i]+=h*(lap(q[0],i,x,y,z)-f*a);p[1][i]+=h*(lap(q[1],i,x,y,z)-f*b);p[2][i]+=h*(lap(q[2],i,x,y,z)-(lp*(s*s-nu*nu)+g2*r2)*s+3*eps*nu*(s-nu)*(s-nu));
  }
 };
 auto drift=[&](double h){
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<ns;i++)for(int j=0;j<3;j++)q[j][i]+=h*p[j][i];
 };
 auto step=[&](double h){kick(h/2);drift(h);kick(h/2);};
 std::ofstream out(std::string(argv[2])+".tsv");out<<std::setprecision(17)<<"t\trho\tnQ\tkinetic\tgradient\tpotential\tphi_mean\tphi_rms\tphi_negative\tlow_phi\tcharge_rms\tcharge_min\tcharge_max\tseconds\n";
 auto start=std::chrono::steady_clock::now();
 auto diagnostics=[&](double t){double ke=0,ge=0,ve=0,Q=0,Q2=0,ph=0,ph2=0;long neg=0,low=0;double qmin=1e300,qmax=-1e300;
  #pragma omp parallel for reduction(+:ke,ge,ve,Q,Q2,ph,ph2,neg,low) reduction(min:qmin) reduction(max:qmax) schedule(static)
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int z=0;z<n;z++){
   size_t i=((size_t)x*n+y)*n+z;double a=q[0][i],b=q[1][i],s=q[2][i],r2=a*a+b*b;double j=a*p[1][i]-b*p[0][i];
   for(int c=0;c<3;c++){
    ke+=.5*p[c][i]*p[c][i];const V&f=q[c];size_t xx=(size_t)x*n*n,yy=(size_t)y*n;
    double d1x=f[(size_t)plus[x]*n*n+yy+z]-f[i],d1y=f[xx+(size_t)plus[y]*n+z]-f[i],d1z=f[xx+yy+plus[z]]-f[i];
    double d2x=f[(size_t)p2[x]*n*n+yy+z]-f[i],d2y=f[xx+(size_t)p2[y]*n+z]-f[i],d2z=f[xx+yy+p2[z]]-f[i];
    ge+=.5*idx2*((4./3)*(d1x*d1x+d1y*d1y+d1z*d1z)-(d2x*d2x+d2y*d2y+d2z*d2z)/12.);
   }ve+=lp/4*(s*s-nu*nu)*(s*s-nu*nu)-eps*nu*(s-nu)*(s-nu)*(s-nu)+g2/2*s*s*r2+lx/4*r2*r2;
   Q+=j;Q2+=j*j;ph+=s;ph2+=s*s;neg+=s<0;low+=std::abs(s)<.5*nu;qmin=std::min(qmin,j);qmax=std::max(qmax,j);
  }
  double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  out<<t<<'\t'<<(ke+ge+ve)/ns<<'\t'<<Q/ns<<'\t'<<ke/ns<<'\t'<<ge/ns<<'\t'<<ve/ns<<'\t'<<ph/ns<<'\t'<<std::sqrt(ph2/ns)<<'\t'<<(double)neg/ns<<'\t'<<(double)low/ns<<'\t'<<std::sqrt(Q2/ns)<<'\t'<<qmin<<'\t'<<qmax<<'\t'<<wall<<'\n';out.flush();
  if(!std::isfinite(ke+ge+ve)||std::abs(sqrt(Q2/ns))>1e12)throw std::runtime_error("Nonfinite/divergent state");
 };
 auto snapshot=[&](long s){char tail[64];std::snprintf(tail,64,"_s%07ld.bin",s);std::ofstream f(std::string(argv[2])+tail,std::ios::binary);double t=s*dt;f.write((char*)&n,4);f.write((char*)&L,8);f.write((char*)&t,8);f.write((char*)&bits,4);std::vector<float> buf(bits==32?ns:0);
  for(auto&array:q){if(bits==64)f.write((char*)array.data(),8*ns);else{for(size_t i=0;i<ns;i++)buf[i]=array[i];f.write((char*)buf.data(),4*ns);}}for(auto&array:p){if(bits==64)f.write((char*)array.data(),8*ns);else{for(size_t i=0;i<ns;i++)buf[i]=array[i];f.write((char*)buf.data(),4*ns);}}
 };
 diagnostics(0);snapshot(0);double w1=1/(2-std::cbrt(2.)),w0=1-2*w1;
 for(long s=1;s<=steps;s++){
  if(order==4){kick(w1*dt/2);drift(w1*dt);kick((w1+w0)*dt/2);drift(w0*dt);kick((w0+w1)*dt/2);drift(w1*dt);kick(w1*dt/2);}else step(dt);
  if(s%every==0||s==steps)diagnostics(s*dt);
  if(s%sv==0||s==steps)snapshot(s);
 }
 std::cout<<"done n="<<n<<" steps="<<steps<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 return 0;
}
