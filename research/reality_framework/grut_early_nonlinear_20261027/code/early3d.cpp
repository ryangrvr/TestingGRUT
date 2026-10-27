// Early conformal quartic scalar-limit pilot, adapted from verified8159556.
// Canonical chi=(sqrt2 Re(aX),sqrt2 Im(aX),a phi), H_i units.
// Fixed radiation clock; nonlinear X/phi, no heavy/gauge/quantum kinetics.
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
 if(argc!=13){std::cerr<<"init prefix hmax eta_end sample_tau snapshot_tau order threads bits tau0 a0 C\n";return 2;}
 std::ifstream in(argv[1],std::ios::binary);char magic[8];int n;double L,par[4];
 in.read(magic,8);in.read((char*)&n,4);in.read((char*)&L,8);in.read((char*)par,32);
 if(std::string(magic,8)!="GRUT3D1\n"||!in||n<4)return 3;
 const size_t ns=(size_t)n*n*n;double dx=L/n,idx2=1/(dx*dx),nu=par[0],g2=par[1],lx=par[2],lp=par[3];
 std::array<V,3> q,p;for(auto&f:q){f.resize(ns);in.read((char*)f.data(),8*ns);}for(auto&f:p){f.resize(ns);in.read((char*)f.data(),8*ns);}if(!in)return 4;
 double dt=std::atof(argv[3]),end=std::atof(argv[4]),sample=std::atof(argv[5]),snap=std::atof(argv[6]);int order=std::atoi(argv[7]);omp_set_num_threads(std::atoi(argv[8]));
 double tau0=std::atof(argv[10]),a0=std::atof(argv[11]),C=std::atof(argv[12]);
 auto clock=[&](double s){return tau0+a0*s+C*s*s/2;};
 auto eta_of_tau=[&](double t){return 2*(t-tau0)/(a0+std::sqrt(a0*a0+2*C*(t-tau0)));};
 int bits=std::atoi(argv[9]);if(bits!=32&&bits!=64)return 5;
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
   p[0][i]+=h*(lap(q[0],i,x,y,z)-f*a);p[1][i]+=h*(lap(q[1],i,x,y,z)-f*b);p[2][i]+=h*(lap(q[2],i,x,y,z)-(lp*(s*s-nu*nu)+g2*r2)*s);
  }
 };
 auto drift=[&](double h){
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<ns;i++)for(int j=0;j<3;j++)q[j][i]+=h*p[j][i];
 };
 auto step=[&](double h){kick(h/2);drift(h);kick(h/2);};
 std::ofstream out(std::string(argv[2])+".tsv");
 out<<std::setprecision(17)<<"eta\ttau\ta\trho\tnQ\tKX\tKphi\tGX\tGphi\tUX\tUphi\tUmix\tq0\tq1\tqphi\tp0\tp1\tpphi\tX2\tphi2\tcharge_RMS\tmax_X2\twork_phi_to_X\tseconds\n";
 auto start=std::chrono::steady_clock::now();
 auto diagnostics=[&](double s){
  double kx=0,kp=0,gx=0,gp=0,ux=0,up=0,um=0,Q=0,Q2=0,q0=0,q1=0,qp=0,p0=0,p1=0,pp=0,rmsx=0,rmsp=0,maxx=0,work=0;
  #pragma omp parallel for reduction(+:kx,kp,gx,gp,ux,up,um,Q,Q2,q0,q1,qp,p0,p1,pp,rmsx,rmsp,work) reduction(max:maxx) schedule(static)
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int z=0;z<n;z++){
   size_t i=((size_t)x*n+y)*n+z;double a=q[0][i],b=q[1][i],f=q[2][i],r2=a*a+b*b,j=a*p[1][i]-b*p[0][i];
   kx+=.5*(p[0][i]*p[0][i]+p[1][i]*p[1][i]);kp+=.5*p[2][i]*p[2][i];
   for(int c=0;c<3;c++){
    const V&v=q[c];size_t xx=(size_t)x*n*n,yy=(size_t)y*n;
    double d1x=v[(size_t)plus[x]*n*n+yy+z]-v[i],d1y=v[xx+(size_t)plus[y]*n+z]-v[i],d1z=v[xx+yy+plus[z]]-v[i];
    double d2x=v[(size_t)p2[x]*n*n+yy+z]-v[i],d2y=v[xx+(size_t)p2[y]*n+z]-v[i],d2z=v[xx+yy+p2[z]]-v[i];
    double ge=.5*idx2*((4./3)*(d1x*d1x+d1y*d1y+d1z*d1z)-(d2x*d2x+d2y*d2y+d2z*d2z)/12.);
    if(c<2)gx+=ge;else gp+=ge;
   }
   ux+=lx*r2*r2/4;up+=lp*f*f*f*f/4;um+=g2*f*f*r2/2;
   Q+=j;Q2+=j*j;q0+=a;q1+=b;qp+=f;p0+=p[0][i];p1+=p[1][i];pp+=p[2][i];rmsx+=r2;rmsp+=f*f;maxx=std::max(maxx,r2);work+=g2*f*p[2][i]*r2;
  }
  double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  out<<s<<'\t'<<clock(s)<<'\t'<<a0+C*s<<'\t'<<(kx+kp+gx+gp+ux+up+um)/ns<<'\t'<<Q/ns<<'\t'<<kx/ns<<'\t'<<kp/ns<<'\t'<<gx/ns<<'\t'<<gp/ns<<'\t'<<ux/ns<<'\t'<<up/ns<<'\t'<<um/ns<<'\t'<<q0/ns<<'\t'<<q1/ns<<'\t'<<qp/ns<<'\t'<<p0/ns<<'\t'<<p1/ns<<'\t'<<pp/ns<<'\t'<<rmsx/ns<<'\t'<<rmsp/ns<<'\t'<<std::sqrt(Q2/ns)<<'\t'<<maxx<<'\t'<<work/ns<<'\t'<<wall<<'\n';out.flush();
  if(!std::isfinite(kx+kp+gx+gp+ux+up+um)||maxx>1e12)throw std::runtime_error("nonfinite/divergent state");
 };
 auto snapshot=[&](long index,double elapsed){char tail[64];std::snprintf(tail,64,"_s%07ld.bin",index);std::ofstream f(std::string(argv[2])+tail,std::ios::binary);double t=elapsed;f.write((char*)&n,4);f.write((char*)&L,8);f.write((char*)&t,8);f.write((char*)&bits,4);std::vector<float> buf(bits==32?ns:0);
  for(auto&array:q){if(bits==64)f.write((char*)array.data(),8*ns);else{for(size_t i=0;i<ns;i++)buf[i]=array[i];f.write((char*)buf.data(),4*ns);}}for(auto&array:p){if(bits==64)f.write((char*)array.data(),8*ns);else{for(size_t i=0;i<ns;i++)buf[i]=array[i];f.write((char*)buf.data(),4*ns);}}
 };
 diagnostics(0);snapshot(0,0);double w1=1/(2-std::cbrt(2.)),w0=1-2*w1;
 double s=0,next_tau=tau0+sample,next_snap=tau0+snap;long steps=0,snapshot_index=0;
 while(s<end-1e-13){
  double target=std::min(end,eta_of_tau(next_tau));
  while(s<target-1e-13){double h=std::min(dt,target-s);
   if(order==4){kick(w1*h/2);drift(w1*h);kick((w1+w0)*h/2);drift(w0*h);kick((w0+w1)*h/2);drift(w1*h);kick(w1*h/2);}else step(h);
   s+=h;steps++;
  }
  s=target;diagnostics(s);next_tau+=sample;
  if(clock(s)>=next_snap-1e-8||s>=end-1e-13){snapshot(++snapshot_index,s);while(next_snap<=clock(s)+1e-8)next_snap+=snap;}
 }
 std::cout<<"done n="<<n<<" steps="<<steps<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 return 0;
}
