// Reviewed FD4/Yoshida scalar Hamiltonian exposed to Python for streaming FFTs.
// A tangent is the exact derivative of the same discrete symplectic map.
#include <cmath>
#include <vector>
#include <omp.h>
extern "C" int evolve(double*z,double*d,int nt,int n,double L,double h,int steps,int threads){
 if(n<4||L<=0||steps<0||nt<0)return 1;
 omp_set_num_threads(threads);
 const size_t ns=(size_t)n*n*n;
 const double idx2=(n/L)*(n/L),lx=.001,g2=.01,lp=.03;
 std::vector<int> plus(n),minus(n),p2(n),m2(n);
 for(int i=0;i<n;i++){plus[i]=(i+1)%n;minus[i]=(i+n-1)%n;p2[i]=(i+2)%n;m2[i]=(i+n-2)%n;}
 auto lap=[&](const double*f,size_t at,int x,int y,int zz){
  size_t xx=(size_t)x*n*n,yy=(size_t)y*n;
  double a=f[(size_t)plus[x]*n*n+yy+zz]+f[(size_t)minus[x]*n*n+yy+zz]+f[xx+(size_t)plus[y]*n+zz]+f[xx+(size_t)minus[y]*n+zz]+f[xx+yy+plus[zz]]+f[xx+yy+minus[zz]];
  double b=f[(size_t)p2[x]*n*n+yy+zz]+f[(size_t)m2[x]*n*n+yy+zz]+f[xx+(size_t)p2[y]*n+zz]+f[xx+(size_t)m2[y]*n+zz]+f[xx+yy+p2[zz]]+f[xx+yy+m2[zz]];
  return idx2*((4./3)*a-b/12.-7.5*f[at]);
 };
 auto kick=[&](double dt){
  #pragma omp parallel for schedule(static)
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int zz=0;zz<n;zz++){
   size_t i=((size_t)x*n+y)*n+zz;
   double a=z[i],b=z[ns+i],s=z[2*ns+i],r2=a*a+b*b,f=g2*s*s+lx*r2;
   z[3*ns+i]+=dt*(lap(z,i,x,y,zz)-f*a);
   z[4*ns+i]+=dt*(lap(z+ns,i,x,y,zz)-f*b);
   z[5*ns+i]+=dt*(lap(z+2*ns,i,x,y,zz)-(lp*s*s+g2*r2)*s);
   for(int t=0;t<nt;t++){
    double*v=d+(size_t)t*6*ns;
    double da=v[i],db=v[ns+i],ds=v[2*ns+i],dot=a*da+b*db;
    v[3*ns+i]+=dt*(lap(v,i,x,y,zz)-f*da-2*lx*a*dot-2*g2*a*s*ds);
    v[4*ns+i]+=dt*(lap(v+ns,i,x,y,zz)-f*db-2*lx*b*dot-2*g2*b*s*ds);
    v[5*ns+i]+=dt*(lap(v+2*ns,i,x,y,zz)-(3*lp*s*s+g2*r2)*ds-2*g2*s*dot);
   }
  }
 };
 auto drift=[&](double dt){
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<3*ns;i++){
   z[i]+=dt*z[3*ns+i];
   for(int t=0;t<nt;t++){double*v=d+(size_t)t*6*ns;v[i]+=dt*v[3*ns+i];}
  }
 };
 const double w1=1/(2-std::cbrt(2.)),w0=1-2*w1;
 for(int j=0;j<steps;j++){
  kick(w1*h/2);drift(w1*h);kick((w1+w0)*h/2);drift(w0*h);kick((w0+w1)*h/2);drift(w1*h);kick(w1*h/2);
 }
 return 0;
}
