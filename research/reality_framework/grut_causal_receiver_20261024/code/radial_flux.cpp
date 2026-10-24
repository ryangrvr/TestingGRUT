// Causal CFS massive-wave layer. Physical u,v retain the quartic action.
// psi lives on faces, q on nodes; their causal initial memory is zero.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
using Vec=std::vector<double>;
constexpr double PI=3.14159265358979323846;
struct State{std::array<Vec,12>a;double er=0,qr=0,ef=0,qf=0;explicit State(int n){for(auto&v:a)v.assign(n+1,0.);}};
int probe_cut=0;
int n;double h,ab,width,smax,alpha,linear_mass;bool nonlinear;Vec sn,sf;
double node(const Vec&u,int i){if(i<0)return -u[-i];if(i>n)return -u[2*n-i];return u[i];}
double face(const Vec&u,int i){if(i<0)return u[-1-i];if(i>=n)return u[2*n-1-i];return u[i];}
double grad(const Vec&u,int i){return ((9./8)*(node(u,i+1)-node(u,i))-(1./24)*(node(u,i+2)-node(u,i-1)))/h;}
double divg(const Vec&u,int i){return ((9./8)*(face(u,i)-face(u,i-1))-(1./24)*(face(u,i+1)-face(u,i-2)))/h;}
double lap(const Vec&u,int i){return ((-365./144)*u[i]+(87./64)*(node(u,i-1)+node(u,i+1))
    -(3./32)*(node(u,i-2)+node(u,i+2))+(1./576)*(node(u,i-3)+node(u,i+3)))/(h*h);}
void rhs(const State&u,State&d){
 d.er=0;d.qr=0;d.ef=0;d.qf=0;
 for(auto&v:d.a)std::fill(v.begin(),v.end(),0.);
 for(int i=1;i<n;i++){
  double r=i*h,xx=u.a[0][i],yy=u.a[1][i],zz=u.a[2][i],density=(xx*xx+yy*yy)/(r*r),S=1+zz/r;
  std::array<double,3>force;
  if(nonlinear)force={(S*S+.1*density)*xx,(S*S+.1*density)*yy,r*(3*(S*S-1)+density)*S};
  else force={linear_mass*linear_mass*xx,linear_mass*linear_mass*yy,6*zz};
  double cx=0,cy=0;
  for(int f=0;f<3;f++){
   double lu=lap(u.a[f],i),extra=divg(u.a[6+f],i),corr=extra+u.a[9+f][i];
   d.a[f][i]=u.a[3+f][i];d.a[3+f][i]=lu-force[f]+corr;
   d.a[9+f][i]=-(alpha+sn[i])*u.a[9+f][i]-sn[i]*(lu+extra);
   d.er-=4*PI*h*u.a[3+f][i]*corr;
   if(f==0)cx=corr;if(f==1)cy=corr;
  }
  d.qr+=4*PI*h*(xx*cy-yy*cx);
 }
 for(int i=0;i<n;i++)for(int f=0;f<3;f++){
  d.a[6+f][i]=-(alpha+sf[i])*u.a[6+f][i]-sf[i]*grad(u.a[f],i);
 }
 // Exact currents for the symmetric seven-point Hamiltonian, across a node cut.
 const std::array<double,3>coupling={87./64,-3./32,1./576};
 for(int i=std::max(1,probe_cut-2);i<=probe_cut;i++)for(int ell=1;ell<=3;ell++){
  int j=i+ell;if(j<=probe_cut||j>=n)continue;
  double coefficient=4*PI/h*coupling[ell-1];
  for(int f=0;f<3;f++)d.ef+=.5*coefficient*(u.a[f][i]*u.a[3+f][j]-u.a[3+f][i]*u.a[f][j]);
  d.qf+=coefficient*(u.a[0][i]*u.a[1][j]-u.a[1][i]*u.a[0][j]);
 }
}
void stage(State&out,const State&u,const State&d,double dt){
 for(int k=0;k<12;k++)for(int i=0;i<=n;i++)out.a[k][i]=u.a[k][i]+dt*d.a[k][i];
}
void step(State&u,double dt,State&v,State&k1,State&k2,State&k3,State&k4){
 rhs(u,k1);stage(v,u,k1,dt/2);rhs(v,k2);stage(v,u,k2,dt/2);rhs(v,k3);stage(v,u,k3,dt);rhs(v,k4);
 for(int k=0;k<12;k++)for(int i=0;i<=n;i++)u.a[k][i]+=dt/6*(k1.a[k][i]+2*k2.a[k][i]+2*k3.a[k][i]+k4.a[k][i]);
 u.er+=dt/6*(k1.er+2*k2.er+2*k3.er+k4.er);u.qr+=dt/6*(k1.qr+2*k2.qr+2*k3.qr+k4.qr);
 u.ef+=dt/6*(k1.ef+2*k2.ef+2*k3.ef+k4.ef);u.qf+=dt/6*(k1.qf+2*k2.qf+2*k3.qf+k4.qf);
}
struct Budget{double E=0,Q=0,core=0,massdelta=0,Fmax=0,Eprobe=0,Qprobe=0;};
Budget budget(const State&u){Budget b;double weight=4*PI*h;
 for(int i=1;i<n;i++){
  double r=i*h,xx=u.a[0][i],yy=u.a[1][i],zz=u.a[2][i],dens=(xx*xx+yy*yy)/(r*r),S=1+zz/r;
  double pot=nonlinear?r*r*(.75*(S*S-1)*(S*S-1)+.5*S*S*dens+.025*dens*dens)
                      :.5*linear_mass*linear_mass*(xx*xx+yy*yy)+3*zz*zz;
  double e=pot;for(int f=0;f<3;f++)e+=.5*u.a[3+f][i]*u.a[3+f][i]-.5*u.a[f][i]*lap(u.a[f],i);
  double q=yy*u.a[3][i]-xx*u.a[4][i];b.E+=weight*e;b.Q+=weight*q;if(i<=probe_cut){b.Eprobe+=weight*e;b.Qprobe+=weight*q;}
  if(r<60-h*.1)b.core+=weight*q;else if(std::abs(r-60)<h*.1)b.core+=.5*weight*q;
  if(r>=ab)b.massdelta=std::max(b.massdelta,std::abs(S*S-1+.1*dens));
  b.Fmax=std::max(b.Fmax,std::sqrt(dens));
 }return b;
}
int main(int argc,char**argv){try{
 if(argc!=13)throw std::runtime_error("usage INITIAL PREFIX h dt L tmax probe ab sigma alpha nonlinear linear_mass");
 std::string prefix=argv[2];h=std::stod(argv[3]);double dt=std::stod(argv[4]),L=std::stod(argv[5]),end=std::stod(argv[6]),probe=std::stod(argv[7]);
 ab=std::stod(argv[8]);width=L-ab;smax=std::stod(argv[9]);alpha=std::stod(argv[10]);nonlinear=std::stoi(argv[11]);linear_mass=std::stod(argv[12]);
 n=std::lround(L/h);sn.assign(n+1,0.);sf.assign(n+1,0.);
 for(int i=0;i<=n;i++){double z=std::clamp((i*h-ab)/width,0.,1.);sn[i]=smax*z*z*z;
  z=std::clamp(((i+.5)*h-ab)/width,0.,1.);sf[i]=smax*z*z*z;}
 State u(n),tmp(n),k1(n),k2(n),k3(n),k4(n);std::ifstream in(argv[1]);std::string line;std::getline(in,line);
 for(int i=1;i<n;i++){
  if(!std::getline(in,line))throw std::runtime_error("short initial grid");std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
  double r,x,y,S,px,py,ps;row>>r>>x>>y>>S>>px>>py>>ps;if(std::abs(r-i*h)>1e-8)throw std::runtime_error("grid mismatch");
  u.a[0][i]=r*x;u.a[1][i]=r*y;u.a[2][i]=r*(S-1);u.a[3][i]=r*px;u.a[4][i]=r*py;u.a[5][i]=r*ps;
 }
 std::ofstream hist(prefix+"_history.csv"),trace(prefix+"_trace.csv");hist<<std::setprecision(17);trace<<std::setprecision(17);
 hist<<"tau,E,Q,E_receiver,Q_receiver,Q60,max_layer_mass2_shift,Fmax,E_inside_cut,Q_inside_cut,E_flux_cut,Q_flux_cut\n";
 trace<<"tau,x,y,s_minus1,dx,dy,ds,px,py,ps\n";
 int ip=std::clamp((int)std::lround(probe/h),2,n-2),steps=std::lround(end/dt),hs=std::max(1,(int)std::lround(10/dt)),ts=std::max(1,(int)std::lround(.25/dt));
 probe_cut=ip;
 auto b0=budget(u);
 for(int j=0;j<=steps;j++){
  if(j%hs==0||j==steps){auto b=budget(u);if(!std::isfinite(b.E)||b.Fmax>100)throw std::runtime_error("runaway");
   hist<<j*dt<<','<<b.E<<','<<b.Q<<','<<u.er<<','<<u.qr<<','<<b.core<<','<<b.massdelta<<','<<b.Fmax<<','<<b.Eprobe<<','<<b.Qprobe<<','<<u.ef<<','<<u.qf<<'\n';
   if(j%std::max(hs,steps/12/hs*hs)==0)std::cout<<std::setprecision(8)<<"PROGRESS "<<j*dt<<" Eerr "<<(b.E+u.er)/b0.E-1<<" Qerr "<<(b.Q+u.qr)/b0.Q-1<<'\n'<<std::flush;
  }
  if(j%ts==0||j==steps){double r=ip*h;trace<<j*dt;
   for(int f=0;f<3;f++)trace<<','<<u.a[f][ip]/r;
   for(int f=0;f<3;f++)trace<<','<<(node(u.a[f],ip+1)/(r+h)-node(u.a[f],ip-1)/(r-h))/(2*h);
   for(int f=0;f<3;f++)trace<<','<<u.a[3+f][ip]/r;trace<<'\n';
  }
  if(j<steps)step(u,dt,tmp,k1,k2,k3,k4);
 }
 std::ofstream final(prefix+"_final.csv"),memory(prefix+"_memory.csv");final<<std::setprecision(17)<<"rho,x,y,s,px,py,ps\n";
 memory<<std::setprecision(17)<<"rho_node,psi_x_at_next_face,psi_y_at_next_face,psi_s_at_next_face,q_x,q_y,q_s\n";
 for(int i=1;i<n;i++){double r=i*h;final<<r<<','<<u.a[0][i]/r<<','<<u.a[1][i]/r<<','<<1+u.a[2][i]/r;
  for(int f=0;f<3;f++)final<<','<<u.a[3+f][i]/r;final<<'\n';
  memory<<r;for(int f=6;f<12;f++)memory<<','<<u.a[f][i];memory<<'\n';}
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
