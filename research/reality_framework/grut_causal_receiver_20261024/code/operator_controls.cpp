#define main production_main
#include "radial_cfs.cpp"
#undef main
#include <string>

void setup(int count,double spacing,double mass){
 n=count;h=spacing;ab=n*h/2;width=n*h-ab;smax=20;alpha=.1;linear_mass=mass;
 sn.assign(n+1,0);sf.assign(n+1,0);
 for(int i=0;i<=n;i++){
  sn[i]=smax*std::pow(std::clamp((i*h-ab)/width,0.,1.),3);
  sf[i]=smax*std::pow(std::clamp(((i+.5)*h-ab)/width,0.,1.),3);
 }
}
int main(int argc,char**argv){
 std::cout<<std::setprecision(17);
 if(argc>1&&std::string(argv[1])=="matrix"){
  setup(std::stoi(argv[2]),std::stod(argv[3]),std::stod(argv[4]));nonlinear=false;
  std::vector<std::pair<int,int>>indices;
  for(int k:{0,3})for(int i=1;i<n;i++)indices.push_back({k,i});
  for(int i=n/2;i<n;i++)indices.push_back({6,i});
  for(int i=n/2+1;i<n;i++)indices.push_back({9,i});
  int size=indices.size();std::vector<Vec>matrix(size,Vec(size));
  for(int j=0;j<size;j++){
   State u(n),d(n);u.a[indices[j].first][indices[j].second]=1;rhs(u,d);
   for(int i=0;i<size;i++)matrix[i][j]=d.a[indices[i].first][indices[i].second];
  }
  for(auto&row:matrix){for(int j=0;j<size;j++){if(j)std::cout<<',';std::cout<<row[j];}std::cout<<'\n';}
  return 0;
 }
 setup(48,.4,1);nonlinear=true;State u(n),d(n),plus(n),minus(n);
 for(int i=1;i<n;i++){
  double r=i*h,env=std::exp(-r*r/70);
  u.a[0][i]=r*.12*env*std::cos(.2*r);u.a[1][i]=r*.1*env*std::sin(.2*r);
  u.a[2][i]=r*.012*env;u.a[3][i]=r*.03*env;u.a[4][i]=-r*.09*env;u.a[5][i]=r*.008*env;
  if(sn[i]>0)for(int f=0;f<3;f++)u.a[9+f][i]=.003*std::cos(r+f);
 }
 for(int i=0;i<n;i++)if(sf[i]>0)for(int f=0;f<3;f++)u.a[6+f][i]=.002*std::sin((i+.5)*h+f);
 rhs(u,d);double eps=1e-6;stage(plus,u,d,eps);stage(minus,u,d,-eps);
 auto b=budget(u),bp=budget(plus),bm=budget(minus);double quadratic=0,gradient=0;
 for(int f=0;f<3;f++){
  for(int i=1;i<n;i++)quadratic-=u.a[f][i]*lap(u.a[f],i);
  for(int i=0;i<n;i++)gradient+=std::pow(grad(u.a[f],i),2);
 }
 std::cout<<"{\"gradient_energy_identity_relative_error\":"<<std::abs(gradient/quadratic-1)
  <<",\"instantaneous_energy_plus_receiver_residual_over_E\":"<<std::abs((bp.E-bm.E)/(2*eps)+d.er)/b.E
  <<",\"instantaneous_charge_plus_receiver_residual_over_Q\":"<<std::abs((bp.Q-bm.Q)/(2*eps)+d.qr)/std::abs(b.Q)<<"}\n";
}
