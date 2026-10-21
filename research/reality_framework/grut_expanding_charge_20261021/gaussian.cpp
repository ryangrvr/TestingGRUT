// Exploratory preparation wrapper around the unchanged validated integrator.
#define main frozen_entrypoint
#include "expanding.cpp"
#undef main
int main(int argc,char** argv){try{
 std::map<std::string,std::string> kv;for(int i=1;i<argc;i+=2){if(i+1==argc)throw std::runtime_error("key/value arguments required");kv[argv[i]]=argv[i+1];}
 auto val=[&](std::string k,std::string d){return kv.count(k)?kv[k]:d;};
 int N=std::stoi(val("--N","32")),threads=std::stoi(val("--threads","8")),mirror=std::stoi(val("--mirror","1"));
 double dx=std::stod(val("--dx","0.4")),dt=std::stod(val("--dt","0.06")),T=std::stod(val("--T","48")),rate=.025,kap=2.;
 omp_set_num_threads(threads);Grid g(N,dx);State b(g.v);auto start=std::chrono::steady_clock::now();
 std::ifstream input(kv.at("--noise-file"),std::ios::binary);
 for(int c=0;c<3;c++)input.read(reinterpret_cast<char*>(b.q[c].data()),sizeof(double)*g.v);
 if(!input)throw std::runtime_error("initial-noise length");
 double F=std::pow(kap,1./3),omega=std::sqrt(.1)*F,mean=0;
 for(size_t z=0;z<g.v;z++){b.q[0][z]=F*(1+.001*b.q[0][z]);b.q[1][z]=F*.001*b.q[1][z];b.q[2][z]=mirror*.01*b.q[2][z];mean+=b.q[0][z]*b.q[0][z]+b.q[1][z]*b.q[1][z];}
 double norm=std::sqrt(F*F*g.v/mean);
 for(size_t z=0;z<g.v;z++){b.q[0][z]*=norm;b.q[1][z]*=norm;b.p[0][z]=omega*b.q[1][z];b.p[1][z]=-omega*b.q[0][z];b.p[2][z]=rate*b.q[2][z];}
 std::string prefix=kv.at("--out");
 // A raw state can be supplied to check exact evolution identity with expanding.cpp.
 if(kv.count("--state-in")){std::ifstream in(kv["--state-in"],std::ios::binary);for(auto* z:{&b.q,&b.p})for(auto& f:*z)in.read(reinterpret_cast<char*>(f.data()),sizeof(double)*g.v);if(!in)throw std::runtime_error("initial-state length");}
 std::ofstream csv(prefix+".csv"),top(prefix+".topology.jsonl");
 csv<<"t,a,H_conf,work_conf,source_conf,Q_com,Q_negative,rho,p,K_phys,G_phys,U_phys,K_phi,G_phi,U_phi,U_cross,U_X,mean_S,min_S,max_S,fplus_04,fminus_04,fplus_05,fminus_05,fplus_06,fminus_06,genuine_plus,genuine_minus,area_core_com_over_V,phi_core_energy_proxy,v2_core_proxy,q_fraction_core,q_fraction_F2gt1,linear_chi_projection,area_smooth_com_over_V,phi_smooth_energy_proxy,v2_smooth_proxy\n";
 scalar(g,b,0,1,rate,0,csv);bool topology_on=val("--topology","1")=="1";if(topology_on)topo_diag(g,b,0,1,top);
 int steps=std::llround(T/dt),stride=std::llround(1.2/dt),topstride=std::llround(4.8/dt);if(std::abs(steps*dt-T)>1e-8)throw std::runtime_error("T/dt");
 double source=g.kick(b,0,1,rate),work=0,abswork=0;
 for(int it=1;it<=steps;it++){
  double t0=(it-1)*dt,t=it*dt,a0=1+rate*t0,a=1+rate*t;
  g.kick(b,.5*dt,a0,rate);g.drift(b,dt);double src=g.kick(b,.5*dt,a,rate);
  work+=.5*dt*(source+src);abswork+=.5*dt*(std::abs(source)+std::abs(src));source=src;
  if(it%stride==0||it==steps)scalar(g,b,t,a,rate,work,csv);
  if((it%topstride==0||it==steps)&&topology_on)topo_diag(g,b,t,a,top);
 }
 slice(g,b,1+rate*T,prefix+".slice.f32");
 if(kv.count("--state")){std::ofstream f(kv["--state"],std::ios::binary);for(auto* z:{&b.q,&b.p})for(auto& a:*z)f.write(reinterpret_cast<const char*>(a.data()),sizeof(double)*g.v);}
 double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::ofstream meta(prefix+".meta.json");meta<<std::setprecision(17)<<"{\"N\":"<<N<<",\"dx\":"<<dx<<",\"dt\":"<<dt<<",\"T\":"<<T<<",\"rate\":"<<rate<<",\"kappa\":2,\"kind\":\"gaussian_charge\",\"steps\":"<<steps<<",\"seconds\":"<<seconds<<",\"abs_work\":"<<abswork<<",\"slice_layout\":\"C-order [N,N,3] little-endian float32, central k=N/2, S,F2,q_phys\"}\n";
 std::cout<<prefix<<" N="<<N<<" steps="<<steps<<" seconds="<<seconds<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
