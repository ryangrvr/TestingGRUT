/* Diagonal vacuum Born kernel in a locally constant mass-phase rotation. */
#define main x5_original_main
#include "../../x5_coupled/code/coupled.c"
#undef main

typedef struct {
    double omega, nx, nz;
    double plus, minus, upper_energy, lower_energy;
} Rates;

typedef struct {
    double energy_error, charge_error;
    int bisected;
} Audit;

static inline Rates rates(Species sp,double k,double delta) {
    double z=k+delta, omega=hypot(sp.mass,z), p=fabs(k);
    Rates out={0};
    if (omega==0) return out;
    double scale=sp.width*sp.mass/(2*omega);
    out.omega=omega; out.nx=sp.mass/omega; out.nz=z/omega;
    out.upper_energy=omega-delta; out.lower_energy=omega+delta;
    out.plus=(out.upper_energy>p) ? scale : 0;
    out.minus=(out.lower_energy>p) ? scale : 0;
    return out;
}

static inline Mode damp(Mode r,Rates R,double dt,double *lost) {
    double b=R.nx*r.x+R.nz*r.z;
    double loss=-expm1(-(R.plus+R.minus)*dt), t=sqrt(1-loss);
    double next=b-loss*(1+b);
    *lost=loss*(1+b)/2;
    return (Mode){t*(r.x-b*R.nx)+next*R.nx,t*r.y,t*(r.z-b*R.nz)+next*R.nz};
}

static void rotating_step(const Params *p,Grid *g,State *s,double dt,Audit *audit,int depth) {
    Matching match=matching(p,s->a);
    double a3=cube(s->a), inv[NF], v[NF];
    for (int j=0;j<NF;j++) {
        inv[j]=1/(a3*match.inertia[j]);
        v[j]=s->mom[j]*inv[j];
    }
    for (int it=0;it<16;it++) {
        double moment[NF]={0}, sterile[NF]={0}, daughters[NF]={0};
        double radiation=0, fermion_energy=0, norm=0;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:moment[:NF],sterile[:NF],daughters[:NF],radiation,fermion_energy,norm)
        for (int j=0;j<p->ns;j++) {
            for (int i=0;i<2*p->nq;i++) {
                Species sp=p->species[j];
                int bin=i%p->nq, id=j*2*p->nq+i;
                double k=(i<p->nq ? -1 : 1)*g->q[bin]/s->a, w=g->w[bin];
                double delta=v[sp.parent]/(2*sp.f);
                Rates R=rates(sp,k,delta);
                Mode old=g->mode[id];
                g->trial[id]=old;
                double total=R.plus+R.minus;
                if (total==0) continue;
                double lost;
                Mode next=damp(old,R,dt,&lost);
                g->trial[id]=next;
                double b=(R.plus-R.minus)/total;
                double ds=w*(next.z-old.z)/(4*sp.f);
                double dr=w*lost*b/(2*sp.f);
                double add=w*lost*(R.omega-b*delta);
                double de=w*(sp.mass*(next.x-old.x)+k*(next.z-old.z))/2;
                moment[sp.parent]+=ds+dr;
                sterile[sp.parent]+=ds;
                daughters[sp.parent]+=dr;
                radiation+=add; fermion_energy+=de;
                norm+=fabs(add)+fabs(de)+fabs(v[sp.parent]*(ds+dr));
            }
        }
        double error=0, next_v[NF];
        for (int j=0;j<NF;j++) {
            next_v[j]=(s->mom[j]+0.5*moment[j])*inv[j];
            error=fmax(error,fabs(next_v[j]-v[j])/(1+fabs(next_v[j])+fabs(v[j])));
        }
        if (error<2e-13) {
            double energy=fermion_energy+radiation;
            for (int j=0;j<NF;j++) {
                energy+=moment[j]*next_v[j];
                double charge=moment[j]-sterile[j]-daughters[j];
                double scale=fabs(moment[j])+fabs(sterile[j])+fabs(daughters[j])+1;
                audit->charge_error=fmax(audit->charge_error,fabs(charge)/scale);
                s->mom[j]+=moment[j];
            }
            audit->energy_error=fmax(audit->energy_error,fabs(energy)/(1+norm));
            s->rad_sm+=s->a*radiation;
            Mode *swap=g->mode;g->mode=g->trial;g->trial=swap;
            return;
        }
        memcpy(v,next_v,sizeof v);
    }
    if (depth==12) fail("Rotating collision midpoint did not converge");
    audit->bisected++;
    rotating_step(p,g,s,dt/2,audit,depth+1);
    rotating_step(p,g,s,dt/2,audit,depth+1);
}

static void project(FILE *out,const Params *p,const Grid *g,const State *s) {
    Stress x=stress(p,g,s,1);
    Matching m=matching(p,s->a);
    double a3=cube(s->a), gradient[NF], force[NF]={0}, acceleration[NF];
    potential_gradient(p,s->phi,s->phi,gradient);
    for (int j=0;j<p->ns;j++) force[p->species[j].parent]+=x.force[j];
    double work=0, pressure_ct=0;
    for (int j=0;j<NF;j++) {
        double A=m.inertia[j], slope=m.derivative[j], v=x.velocity[j];
        acceleration[j]=-(gradient[j]+force[j]+x.hubble*(3*A+slope)*v)/A;
        double matched=force[j]-(1-A)*acceleration[j]-x.hubble*(3*(1-A)-slope)*v;
        work+=matched*v;
        pressure_ct+=(0.5*(1-A)-slope/6)*v*v;
    }
    double source=0, charge=0, safe_eigen=0, safe_threshold=0, safe_both=0;
    double second_open=0, dressed_density=0, slowest=0;
    for (int j=0;j<p->ns;j++) {
        Species sp=p->species[j];
        double delta=x.velocity[sp.parent]/(2*sp.f), ddelta=acceleration[sp.parent]/(2*sp.f);
        double gas_p=0;
        for (int i=0;i<2*p->nq;i++) {
            int bin=i%p->nq;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]/s->a, w=g->w[bin]/a3;
            double E=hypot(k,sp.mass);
            Mode r=g->mode[j*2*p->nq+i];
            double f=0.5*(1+(sp.mass*r.x+k*r.z)/E);
            gas_p+=w*k*k*f/(3*E);
            Rates R=rates(sp,k,delta);
            double dressed=0.5*(1+R.nx*r.x+R.nz*r.z);
            double total=R.plus+R.minus;
            double zdot=ddelta-x.hubble*k, omegadot=R.nz*zdot;
            double eigen=sp.mass*fabs(zdot)/(2*R.omega*R.omega*R.omega);
            double secular=total/(2*R.omega);
            dressed_density+=w*R.omega*dressed;
            slowest=fmax(slowest,total);
            for (int b=-1;b<=1;b+=2) {
                double rate=b>0 ? R.plus : R.minus;
                if (rate==0) continue;
                double energy=b>0 ? R.upper_energy : R.lower_energy;
                double power=w*rate*dressed*energy;
                double gap=energy-fabs(k);
                double gapdot=omegadot-b*ddelta+x.hubble*fabs(k);
                double threshold=fmax(fabs(gapdot)/(gap*gap),fmax(x.hubble,total)/gap);
                source+=power;
                charge+=b*w*rate*dressed;
                if (R.plus>0 && R.minus>0) second_open+=power;
                if (eigen<=0.1 && secular<=0.1) safe_eigen+=power;
                if (threshold<=0.1) safe_threshold+=power;
                if (eigen<=0.1 && secular<=0.1 && threshold<=0.1) safe_both+=power;
            }
        }
        fprintf(out,"%.17g ",gas_p);
    }
    fprintf(out,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            work/x.rho,pressure_ct/x.rho,source,charge,safe_eigen,safe_threshold,
            safe_both,second_open,dressed_density,slowest);
}

static int controls(void) {
    Species sp={0,2,3,0.04};
    const double deltas[]={-2.,-0.3,0.,0.3,2.};
    const double momenta[]={0.,0.2,0.8,4.,20.};
    for (size_t j=0;j<sizeof deltas/sizeof *deltas;j++) {
        for (size_t i=0;i<sizeof momenta/sizeof *momenta;i++) {
            for (int h=-1;h<=1;h+=2) {
                double k=h*momenta[i], d=deltas[j], lost;
                Rates R=rates(sp,k,d);
                Mode r={0.2,-0.3,0.6}, n=damp(r,R,2.,&lost);
                printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                       k,d,R.omega,R.plus,R.minus,R.upper_energy,R.lower_energy,n.x,n.y,n.z,lost);
            }
        }
    }
    return 0;
}

int main(int argc,char **argv) {
    if (argc==2 && strcmp(argv[1],"--controls")==0) return controls();
    if (argc!=6) fail("Usage: rotating INPUT CHECKPOINT HISTORY PROJECTION AUDIT");
    _Static_assert(sizeof(Mode)==3*sizeof(double),"Checkpoint Mode must be packed");
    Params p={0}; State s={0}; Grid g={0}; Audit audit={0};
    read_input(argv[1],&p,&s,&g);
    size_t n=p.nq, states=2*n*p.ns;
    FILE *checkpoint=fopen(argv[2],"rb");
    if (!checkpoint || fread(g.q,sizeof *g.q,n,checkpoint)!=n ||
        fread(g.w,sizeof *g.w,n,checkpoint)!=n || fread(g.mode,sizeof *g.mode,states,checkpoint)!=states ||
        fgetc(checkpoint)!=EOF) fail("Invalid archived covariance");
    fclose(checkpoint);

    FILE *out=fopen(argv[3],"w"), *projection=fopen(argv[4],"w"), *ledger=fopen(argv[5],"w");
    if (!out || !projection || !ledger) fail("Cannot open replay outputs");
    fprintf(out,"# t a H rho pressure work residual sigma chi s dsigma dchi ds rho_s rho_sm rho_dr"
                " scalar_ct gravity_ct max_iterations norm_max f_min f_max Q_sigma Q_chi Q_s");
    for (int j=1;j<=NS;j++) fprintf(out," rho_%d n_%d p_%d dressed_rho_%d dressed_n_%d J_%d",j,j,j,j,j,j);
    fputc('\n',out);

    double dt=(p.t1-p.t0)/p.nsteps,work=0,E0=stress(&p,&g,&s,0).rho*cube(s.a);
    int max_it=0;
    record(out,NULL,&p,&g,&s,p.t0,work,E0,0);
    project(projection,&p,&g,&s);
    for (int step=1;step<=p.nsteps;step++) {
        expansion(&p,&g,&s,dt/2,&work);
        for (int half=0;half<2;half++) {
            double loss=-expm1(-p.gamma_s*dt);
            s.rad_dark+=0.5*s.mom[2]*s.mom[2]*loss/(s.a*s.a);
            s.mom[2]*=exp(-p.gamma_s*dt/2);
            rotating_step(&p,&g,&s,dt/2,&audit,0);
            if (!half) {
                int it=matter(&p,&g,&s,dt);
                if (it>max_it) max_it=it;
            }
        }
        expansion(&p,&g,&s,dt/2,&work);
        if (step%p.stride==0 || step==p.nsteps) {
            record(out,NULL,&p,&g,&s,p.t0+step*dt,work,E0,max_it);
            project(projection,&p,&g,&s);
        }
    }
    fprintf(ledger,"%.17g %.17g %d\n",audit.energy_error,audit.charge_error,audit.bisected);
    fclose(out);fclose(projection);fclose(ledger);
    free(g.q);free(g.w);free(g.mode);free(g.trial);
    return 0;
}
