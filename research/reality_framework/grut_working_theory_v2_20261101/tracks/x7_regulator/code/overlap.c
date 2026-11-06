/* Continue X5's archived covariance without reinitializing its orbitals. */
#define main x5_original_main
#include "../../x5_coupled/code/coupled.c"
#undef main

static void load_checkpoint(const char *name, const Params *p, Grid *g) {
    FILE *f=fopen(name,"rb");
    if (!f) fail("Cannot open checkpoint");
    size_t bins=p->nq, modes=2*bins*p->ns;
    if (fread(g->q,sizeof *g->q,bins,f)!=bins || fread(g->w,sizeof *g->w,bins,f)!=bins
        || fread(g->mode,sizeof *g->mode,modes,f)!=modes || fgetc(f)!=EOF)
        fail("Checkpoint length does not match model");
    fclose(f);
}

static void project(FILE *f, const Params *p, const Grid *g, const State *s) {
    double gas_p[NS]={0}, a3=cube(s->a);
    for (int j=0; j<p->ns; j++) {
        const Species sp=p->species[j];
        for (int i=0; i<2*p->nq; i++) {
            int bin=i%p->nq;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]/s->a;
            double E=sqrt(k*k+sp.mass*sp.mass);
            Mode r=g->mode[j*2*p->nq+i];
            double occupation=0.5*(1+(sp.mass*r.x+k*r.z)/E);
            gas_p[j]+=g->w[bin]*k*k*occupation/(3*E);
        }
        fprintf(f,"%.17g ",gas_p[j]/a3);
    }
    Stress x=stress(p,g,s,1);
    Matching m=matching(p,s->a);
    double scalar_force[NF], fermion_force[NF]={0};
    potential_gradient(p,s->phi,s->phi,scalar_force);
    for (int j=0; j<p->ns; j++)
        fermion_force[p->species[j].parent]+=x.force[j];
    double work=0, pressure_ct=0;
    for (int j=0; j<NF; j++) {
        double A=m.inertia[j], slope=m.derivative[j], v=x.velocity[j];
        double acceleration=-(scalar_force[j]+fermion_force[j]+x.hubble*(3*A+slope)*v)/A;
        double matched=fermion_force[j]-(1-A)*acceleration-x.hubble*(3*(1-A)-slope)*v;
        work+=matched*v;
        pressure_ct+=(0.5*(1-A)-slope/6)*v*v;
    }
    fprintf(f,"%.17g %.17g\n",work/x.rho,pressure_ct/x.rho);
}

int main(int argc, char **argv) {
    if (argc!=5) fail("Usage: overlap INPUT CHECKPOINT HISTORY PROJECTION");
    _Static_assert(sizeof(Mode)==3*sizeof(double),"Mode checkpoint requires packed doubles");
    Params p={0}; State s={0}; Grid g={0};
    read_input(argv[1],&p,&s,&g);
    load_checkpoint(argv[2],&p,&g);
    FILE *out=fopen(argv[3],"w"), *projection=fopen(argv[4],"w");
    if (!out || !projection) fail("Cannot open overlap output");
    fprintf(out,"# t a H rho pressure work residual sigma chi s dsigma dchi ds rho_s rho_sm rho_dr"
                " scalar_ct gravity_ct max_iterations norm_max f_min f_max Q_sigma Q_chi Q_s");
    for (int j=1; j<=NS; j++)
        fprintf(out," rho_%d n_%d p_%d dressed_rho_%d dressed_n_%d J_%d",j,j,j,j,j,j);
    fputc('\n',out);

    const double dt=(p.t1-p.t0)/p.nsteps, E0=stress(&p,&g,&s,0).rho*cube(s.a);
    double work=0;
    int iterations=0;
    record(out,NULL,&p,&g,&s,p.t0,work,E0,iterations);
    project(projection,&p,&g,&s);
    for (int step=1; step<=p.nsteps; step++) {
        expansion(&p,&g,&s,dt/2,&work);
        decays(&p,&g,&s,dt/2);
        int n=matter(&p,&g,&s,dt);
        if (n>iterations) iterations=n;
        decays(&p,&g,&s,dt/2);
        expansion(&p,&g,&s,dt/2,&work);
        if (step%p.stride==0 || step==p.nsteps) {
            record(out,NULL,&p,&g,&s,p.t0+step*dt,work,E0,iterations);
            project(projection,&p,&g,&s);
        }
    }
    fclose(out); fclose(projection);
    free(g.q); free(g.w); free(g.mode); free(g.trial);
    return 0;
}
