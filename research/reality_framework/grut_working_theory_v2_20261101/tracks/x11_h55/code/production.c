/* Shared-FRW scalar lattice with site-local coherent Majorana continua. */
#define main x5_main
#include "../../x5_coupled/code/coupled.c"
#undef main

typedef struct {
    int n, cells;
    double inv_dx2, a;
    int *links;
    State *fields;
    Grid *grids;
} Box;

typedef struct {
    double rho, pressure, hubble, gradient, scalar_ct, dark, visible, quantum_e, quantum_p;
} Average;

static void gradient(const Box *box, int i, double a, Matching match,
                     double *energy, double *pressure, double lap[NF]) {
    for (int j=0; j<NF; j++) {
        double center=box->fields[i].phi[j], sum=0, edges=0;
        for (int l=0; l<6; l++) {
            double d=box->fields[box->links[6*i+l]].phi[j]-center;
            sum+=d;
            if (l%2==0) edges+=d*d;
        }
        double g=0.5*box->inv_dx2*edges/(a*a);
        *energy+=match.inertia[j]*g;
        *pressure-=(match.inertia[j]+match.derivative[j])*g/3;
        lap[j]=box->inv_dx2*sum;
    }
}

static void kick(const Params *p, Box *box, double dt) {
    double a=box->a;
    Matching match=matching(p,a);
    #pragma omp parallel for schedule(static)
    for (int i=0; i<box->cells; i++) {
        double lap[NF], g=0, pressure=0;
        gradient(box,i,a,match,&g,&pressure,lap);
        for (int j=0; j<NF; j++)
            box->fields[i].mom[j]+=dt*a*match.inertia[j]*lap[j];
    }
}

/* Hold canonical momenta and covariances fixed in an expansion substep. */
static Average average(const Params *p, const Box *box, double a) {
    Matching match=matching(p,a);
    double rho=0, pressure=0, gradients=0, scalar_ct=0, dark=0, visible=0, quantum_e=0, quantum_p=0;
    double a4=cube(a)*a;
    #pragma omp parallel for schedule(static) reduction(+:rho,pressure,gradients,scalar_ct,dark,visible,quantum_e,quantum_p)
    for (int i=0; i<box->cells; i++) {
        State s=box->fields[i];
        s.a=a;
        Stress local=stress(p,&box->grids[i],&s,0);
        double g=0, pg=0, lap[NF];
        gradient(box,i,a,match,&g,&pg,lap);
        rho+=local.rho+g;
        pressure+=local.pressure+pg;
        gradients+=g;
        scalar_ct+=local.scalar_ct;
        dark+=s.rad_dark/a4;
        visible+=s.rad_sm/a4;
        for (int j=0; j<p->ns; j++) {
            quantum_e+=local.e[j];
            quantum_p+=local.p[j];
        }
    }
    double mean=1.0/box->cells;
    rho*=mean;
    return (Average){rho,pressure*mean,p->expand ? sqrt(rho/(3*p->mp*p->mp+match.gravity/2)) : 0,
                     gradients*mean,scalar_ct*mean,dark*mean,visible*mean,quantum_e*mean,quantum_p*mean};
}

static void expand_box(const Params *p, Box *box, double dt, double *work, double *qwork) {
    if (!p->expand) return;
    double a=box->a;
    Average start=average(p,box,a);
    Average middle=average(p,box,a*exp(0.5*dt*start.hubble));
    box->a=a*exp(dt*middle.hubble);
    double dv=cube(a)*expm1(3*dt*middle.hubble);
    *work+=middle.pressure*dv;
    *qwork+=middle.quantum_p*dv;
}

static int exchange(const Params *p, Box *box, double dt, double *work) {
    double a3=cube(box->a), power=0;
    Matching match=matching(p,box->a);
    int iterations=0;
    #pragma omp parallel for schedule(static) reduction(max:iterations) reduction(+:power)
    for (int i=0; i<box->cells; i++) {
        State *s=&box->fields[i];
        s->a=box->a;
        decays(p,&box->grids[i],s,dt/2);
        State old=*s;
        int steps=matter(p,&box->grids[i],s,dt);
        if (steps>iterations) iterations=steps;
        if (p->ns) {
            double force[NF];
            potential_gradient(p,old.phi,s->phi,force);
            for (int j=0; j<NF-1; j++) {
                double impulse=s->mom[j]-old.mom[j]+a3*dt*force[j];
                double midpoint=0.5*(old.mom[j]+s->mom[j])/(a3*match.inertia[j]);
                power-=midpoint*impulse;
            }
        }
        decays(p,&box->grids[i],s,dt/2);
    }
    *work+=power/box->cells;
    return iterations;
}

static int index_at(int n, int x, int y, int z) {
    return ((x+n)%n*n+(y+n)%n)*n+(z+n)%n;
}

static Box prepare_box(const Params *p, State base, Grid grid, const char *file) {
    FILE *fp=fopen(file,"r");
    if (!fp) fail("Cannot read scalar seed");
    int n=read_int(fp);
    double L=read_real(fp), H0=read_real(fp);
    if (n<1 || n>64 || L<=0) fail("Invalid spatial grid");
    int cells=n*n*n;
    Box box={n,cells,pow(n/L,2),base.a,malloc(6*cells*sizeof(int)),
             malloc(cells*sizeof(State)),malloc(cells*sizeof(Grid))};
    if (!box.links || !box.fields || !box.grids) fail("Cannot allocate cells");
    Matching match=matching(p,base.a);
    size_t modes=(size_t)(p->ns ? p->ns : 1)*2*p->nq;
    for (int x=0; x<n; x++) for (int y=0; y<n; y++) for (int z=0; z<n; z++) {
        int i=index_at(n,x,y,z);
        int neighbors[6]={index_at(n,x+1,y,z),index_at(n,x-1,y,z),
                          index_at(n,x,y+1,z),index_at(n,x,y-1,z),
                          index_at(n,x,y,z+1),index_at(n,x,y,z-1)};
        memcpy(&box.links[6*i],neighbors,sizeof neighbors);
        State *s=&box.fields[i];
        *s=base;
        for (int j=0; j<NF; j++) s->phi[j]+=read_real(fp);
        for (int j=0; j<NF; j++) s->mom[j]+=cube(base.a)*match.inertia[j]*read_real(fp);
        Grid *g=&box.grids[i];
        *g=grid;
        g->mode=malloc(modes*sizeof(Mode));
        g->trial=malloc(modes*sizeof(Mode));
        if (!g->mode || !g->trial) fail("Cannot allocate cell modes");
    }
    fclose(fp);
    for (int i=0; i<cells; i++) {
        State s=box.fields[i];
        Grid *g=&box.grids[i];
        double dV[NF], lap[NF], G=0, P=0;
        potential_gradient(p,s.phi,s.phi,dV);
        gradient(&box,i,s.a,match,&G,&P,lap);
        for (int j=0; j<p->ns; j++) {
            Species sp=p->species[j];
            int field=sp.parent;
            double v=s.mom[field]/(cube(s.a)*match.inertia[field]);
            double acceleration=-3*H0*v-dV[field]+lap[field]/(s.a*s.a);
            for (int l=0; l<2*p->nq; l++) {
                double q=g->q[l%p->nq]*(l<p->nq ? -1 : 1)/s.a;
                g->mode[j*2*p->nq+l]=vacuum(sp,q,v,acceleration,H0,p->dressed_start);
            }
        }
    }
    return box;
}

static void snapshot(const Params *p, const Box *box, double t, FILE *out) {
    const Grid *grid=box->grids;
    for (int j=0; j<p->ns; j++) {
        Species sp=p->species[j];
        for (int l=0; l<2*p->nq; l++) {
            int b=l%p->nq, h=l<p->nq ? -1 : 1;
            double k=h*grid->q[b]/box->a, E=hypot(k,sp.mass), occ=0;
            for (int i=0; i<box->cells; i++) {
                Mode r=box->grids[i].mode[j*2*p->nq+l];
                occ+=0.5*(1+(sp.mass*r.x+k*r.z)/E);
            }
            fprintf(out,"%.17g %d %d %.17g %.17g %.17g\n",t,j+1,h,
                    grid->q[b],grid->w[b],occ/box->cells);
        }
    }
}

static void columns(FILE *out) {
    fprintf(out,"# t a H rho P volume_work residual G scalar_ct rho_dr rho_sm rho_s chi2"
                " sigma chi s rms_sigma rms_chi rms_s phase_sigma phase_chi nmin nmax norm_error iterations");
    for (int j=1; j<=NS; j++) fprintf(out," rho_%d n_%d p_%d p_free_%d work_%d drive_%d",j,j,j,j,j,j);
    fprintf(out," Q_pressure_work Q_scalar_work Q_ledger_error");
    fputc('\n',out);
}

static void record_box(const Params *p, Box *box, double t, double work, double E0,
                       double qpressure, double qscalar, double qinitial, int iterations, FILE *out) {
    Average avg=average(p,box,box->a);
    double a=box->a, a3=cube(a), inv=1.0/box->cells;
    double rho[NS]={0}, count[NS]={0}, pressure[NS]={0}, free_p[NS]={0}, power[NS]={0}, drive[NS]={0};
    double phi[NF]={0}, var[NF]={0}, ca=0, sa=0, cb=0, sb=0, rho_s=0, chi2=0;
    double occ_min=1, occ_max=0, norm_err=0;
    #pragma omp parallel for schedule(static) reduction(+:rho[:NS],count[:NS],pressure[:NS],free_p[:NS],power[:NS],drive[:NS],phi[:NF],var[:NF],ca,sa,cb,sb,rho_s,chi2) reduction(min:occ_min) reduction(max:occ_max,norm_err)
    for (int i=0; i<box->cells; i++) {
        State s=box->fields[i];
        s.a=a;
        const Grid *g=&box->grids[i];
        Stress q=stress(p,g,&s,1);
        for (int j=0; j<NF; j++) {
            phi[j]+=s.phi[j]*inv;
            var[j]+=s.phi[j]*s.phi[j]*inv;
        }
        chi2+=s.phi[1]*s.phi[1]*inv;
        rho_s+=q.rho_s*inv;
        if (p->ns) {
            double first=s.phi[0]/p->species[0].f, second=s.phi[1]/p->species[1].f;
            ca+=inv*cos(first); sa+=inv*sin(first);
            cb+=inv*cos(second); sb+=inv*sin(second);
        }
        for (int j=0; j<p->ns; j++) {
            Species sp=p->species[j];
            rho[j]+=q.e[j]*inv;
            count[j]+=q.n[j]*inv;
            pressure[j]+=q.p[j]*inv;
            power[j]+=q.force[j]*q.velocity[sp.parent]*inv;
            drive[j]+=pow(q.velocity[sp.parent]/(2*sp.f*sp.mass),2)*inv;
            for (int l=0; l<2*p->nq; l++) {
                int b=l%p->nq;
                double k=g->q[b]*(l<p->nq ? -1 : 1)/a, E=hypot(k,sp.mass);
                Mode r=g->mode[j*2*p->nq+l];
                double f=0.5*(1+(sp.mass*r.x+k*r.z)/E);
                free_p[j]+=inv*g->w[b]*f*k*k/(3*a3*E);
                occ_min=fmin(occ_min,f); occ_max=fmax(occ_max,f);
                norm_err=fmax(norm_err,fabs(r.x*r.x+r.y*r.y+r.z*r.z-1));
            }
        }
    }
    if (!p->ns) occ_min=0;
    double D0=0, D1=0;
    if (p->ns) {
        double m0=phi[0]/p->species[0].f, m1=phi[1]/p->species[1].f;
        D0=sqrt(fmax(0,2*(1-ca*cos(m0)-sa*sin(m0))));
        D1=sqrt(fmax(0,2*(1-cb*cos(m1)-sb*sin(m1))));
    }
    fprintf(out,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g",
            t,a,avg.hubble,avg.rho,avg.pressure,work,(avg.rho*a3+work)/E0-1,avg.gradient,
            avg.scalar_ct,avg.dark,avg.visible,rho_s,chi2);
    for (int j=0; j<NF; j++) fprintf(out," %.17g",phi[j]);
    for (int j=0; j<NF; j++) fprintf(out," %.17g",sqrt(fmax(0,var[j]-phi[j]*phi[j])));
    fprintf(out," %.17g %.17g %.17g %.17g %.17g %d",D0,D1,occ_min,occ_max,norm_err,iterations);
    for (int j=0; j<NS; j++) fprintf(out," %.17g %.17g %.17g %.17g %.17g %.17g",
            rho[j],count[j],pressure[j],free_p[j],power[j],sqrt(drive[j]));
    fprintf(out," %.17g %.17g %.17g",qpressure,qscalar,(avg.quantum_e*a3+qpressure-qinitial-qscalar)/E0);
    fputc('\n',out);
    fflush(out);
}

int main(int argc,char **argv) {
    if (argc!=5) fail("Usage: production INPUT SEED HISTORY SPECTRA");
    Params p={0};
    State base={0};
    Grid grid={0};
    read_input(argv[1],&p,&base,&grid);
    for (int j=0; j<p.ns; j++) if (p.species[j].width!=0) fail("Collisionless sterile ledger required");
    Box box=prepare_box(&p,base,grid,argv[2]);
    FILE *out=fopen(argv[3],"w"), *spectra=fopen(argv[4],"w");
    if (!out || !spectra) fail("Cannot open records");
    columns(out);
    fprintf(spectra,"# t species h q w f_free\n");
    double work=0, qpressure=0, qscalar=0, dt=(p.t1-p.t0)/p.nsteps;
    Average initial=average(&p,&box,box.a);
    double E0=initial.rho*cube(box.a), qinitial=initial.quantum_e*cube(box.a);
    int iterations=0;
    record_box(&p,&box,p.t0,work,E0,qpressure,qscalar,qinitial,iterations,out);
    for (int k=1; k<=p.nsteps; k++) {
        expand_box(&p,&box,dt/2,&work,&qpressure);
        kick(&p,&box,dt/2);
        int steps=exchange(&p,&box,dt,&qscalar);
        if (steps>iterations) iterations=steps;
        kick(&p,&box,dt/2);
        expand_box(&p,&box,dt/2,&work,&qpressure);
        if (k%p.stride==0 || k==p.nsteps)
            record_box(&p,&box,p.t0+k*dt,work,E0,qpressure,qscalar,qinitial,iterations,out);
    }
    snapshot(&p,&box,p.t1,spectra);
    fclose(out); fclose(spectra);
    for (int i=0; i<box.cells; i++) {
        free(box.grids[i].mode);
        free(box.grids[i].trial);
    }
    free(box.fields); free(box.grids); free(box.links);
    free(grid.q); free(grid.w); free(grid.mode); free(grid.trial);
    return 0;
}
