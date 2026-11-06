/* Finite comoving regulator for homogeneous scalar/Majorana backreaction. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { NF = 3, NS = 4 };
static const double pi = 3.1415926535897932384626433832795;

typedef struct { double x, y, z; } Mode;
typedef struct { int parent; double mass, f, width; } Species;
typedef struct {
    int expand, subtract, dressed_start, nq, ns, nsteps, stride;
    double t0, t1, qmax, mp, lambda, v, g2, sc, c1, m2, mu2, cb, offset, gs2, gamma_s;
    Species species[NS];
} Params;
typedef struct {
    double a, phi[NF], mom[NF], rad_sm, rad_dark;
} State;
typedef struct {
    double inertia[NF], derivative[NF], gravity;
} Matching;
typedef struct {
    double rho, pressure, hubble, velocity[NF], rho_s, scalar_ct, gravity_ct;
    double e[NS], n[NS], p[NS], force[NS];
} Stress;
typedef struct {
    double *q, *w;
    Mode *mode, *trial;
} Grid;

static void fail(const char *message) {
    fprintf(stderr, "%s\n", message);
    exit(1);
}

static double read_real(FILE *f) {
    double v;
    if (fscanf(f, "%lf", &v) != 1 || !isfinite(v)) fail("Invalid input scalar");
    return v;
}

static int read_int(FILE *f) {
    int v;
    if (fscanf(f, "%d", &v) != 1) fail("Invalid input integer");
    return v;
}

static double cube(double x) { return x*x*x; }

static double loop_potential(const Params *p, double sigma) {
    if (sigma <= p->sc) return 0;
    double u = p->g2*sigma*sigma - p->lambda*p->v*p->v;
    if (u == 0) return 0;
    return u*u*(log(fabs(u)/p->mu2)-1.5)/(64*pi*pi);
}

static double loop_derivative(const Params *p, double sigma) {
    if (sigma <= p->sc) return 0;
    double u = p->g2*sigma*sigma - p->lambda*p->v*p->v;
    if (u == 0) return 0;
    return u*2*p->g2*sigma*(log(fabs(u)/p->mu2)-1)/(32*pi*pi);
}

static double potential(const Params *p, const double field[NF]) {
    double s = field[0], x = field[1], y = field[2], d = s-p->sc, x2 = x*x;
    double u = x2-p->v*p->v;
    return p->lambda*u*u/4 + p->g2*s*s*x2/2 - p->cb*x + p->c1*d
           - p->m2*d*d/2 + p->offset + p->gs2*x2*y*y/2 + loop_potential(p,s);
}

/* Symmetric discrete gradient for all polynomial terms. Its contraction
 * with the displacement gives the potential change, including the mixed
 * spectator coupling. The one-field loop term uses a divided difference. */
static void potential_gradient(const Params *p, const double a[NF],
                               const double b[NF], double result[NF]) {
    double mid[NF], square[NF];
    for (int j=0; j<NF; j++) {
        mid[j] = (a[j]+b[j])/2;
        square[j] = (a[j]*a[j]+b[j]*b[j])/2;
    }
    double ds = b[0]-a[0];
    double loop = fabs(ds) > 1e-9*(1+fabs(mid[0]))
        ? (loop_potential(p,b[0])-loop_potential(p,a[0]))/ds
        : loop_derivative(p,mid[0]);
    result[0] = p->g2*mid[0]*square[1] + p->c1 - p->m2*(mid[0]-p->sc) + loop;
    result[1] = mid[1]*(p->lambda*(square[1]-p->v*p->v) + p->g2*square[0]
                         + p->gs2*square[2]) - p->cb;
    result[2] = p->gs2*square[1]*mid[2];
}

static Matching matching(const Params *p, double a) {
    Matching m = {{1,1,1}, {0,0,0}, 0};
    if (!p->subtract) return m;
    double cutoff = p->qmax/a;
    for (int j=0; j<p->ns; j++) {
        Species sp = p->species[j];
        if (sp.mass == 0) continue;
        double beta = cutoff/hypot(cutoff, sp.mass);
        double b3 = cube(beta);
        double integral = asinh(cutoff/sp.mass)-beta;
        double factor = sp.mass*sp.mass/(8*pi*pi);
        double wave = factor/(sp.f*sp.f);
        m.inertia[sp.parent] -= wave*integral;
        m.derivative[sp.parent] += wave*b3;
        m.gravity += factor*(integral-b3/3);
    }
    if (m.inertia[0] <= 0 || m.inertia[1] <= 0) fail("Nonpositive scalar inertia");
    return m;
}

static Stress stress(const Params *p, const Grid *g, const State *s, int detail) {
    Stress out = {0};
    Matching match = matching(p,s->a);
    double a3=cube(s->a), a4=a3*s->a, ia=1/s->a;
    double V = potential(p,s->phi);
    out.rho = V + (s->rad_sm+s->rad_dark)/a4;
    out.pressure = -V + (s->rad_sm+s->rad_dark)/(3*a4);
    for (int j=0; j<NF; j++) {
        double v=s->mom[j]/(a3*match.inertia[j]);
        double kin=0.5*match.inertia[j]*v*v;
        out.velocity[j]=v;
        out.rho+=kin;
        out.pressure+=kin+match.derivative[j]*v*v/6;
        out.scalar_ct+=0.5*(1-match.inertia[j])*v*v;
    }
    out.rho_s = 0.5*(out.velocity[2]*out.velocity[2]
                    + p->gs2*s->phi[1]*s->phi[1]*s->phi[2]*s->phi[2]);

    double energies[NS]={0}, numbers[NS]={0}, pressures[NS]={0}, forces[NS]={0};
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:energies[:NS],numbers[:NS],pressures[:NS],forces[:NS])
    for (int j=0; j<p->ns; j++) {
        for (int i=0; i<2*p->nq; i++) {
            Species sp = p->species[j];
            int bin=i%p->nq;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]*ia, w=g->w[bin];
            double E=sqrt(k*k+sp.mass*sp.mass);
            Mode r=g->mode[j*2*p->nq+i];
            double f=0.5*(1+(sp.mass*r.x+k*r.z)/E);
            energies[j]+=w*E*f;
            pressures[j]+=w*k*(r.z+k/E)/6;
            if (detail) {
                numbers[j]+=w*f;
                forces[j]-=w*sp.mass*r.y/(2*sp.f);
            }
        }
    }
    for (int j=0; j<p->ns; j++) {
        out.e[j]=energies[j]/a3;
        out.n[j]=numbers[j]/a3;
        out.p[j]=pressures[j]/a3;
        out.force[j]=forces[j]/a3;
        out.rho+=out.e[j];
        out.pressure+=out.p[j];
    }
    if (!(out.rho>0) || !isfinite(out.pressure)) fail("Invalid density or pressure");
    out.hubble=p->expand ? sqrt(out.rho/(3*p->mp*p->mp+match.gravity/2)) : 0;
    out.gravity_ct=0.5*match.gravity*out.hubble*out.hubble;
    return out;
}

/* Expansion with canonical scalar momenta and the spin state held fixed.
 * The explicit midpoint uses a separate integral of the analytic pressure. */
static void expansion(const Params *p, const Grid *g, State *s, double dt, double *work) {
    if (!p->expand) return;
    double before=s->a;
    Stress begin=stress(p,g,s,0);
    s->a=before*exp(0.5*dt*begin.hubble);
    Stress middle=stress(p,g,s,0);
    s->a=before*exp(dt*middle.hubble);
    *work+=middle.pressure*cube(before)*expm1(3*dt*middle.hubble);
}

static inline Mode rotate(Mode r, double x, double z, double dt, double *dz) {
    double omega=sqrt(x*x+z*z);
    if (omega==0) { *dz=0; return r; }
    double ux=x/omega, uz=z/omega, angle=2*omega*dt;
    double c=cos(angle), sn=sin(angle);
    double dot=ux*r.x+uz*r.z;
    double cm1 = fabs(angle)<1e-4 ? -0.5*angle*angle*(1-angle*angle/12) : c-1;
    *dz=cm1*(r.z-uz*dot)+sn*ux*r.y;
    return (Mode){r.x+cm1*(r.x-ux*dot)-sn*uz*r.y,
                  c*r.y+sn*(uz*r.x-ux*r.z), r.z+*dz};
}

static void trial_impulse(const Params *p, Grid *g, double a,
                          const double v[NF], double dt, double impulse[NF]) {
    double sum[NF]={0}, ia=1/a;
    #pragma omp parallel for collapse(2) schedule(static) reduction(+:sum[:NF])
    for (int j=0; j<p->ns; j++) {
        for (int i=0; i<2*p->nq; i++) {
            Species sp=p->species[j];
            int bin=i%p->nq, id=j*2*p->nq+i;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]*ia;
            double dz;
            g->trial[id]=rotate(g->mode[id],sp.mass,k+v[sp.parent]/(2*sp.f),dt,&dz);
            sum[sp.parent]+=g->w[bin]*dz/(4*sp.f);
        }
    }
    memcpy(impulse,sum,sizeof sum);
}

/* Solve the three midpoint velocities. At a converged fixed point the
 * mode rotation supplies both the scalar kick and the opposite work. */
static int matter(const Params *p, Grid *g, State *s, double dt) {
    Matching m=matching(p,s->a);
    double a3=cube(s->a);
    double v[NF], inv_mass[NF], next_field[NF], next_mom[NF];
    for (int j=0; j<NF; j++) {
        inv_mass[j]=1/(a3*m.inertia[j]);
        v[j]=s->mom[j]*inv_mass[j];
    }
    for (int it=1; it<=16; it++) {
        for (int j=0; j<NF; j++) next_field[j]=s->phi[j]+dt*v[j];
        double force[NF], impulse[NF];
        potential_gradient(p,s->phi,next_field,force);
        trial_impulse(p,g,s->a,v,dt,impulse);
        double err=0;
        for (int j=0; j<NF; j++) {
            next_mom[j]=s->mom[j]-a3*dt*force[j]+impulse[j];
            double vnew=0.5*(s->mom[j]+next_mom[j])*inv_mass[j];
            err=fmax(err,fabs(vnew-v[j])/(1+fabs(vnew)+fabs(s->mom[j]*inv_mass[j])));
            v[j]=vnew;
        }
        if (err<2e-12) {
            memcpy(s->phi,next_field,sizeof next_field);
            memcpy(s->mom,next_mom,sizeof next_mom);
            Mode *swap=g->mode; g->mode=g->trial; g->trial=swap;
            return it;
        }
    }
    fail("Midpoint scalar/fermion update did not converge");
    return 0;
}

/* Vacuum decay control in the instantaneous free-mass basis. The same
 * occupation loss goes to visible radiation. Spectator damping uses its
 * exact kinetic-energy loss in this fixed-a split substep. */
static void decays(const Params *p, Grid *g, State *s, double dt) {
    if (p->gamma_s>0) {
        double loss=-expm1(-2*p->gamma_s*dt);
        s->rad_dark+=0.5*s->mom[2]*s->mom[2]*loss/(s->a*s->a);
        s->mom[2]*=exp(-p->gamma_s*dt);
    }
    double sum=0, ia=1/s->a;
    for (int j=0; j<p->ns; j++) {
        Species sp=p->species[j];
        if (sp.width==0 || sp.mass==0) continue;
        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int i=0; i<2*p->nq; i++) {
            int bin=i%p->nq, id=j*2*p->nq+i;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]*ia;
            double E=sqrt(k*k+sp.mass*sp.mass), ux=sp.mass/E, uz=k/E;
            Mode r=g->mode[id];
            double parallel=ux*r.x+uz*r.z;
            double loss=-expm1(-sp.width*sp.mass*dt/E), transverse=sqrt(1-loss);
            double newer=parallel-loss*(parallel+1);
            g->mode[id]=(Mode){transverse*(r.x-parallel*ux)+newer*ux,
                              transverse*r.y,transverse*(r.z-parallel*uz)+newer*uz};
            sum+=0.5*g->w[bin]*E*loss*(1+parallel);
        }
    }
    s->rad_sm+=s->a*sum;
}

static Mode vacuum(Species sp, double k, double v, double acceleration,
                   double H, int dressed) {
    if (sp.mass==0) return (Mode){0,0,k<0 ? 1 : -1};
    double z=k+(dressed ? v/(2*sp.f) : 0);
    double omega=hypot(z,sp.mass);
    if (omega==0) return (Mode){0,0,k<0 ? 1 : -1};
    double y=dressed ? -sp.mass*(acceleration/(2*sp.f)-H*k)/(2*cube(omega)) : 0;
    double norm=hypot(1,y);
    return (Mode){-sp.mass/(omega*norm),y/norm,-z/(omega*norm)};
}

static void read_input(const char *name, Params *p, State *s, Grid *g) {
    FILE *f=fopen(name,"r");
    if (!f) fail("Cannot open input");
    p->expand=read_int(f); p->subtract=read_int(f); p->dressed_start=read_int(f);
    p->nq=read_int(f); p->ns=read_int(f); p->nsteps=read_int(f); p->stride=read_int(f);
    p->t0=read_real(f); p->t1=read_real(f); p->qmax=read_real(f);
    p->mp=read_real(f); p->lambda=read_real(f); p->v=read_real(f); p->g2=read_real(f);
    p->sc=read_real(f); p->c1=read_real(f); p->m2=read_real(f); p->mu2=read_real(f);
    p->cb=read_real(f); p->offset=read_real(f); p->gs2=read_real(f); p->gamma_s=read_real(f);
    s->a=read_real(f); double H0=read_real(f);
    s->rad_sm=cube(s->a)*s->a*read_real(f);
    s->rad_dark=cube(s->a)*s->a*read_real(f);
    if (p->nq<1 || p->ns<0 || p->ns>NS || p->nsteps<1 || p->stride<1
        || p->t1<=p->t0 || p->qmax<=0 || s->a<=0 || p->mu2<=0 || p->mp<=0)
        fail("Invalid run parameters");
    double v[NF];
    for (int j=0; j<NF; j++) s->phi[j]=read_real(f);
    for (int j=0; j<NF; j++) v[j]=read_real(f);
    for (int j=0; j<p->ns; j++) {
        Species *sp=&p->species[j];
        sp->parent=read_int(f); sp->mass=read_real(f); sp->f=read_real(f); sp->width=read_real(f);
        if (sp->parent<0 || sp->parent>=NF-1 || sp->mass<0 || sp->f<=0 || sp->width<0)
            fail("Invalid species");
    }
    fclose(f);

    size_t points=p->nq, nmodes=2*points*(p->ns ? p->ns : 1);
    g->q=malloc(points*sizeof *g->q); g->w=malloc(points*sizeof *g->w);
    g->mode=malloc(nmodes*sizeof *g->mode); g->trial=malloc(nmodes*sizeof *g->trial);
    if (!g->q || !g->w || !g->mode || !g->trial) fail("Grid allocation failed");
    double dq=p->qmax/p->nq;
    for (int i=0; i<p->nq; i++) {
        g->q[i]=dq*(i+0.5);
        g->w[i]=(cube(dq*(i+1))-cube(dq*i))/(6*pi*pi);
    }
    Matching match=matching(p,s->a);
    double force[NF];
    potential_gradient(p,s->phi,s->phi,force);
    for (int j=0; j<NF; j++) {
        s->mom[j]=cube(s->a)*match.inertia[j]*v[j];
        force[j]=-3*H0*v[j]-force[j];
    }
    for (int j=0; j<p->ns; j++) {
        Species sp=p->species[j];
        for (int i=0; i<2*p->nq; i++) {
            double k=(i<p->nq ? -1 : 1)*g->q[i%p->nq]/s->a;
            g->mode[j*2*p->nq+i]=vacuum(sp,k,v[sp.parent],force[sp.parent],H0,p->dressed_start);
        }
    }
}

static void record(FILE *f, FILE *spectra, const Params *p, const Grid *g, const State *s,
                   double time, double work, double E0, int iterations) {
    Stress stress_now=stress(p,g,s,1);
    Matching match=matching(p,s->a);
    double a3=cube(s->a), a4=a3*s->a, H=stress_now.hubble;
    double force[NF];
    potential_gradient(p,s->phi,s->phi,force);
    for (int j=0; j<p->ns; j++) force[p->species[j].parent]+=stress_now.force[j];
    double acceleration[NF];
    for (int j=0; j<NF; j++)
        acceleration[j]=-(force[j]+H*(3*match.inertia[j]+match.derivative[j])*stress_now.velocity[j])
                         /match.inertia[j];

    double norm_max=0, nmin=1, nmax=0, dressed_e[NS]={0}, dressed_n[NS]={0}, shifts[NF]={0};
    for (int j=0; j<p->ns; j++) {
        Species sp=p->species[j];
        for (int i=0; i<2*p->nq; i++) {
            int h=i<p->nq ? -1 : 1, bin=i%p->nq;
            double k=h*g->q[bin]/s->a, w=g->w[bin];
            Mode r=g->mode[j*2*p->nq+i];
            double E=hypot(k,sp.mass);
            double f0=0.5*(1+(sp.mass*r.x+k*r.z)/E);
            Mode ad=vacuum(sp,k,stress_now.velocity[sp.parent],acceleration[sp.parent],H,1);
            double f1=0.5*(1-ad.x*r.x-ad.y*r.y-ad.z*r.z);
            norm_max=fmax(norm_max,sqrt(r.x*r.x+r.y*r.y+r.z*r.z));
            nmin=fmin(nmin,fmin(f0,f1)); nmax=fmax(nmax,fmax(f0,f1));
            dressed_e[j]+=w*E*f1/a3; dressed_n[j]+=w*f1/a3;
            shifts[sp.parent]+=w*r.z/(4*sp.f);
            if (spectra)
                fprintf(spectra,"%d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                        j+1,h,g->q[bin],w,f0,f1,r.x,r.y,r.z);
        }
    }
    if (p->ns==0) nmin=0;
    double residual=(stress_now.rho*a3+work)/E0-1;
    fprintf(f,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g",time,s->a,H,stress_now.rho,
            stress_now.pressure,work,residual);
    for (int j=0; j<NF; j++) fprintf(f," %.17g",s->phi[j]);
    for (int j=0; j<NF; j++) fprintf(f," %.17g",stress_now.velocity[j]);
    fprintf(f," %.17g %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g",stress_now.rho_s,
            s->rad_sm/a4,s->rad_dark/a4,stress_now.scalar_ct,stress_now.gravity_ct,
            iterations,norm_max,nmin,nmax);
    for (int j=0; j<NF; j++) fprintf(f," %.17g",s->mom[j]-shifts[j]);
    for (int j=0; j<NS; j++)
        fprintf(f," %.17g %.17g %.17g %.17g %.17g %.17g",stress_now.e[j],stress_now.n[j],
                stress_now.p[j],dressed_e[j],dressed_n[j],stress_now.force[j]);
    fputc('\n',f);
    fflush(f);
}

int main(int argc, char **argv) {
    if (argc!=4) fail("Usage: coupled INPUT HISTORY SPECTRA");
    Params p={0}; State s={0}; Grid g={0};
    read_input(argv[1],&p,&s,&g);
    FILE *out=fopen(argv[2],"w"), *spectra=fopen(argv[3],"w");
    if (!out || !spectra) fail("Cannot open output");
    fprintf(out,"# t a H rho pressure work residual sigma chi s dsigma dchi ds rho_s rho_sm rho_dr"
                " scalar_ct gravity_ct max_iterations norm_max f_min f_max Q_sigma Q_chi Q_s");
    for (int j=1; j<=NS; j++)
        fprintf(out," rho_%d n_%d p_%d dressed_rho_%d dressed_n_%d J_%d",j,j,j,j,j,j);
    fputc('\n',out);
    fprintf(spectra,"# species h q w f_free f_dressed r_x r_y r_z\n");
    double dt=(p.t1-p.t0)/p.nsteps, work=0;
    double E0=stress(&p,&g,&s,0).rho*cube(s.a);
    int iterations=0;
    record(out,NULL,&p,&g,&s,p.t0,work,E0,0);
    for (int k=1; k<=p.nsteps; k++) {
        expansion(&p,&g,&s,dt/2,&work);
        decays(&p,&g,&s,dt/2);
        int n=matter(&p,&g,&s,dt);
        if (n>iterations) iterations=n;
        decays(&p,&g,&s,dt/2);
        expansion(&p,&g,&s,dt/2,&work);
        if (k%p.stride==0 || k==p.nsteps)
            record(out,k==p.nsteps ? spectra : NULL,&p,&g,&s,p.t0+k*dt,work,E0,iterations);
    }
    fclose(out); fclose(spectra);
    free(g.q); free(g.w); free(g.mode); free(g.trial);
    return 0;
}
