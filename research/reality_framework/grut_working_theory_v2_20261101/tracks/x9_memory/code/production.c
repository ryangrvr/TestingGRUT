/* Full homogeneous sterile-scalar replay with a compact phase preparation. */
#define main x5_original_main
#include "../../x5_coupled/code/coupled.c"
#undef main

static void columns(FILE *out) {
    fprintf(out,"# t a H rho pressure work residual sigma chi s dsigma dchi ds rho_s rho_sm rho_dr"
                " scalar_ct gravity_ct max_iterations norm_max f_min f_max Q_sigma Q_chi Q_s");
    for (int j=1;j<=NS;j++) fprintf(out," rho_%d n_%d p_%d dressed_rho_%d dressed_n_%d J_%d",j,j,j,j,j,j);
    fputc('\n',out);
}

static void prepare(const Params *p,Grid *g,double a) {
    for (int j=0;j<p->ns;j++) {
        for (int i=0;i<2*p->nq;i++) {
            double momentum=g->q[i%p->nq]/a;
            double x=fmin(1.,fmax(0.,(momentum-3000.)/3000.));
            double same=x*x*(3-2*x);
            g->mode[j*2*p->nq+i].y*=2*same-1;
        }
    }
}

static int step(const Params *p,Grid *g,State *s,double dt,double *work) {
    expansion(p,g,s,dt/2,work);
    decays(p,g,s,dt/2);
    int it=matter(p,g,s,dt);
    decays(p,g,s,dt/2);
    expansion(p,g,s,dt/2,work);
    return it;
}

int main(int argc,char **argv) {
    if (argc!=6) fail("Usage: production INPUT HISTORY MATCH REFLECT FINE_STRIDE");
    Params p={0}; State s={0}; Grid g={0};
    read_input(argv[1],&p,&s,&g);
    const int change=atoi(argv[4]), fine_stride=atoi(argv[5]);
    if (fine_stride<1 || change<0 || change>1 || !(p.t0<0.16 && 0.42<p.t1))
        fail("Invalid production interval");
    for (int j=0;j<p.ns;j++) if (p.species[j].width!=0)
        fail("This replay requires the collisionless sterile model");
    FILE *out=fopen(argv[2],"w"), *match=fopen(argv[3],"w");
    if (!out || !match) fail("Cannot write production output");
    columns(out); columns(match);
    double work=0, E0=stress(&p,&g,&s,0).rho*cube(s.a), target=(p.t1-p.t0)/p.nsteps;
    const double split=0.16;
    int iterations=0;
    record(out,NULL,&p,&g,&s,p.t0,work,E0,0);
    for (int part=0;part<2;part++) {
        double begin=part ? split : p.t0, end=part ? p.t1 : split;
        int steps=(int)ceil((end-begin)/target);
        double dt=(end-begin)/steps;
        for (int i=1;i<=steps;i++) {
            int it=step(&p,&g,&s,dt,&work);
            if (it>iterations) iterations=it;
            double t=begin+i*dt;
            int stride=t<=0.423 ? fine_stride : p.stride;
            if (i%stride==0 || i==steps)
                record(out,NULL,&p,&g,&s,t,work,E0,iterations);
        }
        if (!part) {
            record(match,NULL,&p,&g,&s,split,work,E0,iterations);
            if (change) prepare(&p,&g,s.a);
            record(match,NULL,&p,&g,&s,split,work,E0,iterations);
        }
    }
    fclose(out); fclose(match);
    free(g.q); free(g.w); free(g.mode); free(g.trial);
    return 0;
}
