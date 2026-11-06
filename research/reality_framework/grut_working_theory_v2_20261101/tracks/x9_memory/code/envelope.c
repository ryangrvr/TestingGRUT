/* Exact response envelope at fixed X5 midpoint scalar history. */
#define main x5_original_main
#include "../../x5_coupled/code/coupled.c"
#undef main

static double window(double p) {
    double x=fmin(1.,fmax(0.,(p-3000.)/3000.));
    return 1-x*x*(3-2*x);
}

static void coefficient(const Params *p,const Grid *g,const State *s,const Mode *dir,
                        double *out,double *energy,double *number) {
    double a3=cube(s->a);
    *energy=0;*number=0;
    #pragma omp parallel for collapse(2) reduction(+:energy[:1],number[:1]) schedule(static)
    for (int j=0;j<p->ns;j++) {
        for (int i=0;i<2*p->nq;i++) {
            int bin=i%p->nq,id=j*2*p->nq+i;
            double k=(i<p->nq ? -1 : 1)*g->q[bin]/s->a,M=p->species[j].mass;
            double E=hypot(k,M),w=g->w[bin]/a3;
            Mode r=g->mode[id],u=dir[id];
            double f=0.5*(1+(M*r.x+k*r.z)/E);
            double d=w*(M*u.x+k*u.z)/E;
            out[2*id]=d;out[2*id+1]=E*d;
            *number+=w*f;*energy+=w*E*f;
        }
    }
}

int main(int argc,char **argv) {
    if (argc!=3) fail("Usage: envelope INPUT OUTPUT");
    Params p={0}; State s={0}; Grid g={0};
    read_input(argv[1],&p,&s,&g);
    const double split=0.16,start_window=p.t1-0.10,goal=(p.t1-p.t0)/p.nsteps;
    const int modes=p.ns*2*p.nq;
    Mode *response=calloc((size_t)modes,sizeof *response);
    double *amplitude=calloc((size_t)modes,sizeof *amplitude);
    double *last=calloc((size_t)(2*modes),sizeof *last);
    double *current=calloc((size_t)(2*modes),sizeof *current);
    double *integrated=calloc((size_t)(2*modes),sizeof *integrated);
    if (!response || !amplitude || !last || !current || !integrated) fail("Envelope allocation failed");
    double total_work=0,n=0,e=0,mean_n=0,mean_e=0;
    const double cuts[]={p.t0,split,start_window,p.t1};
    for (int part=0;part<3;part++) {
        int nsteps=(int)ceil((cuts[part+1]-cuts[part])/goal);
        double dt=(cuts[part+1]-cuts[part])/nsteps;
        if (part==1) {
            for (int i=0;i<modes;i++) {
                Mode r=g.mode[i];
                double bound=sqrt(fmax(0.,1-r.x*r.x-r.z*r.z));
                amplitude[i]=window(g.q[i%p.nq]/s.a)*bound;
                response[i]=(Mode){0,amplitude[i],0};
            }
        }
        if (part==2) coefficient(&p,&g,&s,response,last,&e,&n);
        for (int step=0;step<nsteps;step++) {
            expansion(&p,&g,&s,dt/2,&total_work);
            decays(&p,&g,&s,dt/2);
            double mom[NF];
            memcpy(mom,s.mom,sizeof mom);
            matter(&p,&g,&s,dt);
            if (part) {
                Matching match=matching(&p,s.a);
                double delta[NS];
                for (int j=0;j<p.ns;j++) {
                    int parent=p.species[j].parent;
                    double velocity=(s.mom[parent]+mom[parent])/(2*cube(s.a)*match.inertia[parent]);
                    delta[j]=velocity/(2*p.species[j].f);
                }
                #pragma omp parallel for collapse(2) schedule(static)
                for (int j=0;j<p.ns;j++) {
                    for (int i=0;i<2*p.nq;i++) {
                        int id=j*2*p.nq+i;
                        double k=(i<p.nq ? -1 : 1)*g.q[i%p.nq]/s.a, dz;
                        if (amplitude[id]>0)
                            response[id]=rotate(response[id],p.species[j].mass,k+delta[j],dt,&dz);
                    }
                }
            }
            decays(&p,&g,&s,dt/2);
            expansion(&p,&g,&s,dt/2,&total_work);
            if (part==2) {
                double next_e,next_n;
                coefficient(&p,&g,&s,response,current,&next_e,&next_n);
                for (int i=0;i<2*modes;i++) integrated[i]+=0.5*dt*(last[i]+current[i]);
                mean_n+=0.5*dt*(n+next_n);mean_e+=0.5*dt*(e+next_e);
                n=next_n;e=next_e;
                double *swap=last;last=current;current=swap;
            }
        }
    }
    FILE *out=fopen(argv[2],"w");
    if (!out) fail("Cannot write envelope");
    fprintf(out,"# species endpoint_number endpoint_energy mean_number mean_energy tangent_error active\n");
    for (int j=0;j<p.ns;j++) {
        double dn=0,de=0,wn=0,we=0,error=0;
        int active=0;
        for (int i=0;i<2*p.nq;i++) {
            int id=j*2*p.nq+i;
            dn+=fabs(last[2*id]); de+=fabs(last[2*id+1]);
            wn+=fabs(integrated[2*id]); we+=fabs(integrated[2*id+1]);
            Mode u=response[id];
            double len=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
            if (amplitude[id]>0) {
                active++;
                error=fmax(error,fabs(len/amplitude[id]-1));
            }
        }
        fprintf(out,"%d %.17g %.17g %.17g %.17g %.17g %d\n",j+1,dn,de,wn/0.10,we/0.10,error,active);
    }
    fprintf(out,"0 %.17g %.17g %.17g %.17g 0 0\n",n,e,mean_n/0.10,mean_e/0.10);
    fclose(out);
    free(g.q);free(g.w);free(g.mode);free(g.trial);
    free(response);free(amplitude);free(last);free(current);free(integrated);
    return 0;
}
