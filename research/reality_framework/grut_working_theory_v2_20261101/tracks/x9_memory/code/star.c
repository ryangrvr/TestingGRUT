/* Unitary split evolution of two preparations in the same Hermitian star. */
#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

enum { COLS=12 };

static inline double norm2(double complex z) {
    return creal(z)*creal(z)+cimag(z)*cimag(z);
}

static inline double x_value(const double complex *v) {
    return 2*creal(conj(v[0])*v[1]);
}

static inline double y_value(const double complex *v) {
    return 2*cimag(conj(v[0])*v[1]);
}

static inline void spin(double complex *v,double complex r0,double complex r1,
                        double complex off) {
    double complex s=v[0], t=v[1];
    v[0]=r0*s+off*t;
    v[1]=off*s+r1*t;
}

static void sample(int n,double k,const double *freq,const double *g,
                   const double *path,const double complex *v,double wx,double wy,double *row) {
    double norm=norm2(v[0])+norm2(v[1]), z=norm2(v[0])-norm2(v[1]);
    double S=path[0]*x_value(v)+k*z, E=hypot(k,path[3]);
    double number=(norm+(k*z+path[3]*x_value(v))/E)/2;
    double eb=0, mix=0, pos=0, neg=0, charge=0;
    for (int b=0;b<2;b++) {
        double complex coupling=0;
        const double complex *bath=v+2+b*n;
        for (int i=0;i<n;i++) {
            double p=norm2(bath[i]);
            coupling+=g[b*n+i]*bath[i];
            eb+=freq[i]*p;
            if (freq[i]>0) {
                pos+=p;
                charge+=(b ? -1 : 1)*p;
            } else {
                neg+=p;
            }
        }
        mix+=2*creal(conj(v[b])*coupling);
    }
    row[0]=norm+pos+neg; row[1]=norm; row[2]=number;
    row[3]=pos; row[4]=neg; row[5]=charge;
    row[6]=S; row[7]=eb; row[8]=mix; row[9]=S+eb+mix;
    row[10]=wx; row[11]=wy;
}

/* path columns are matched mass, phase speed/2, mass derivative and physical mass.
   Each state contains two sterile amplitudes followed by two n-mode baths. */
int propagate(int n,int steps,int stride,double dt,double k,const double *freq,
              const double *g,const double *path,double complex *v,double *out) {
    if (n<1 || steps<1 || stride<1 || dt<=0) return -1;
    double complex *free_phase=malloc((size_t)n*sizeof *free_phase);
    double *direction=malloc((size_t)(2*n)*sizeof *direction);
    if (!free_phase || !direction) {
        free(free_phase);free(direction);
        return -2;
    }
    double cs[2], sn[2];
    for (int b=0;b<2;b++) {
        double size=0;
        for (int i=0;i<n;i++) size+=g[b*n+i]*g[b*n+i];
        size=sqrt(size);
        cs[b]=cos(size*dt); sn[b]=sin(size*dt);
        for (int i=0;i<n;i++) direction[b*n+i]=size>0 ? g[b*n+i]/size : 0.;
    }
    for (int i=0;i<n;i++) {
        double a=-freq[i]*dt/2;
        free_phase[i]=cos(a)+I*sin(a);
    }
    const int dimension=2+2*n;
    double wx[2]={0}, wy[2]={0};
    int at=0;
    for (int j=0;j<2;j++)
        sample(n,k,freq,g,path,v+j*dimension,0,0,out+(at*2+j)*COLS);
    at++;
    for (int p=0;p<steps;p++) {
        const double *old=path+4*p, *next=old+4;
        double m=(old[0]+next[0])/2, delta=(old[1]+next[1])/2;
        double z=k+delta, energy=hypot(m,z);
        double c=cos(energy*dt/2), s=energy>0 ? sin(energy*dt/2)/energy : dt/2;
        double complex r0=c-I*z*s, r1=c+I*z*s, off=-I*m*s;
        double angle=-delta*dt/2;
        double complex axial[2]={cos(angle)+I*sin(angle),cos(angle)-I*sin(angle)};
        for (int j=0;j<2;j++) {
            double complex *u=v+j*dimension;
            double before_x=x_value(u), before_y=y_value(u);
            spin(u,r0,r1,off);
            for (int b=0;b<2;b++) {
                double complex *bath=u+2+b*n, phase=axial[b];
                const double *dir=direction+b*n;
                double complex aligned=0;
                for (int i=0;i<n;i++) {
                    bath[i]*=free_phase[i]*phase;
                    aligned+=dir[i]*bath[i];
                }
                double complex increment=(cs[b]-1)*aligned-I*sn[b]*u[b];
                u[b]=cs[b]*u[b]-I*sn[b]*aligned;
                for (int i=0;i<n;i++)
                    bath[i]=(bath[i]+dir[i]*increment)*free_phase[i]*phase;
            }
            spin(u,r0,r1,off);
            wx[j]+=dt*(old[2]*before_x+next[2]*x_value(u))/2;
            wy[j]-=dt*(old[1]*old[0]*before_y+next[1]*next[0]*y_value(u));
        }
        if ((p+1)%stride==0 || p+1==steps) {
            for (int j=0;j<2;j++)
                sample(n,k,freq,g,next,v+j*dimension,wx[j],wy[j],out+(at*2+j)*COLS);
            at++;
        }
    }
    free(free_phase); free(direction);
    return at;
}

static void pade_factor(int n,double complex z,double mass,double k,double shift,
                        const double *freq,const double *g,double complex *d,double complex *v) {
    const int dim=2+2*n;
    double complex a[2], right[2][2];
    for (int b=0;b<2;b++) {
        double signed_shift=b ? -shift : shift, signed_k=b ? -k : k;
        const double *g_b=g+b*n;
        double complex *d_b=d+b*n, sigma=0;
        double complex *u0=v+2+b*n, *u1=v+dim+2+b*n;
        double complex p0=v[b], p1=v[dim+b];
        for (int i=0;i<n;i++) {
            double complex f=1./(z-freq[i]-signed_shift);
            double complex gf=g_b[i]*f;
            d_b[i]=f;
            sigma+=g_b[i]*gf;
            p0+=gf*u0[i];
            p1+=gf*u1[i];
        }
        a[b]=z-signed_k-signed_shift-sigma;
        right[0][b]=p0;right[1][b]=p1;
    }
    double complex inv=1./(a[0]*a[1]-mass*mass);
    for (int j=0;j<2;j++) {
        double complex *u=v+j*dim;
        double complex s[2]={(a[1]*right[j][0]+mass*right[j][1])*inv,
                             (a[0]*right[j][1]+mass*right[j][0])*inv};
        for (int b=0;b<2;b++) {
            u[b]=2*z*s[b]-u[b];
            double complex *bath=u+2+b*n;
            const double *g_b=g+b*n;
            const double complex *d_b=d+b*n;
            for (int i=0;i<n;i++) bath[i]=2*z*d_b[i]*(bath[i]+g_b[i]*s[b])-bath[i];
        }
    }
}

static double axial_charge(int n,const double complex *v) {
    double charge=norm2(v[0])-norm2(v[1]);
    for (int i=0;i<n;i++) charge+=norm2(v[2+i])-norm2(v[2+n+i]);
    return charge;
}

int propagate_pade(int n,int steps,int stride,double dt,double k,const double *freq,
                   const double *g,const double *path,double complex *v,double *out,double *checks) {
    if (n<1 || steps<1 || stride<1 || dt<=0) return -1;
    double complex *d=malloc((size_t)(2*n)*sizeof *d);
    if (!d) return -2;
    /* Poles of the [3/3] Padé exponential in its dimensionless argument. */
    const double complex roots[3]={4.6443707092521711,
        3.6778146453739096+3.5087619195674429*I,3.6778146453739096-3.5087619195674429*I};
    const int dim=2+2*n;
    double wx[2]={0}, wy[2]={0}, continuous[2]={0};
    double before_z[2]={axial_charge(n,v),axial_charge(n,v+dim)};
    int at=0;
    for (int j=0;j<2;j++) sample(n,k,freq,g,path,v+j*dim,0,0,out+(at*2+j)*COLS);
    at++;
    for (int p=0;p<steps;p++) {
        const double *old=path+4*p, *next=old+4;
        double m=(old[0]+next[0])/2, delta=(old[1]+next[1])/2;
        double x[2]={x_value(v),x_value(v+dim)}, y[2]={y_value(v),y_value(v+dim)};
        for (int i=0;i<3;i++) pade_factor(n,I*roots[i]/dt,m,k,delta,freq,g,d,v);
        for (int j=0;j<2;j++) {
            double complex *u=v+j*dim;
            double z=axial_charge(n,u);
            wx[j]+=(next[0]-old[0])*(x[j]+x_value(u))/2;
            wy[j]-=delta*(z-before_z[j]);
            continuous[j]-=dt*(old[1]*old[0]*y[j]+next[1]*next[0]*y_value(u));
            before_z[j]=z;
        }
        if ((p+1)%stride==0 || p+1==steps) {
            for (int j=0;j<2;j++) sample(n,k,freq,g,next,v+j*dim,wx[j],wy[j],out+(at*2+j)*COLS);
            at++;
        }
    }
    for (int j=0;j<2;j++) checks[j]=continuous[j]-wy[j];
    free(d);
    return at;
}
