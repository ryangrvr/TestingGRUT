/* Full six-real-scalar Coulomb-gauge tangent evolution on the homogeneous
 * nonminimal background. The longitudinal mode is in these six variables.
 * Metric/inflaton/radiation perturbations and transverse vectors excluded.
 */
#define main original_background_main
#include "../../grut_reheating_preparation_20261026/code/reheat.c"
#undef main
#include <stdint.h>
#include <string.h>
static int nk,nvar;
static double eg,*waves;
#ifdef NEUTRAL_SECTOR
#define NMODE 2
#else
#define NMODE 12
#endif
#define BLOCK (NMODE*NMODE)

static void perturbation_matrix(const double *v,double mat[6][6],double g[6],double gp[6],double ward[6]) {
    double complex z[3],vel[3],force[3];double pot,f,fd,R,hd;
    matter(v,z,vel,force,&pot);geometry(v,z,vel,force,pot,&f,&fd,&R,&hd);
    double a=exp(v[15]),a2=a*a,a3=a*a2;
    double complex A[3][3]={{0}},B[3][3]={{0}};
    A[0][0]=4*lx*pow(cabs(z[0]),2)+xi*R;
    B[0][0]=2*lx*z[0]*z[0]+6*kap*z[1]*conj(z[0]);
    A[0][1]=3*kap*conj(z[0])*conj(z[0]);A[1][0]=conj(A[0][1]);
    A[1][1]=4*ly*pow(cabs(z[1]),2)+mm*mm;
    B[1][1]=2*ly*z[1]*z[1]+2*mu*z[2];
    A[1][2]=2*mu*conj(z[1]);A[2][1]=conj(A[1][2]);
    A[2][2]=4*ls*pow(cabs(z[2]),2)-ls*hh*hh;B[2][2]=2*ls*z[2]*z[2];
    for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
        mat[2*i][2*j]=a2*creal(A[i][j]+B[i][j]);
        mat[2*i][2*j+1]=-a2*cimag(A[i][j]-B[i][j]);
        mat[2*i+1][2*j]=a2*cimag(A[i][j]+B[i][j]);
        mat[2*i+1][2*j+1]=a2*creal(A[i][j]-B[i][j]);
    }
    const double charge[]={1,3,6};
    for(int i=0;i<6;i++)mat[i][i]-=a2*R/6;
    for(int i=0;i<3;i++) {
        double complex w=sqrt(2.)*a*z[i],wp=sqrt(2.)*a2*(vel[i]+v[16]*z[i]);
        g[2*i]=-charge[i]*cimag(w);g[2*i+1]=charge[i]*creal(w);
        gp[2*i]=-charge[i]*cimag(wp);gp[2*i+1]=charge[i]*creal(wp);
        double extra=a2*((i==0?xi:0)-1./6)*R;
        ward[2*i]=-sqrt(2.)*charge[i]*a3*cimag(force[i])+extra*g[2*i];
        ward[2*i+1]=sqrt(2.)*charge[i]*a3*creal(force[i])+extra*g[2*i+1];
    }
}

static void full_rhs(double t,const double *v,double *out) {
    rhs(t,v,out);
#ifdef NEUTRAL_SECTOR
    double complex z[3],vel[3],force[3];double pot,f,fd,R,hd;
    matter(v,z,vel,force,&pot);geometry(v,z,vel,force,pot,&f,&fd,&R,&hd);
    double a=exp(v[15]);
    /* xi_phi=0, g_X=.1, lambda_phi=.03. The vacuum scale/H_i is tiny. */
    double vacuum_ratio=.924089637375362/133349619398.33566;
    double mass=a*a*(.02*pow(cabs(z[0]),2)-.03*vacuum_ratio*vacuum_ratio-R/6);
    for(int m=0;m<nk;m++) {
        const double *T=v+19+BLOCK*m;double *d=out+19+BLOCK*m;
        for(int j=0;j<2;j++){d[j]=T[2+j]/a;d[2+j]=-(waves[m]*waves[m]+mass)*T[j]/a;}
    }
#else
    double M[6][6],g[6],gp[6],ward[6],gg=0,a=exp(v[15]);
    perturbation_matrix(v,M,g,gp,ward);for(int i=0;i<6;i++)gg+=g[i]*g[i];
    for(int m=0;m<nk;m++) {
        const double *T=v+19+144*m;double *d=out+19+144*m;
        double k2=waves[m]*waves[m],c=eg?eg*eg/k2:0;
        for(int col=0;col<12;col++) {
            double current=0,common=0;
            for(int i=0;i<6;i++){current+=g[i]*T[(i+6)*12+col]-gp[i]*T[i*12+col];common+=g[i]*T[i*12+col];}
            common/=gg;
            for(int i=0;i<6;i++) {
                double mq=common*ward[i];
                for(int j=0;j<6;j++)mq+=M[i][j]*(T[j*12+col]-common*g[j]);
                d[i*12+col]=(T[(i+6)*12+col]+c*g[i]*current)/a;
                d[(i+6)*12+col]=(-k2*T[i*12+col]-mq+c*gp[i]*current)/a;
            }
        }
    }
#endif
}

int main(int argc,char **argv) {
    if(argc!=3){fprintf(stderr,"modes input.bin output.bin\n");return 1;}
    FILE *in=fopen(argv[1],"rb"),*out=fopen(argv[2],"wb");if(!in||!out)return 2;
    int32_t header[2];if(fread(header,sizeof(int32_t),2,in)!=2)return 3;
    nk=header[0];int nt=header[1];double pars[9];
    if(nk<1||nk>256||nt<2||fread(pars,sizeof(double),9,in)!=9)return 4;
    beta=pars[0];decay=pars[1];eg=pars[2];xi=pars[3];eta2=pars[4];mu=pars[5];kap=pars[6];double end=pars[7],rtol=pars[8],atol=rtol*.001;
    nvar=19+BLOCK*nk;double *v=calloc(nvar,sizeof(double)),*tmp=calloc(nvar,sizeof(double)),*next=calloc(nvar,sizeof(double));
    double *stages=calloc(13*nvar,sizeof(double)),*times=calloc(nt,sizeof(double));waves=calloc(nk,sizeof(double));
    if(fread(v,sizeof(double),19,in)!=19||fread(waves,sizeof(double),nk,in)!=(size_t)nk||fread(times,sizeof(double),nt,in)!=(size_t)nt)return 5;
    fclose(in);for(int m=0;m<nk;m++)for(int i=0;i<NMODE;i++)v[19+BLOCK*m+NMODE*i+i]=1;
    double t=0,step=.0001,maxf=0,maxg=0;long count=0,reject=0;full_rhs(t,v,stages);
    for(int sample=0;sample<nt;sample++) {
        double target=times[sample];
        while(t<target) {
            double dt=fmin(step,target-t);if(dt<1e-14)return 6;
            for(int st=1;st<12;st++) {
                for(int j=0;j<nvar;j++){double sum=0;for(int k=0;k<st;k++)sum+=AA[st][k]*stages[k*nvar+j];tmp[j]=v[j]+dt*sum;}
                full_rhs(t+CC[st]*dt,tmp,stages+st*nvar);
            }
            for(int j=0;j<nvar;j++){double sum=0;for(int k=0;k<12;k++)sum+=BB[k]*stages[k*nvar+j];next[j]=v[j]+dt*sum;}
            full_rhs(t+dt,next,stages+12*nvar);double e3sq=0,e5sq=0;
            for(int j=0;j<nvar;j++) {
                double e3=0,e5=0;for(int k=0;k<13;k++){e3+=EE3[k]*stages[k*nvar+j];e5+=EE5[k]*stages[k*nvar+j];}
                double scale=atol+rtol*fmax(fabs(v[j]),fabs(next[j]));e3sq+=pow(e3/scale,2);e5sq+=pow(e5/scale,2);
            }
            double err=(e3sq==0&&e5sq==0)?0:fabs(dt)*e5sq/sqrt((e5sq+.01*e3sq)*nvar);
            if(!isfinite(err)){fprintf(stderr,"nonfinite %.17g\n",t);return 7;}
            double fac=err>0?.9*pow(err,-.125):5.;fac=fmax(.2,fmin(5.,fac));
            if(err<=1) {
                t+=dt;memcpy(v,next,nvar*sizeof(double));memcpy(stages,stages+12*nvar,nvar*sizeof(double));count++;
                double q[NO];quantities(t,v,q);maxf=fmax(maxf,fabs(q[15]));maxg=fmax(maxg,fabs(q[8]));
            }else reject++;
            step=fmin(.08/fmax(1.,beta),dt*fac);if(count>50000000L)return 8;
        }
        if(fwrite(&t,sizeof(double),1,out)!=1||fwrite(v,sizeof(double),nvar,out)!=(size_t)nvar)return 9;
    }
    fclose(out);printf("CHECK %ld %ld %.17g %.17g %.17g\n",count,reject,maxf,maxg,end);
    free(v);free(tmp);free(next);free(stages);free(waves);free(times);return 0;
}
