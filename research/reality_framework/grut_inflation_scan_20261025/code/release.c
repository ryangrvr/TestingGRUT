/* Full X/Y/Sigma release in radiation: a=conformal time in H_i units.
 * Original Codex C implementation of DOP853 (Hairer; SciPy coefficient table).
 * State is (Re/Im chi_i, Re/Im pi_i, integrated explicit Hamiltonian work).
 * chi=a z/H_i; pi=dchi/da. No physical initial velocity means pi_i=chi_i.
 * The one-field EFT is never used in these evolution equations.
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "dop853.h"

#define ND 13
#define NO 20
static double lx=1e-3,ly=1.,ls=1.,kap=1e-3,mm,hh,mu;

static void quantities(double a,const double *v,double *o) {
    double complex u[3],p[3],w[3];double K[3],n[3];
    for(int i=0;i<3;i++) {
        u[i]=v[2*i]+I*v[2*i+1];p[i]=v[6+2*i]+I*v[7+2*i];
        w[i]=p[i]-u[i]/a;K[i]=creal(w[i]*conj(w[i]));
        n[i]=-2*cimag(conj(u[i])*p[i]);
    }
    double x2=creal(u[0]*conj(u[0])),y2=creal(u[1]*conj(u[1])),s2=creal(u[2]*conj(u[2]));
    double EX=K[0]+lx*x2*x2;
    double EY=K[1]+mm*mm*a*a*y2+ly*y2*y2;
    double ES=K[2]+ls*pow(s2-a*a*hh*hh/2,2);
    double VK=2*creal(kap*conj(u[1])*u[0]*u[0]*u[0]);
    double VM=2*creal(mu*a*conj(u[2])*u[1]*u[1]);
    double VC=lx*x2*x2+ly*y2*y2+ls*pow(s2-a*a*hh*hh/2,2)+mm*mm*a*a*y2+VK+VM;
    double EC=VC;
    for(int i=0;i<3;i++) EC+=creal(p[i]*conj(p[i]));
    double nmax=2*sqrt(2*lx)*pow(EX/(3*lx),.75);
    double fgchi=sqrt(2*(x2+9*y2+36*s2));
    double A=mm*mm*a*a;double complex B=2*mu*a*u[2];
    double D=A*A-creal(B*conj(B));
    double EYexc=-1.,JY=-1.;
    if(D>0) {
        double complex J=kap*u[0]*u[0]*u[0];
        double complex yp=(-A*J+B*conj(J))/D;
        double complex JD=3*kap*u[0]*u[0]*p[0]+2*mm*mm*a*yp+2*mu*(u[2]+a*p[2])*conj(yp);
        double complex ypd=(-A*JD+B*conj(JD))/D;
        double complex dy=u[1]-yp,dw=p[1]-ypd-dy/a;
        EYexc=creal(dw*conj(dw))+A*creal(dy*conj(dy))+creal(conj(B)*dy*dy);
        JY=EYexc/(mm*a);
    }
    o[0]=n[0];o[1]=n[1];o[2]=n[2];o[3]=n[0]+3*n[1];
    o[4]=EX;o[5]=EY;o[6]=ES;o[7]=VK;o[8]=VM;o[9]=n[0]/nmax;
    o[10]=nmax;o[11]=EYexc;o[12]=JY;o[13]=fgchi/a;
    o[14]=1e-7*fgchi/a;o[15]=EC;o[16]=EX+EY+ES+VK+VM;
    o[17]=n[0]+3*n[1]+6*n[2];o[18]=EC-v[12];o[19]=cabs(u[1]);
}

static void rhs(double a,const double *v,double *f) {
    double complex X=v[0]+I*v[1],Y=v[2]+I*v[3],S=v[4]+I*v[5];
    double x2=creal(X*conj(X)),y2=creal(Y*conj(Y)),s2=creal(S*conj(S));
    double complex F[3]={2*lx*x2*X+3*kap*Y*conj(X)*conj(X),
        2*ly*y2*Y+mm*mm*a*a*Y+kap*X*X*X+2*mu*a*S*conj(Y),
        2*ls*(s2-a*a*hh*hh/2)*S+mu*a*Y*Y};
    for(int i=0;i<3;i++) {
        f[2*i]=v[6+2*i];f[2*i+1]=v[7+2*i];
        f[6+2*i]=-creal(F[i]);f[7+2*i]=-cimag(F[i]);
    }
    f[12]=2*mm*mm*a*y2-2*ls*a*hh*hh*(s2-a*a*hh*hh/2)+2*mu*creal(conj(S)*Y*Y);
}

int main(int argc,char **argv) {
    if(argc!=14){fprintf(stderr,"release M h mu end rtol atol x_re x_im y_re y_im s_re s_im nsamp\n");return 1;}
    mm=atof(argv[1]);hh=atof(argv[2]);mu=atof(argv[3]);
    double end=atof(argv[4]),rtol=atof(argv[5]),atol=atof(argv[6]);int ns=atoi(argv[13]);
    double v[ND]={0};for(int j=0;j<6;j++)v[j]=v[j+6]=atof(argv[j+7]);
    double q[NO],old[NO];quantities(1.,v,q);v[12]=q[15];
    double E0=q[15],N0=q[10],maxG=0,maxW=0;
    double stepsize=.001,a=1.;long steps=0,rejects=0;
    double ks[13][ND],temp[ND],next[ND];rhs(a,v,ks[0]);
    const double edges[]={1,100,300,600,1000,2000,3000,10000};
    double sums[7][NO]={{0}},sq[7][3]={{0}},span[7]={0};int wi=0,si=0;
    double nextsample=exp(log(end)/ns);
    printf("S %.17g",a);for(int j=0;j<NO;j++)printf(" %.17g",q[j]);printf("\n");
    while(a<end) {
        double target=fmin(end,nextsample);
        if(wi<7)target=fmin(target,edges[wi+1]);
        double h=fmin(stepsize,target-a);
        if(h<=0){fprintf(stderr,"invalid step %.17g %.17g\n",a,target);return 2;}
        for(int stage=1;stage<12;stage++) {
            for(int j=0;j<ND;j++) {
                double sum=0;for(int k=0;k<stage;k++)sum+=AA[stage][k]*ks[k][j];
                temp[j]=v[j]+h*sum;
            }
            rhs(a+CC[stage]*h,temp,ks[stage]);
        }
        for(int j=0;j<ND;j++) {
            double sum=0;for(int k=0;k<12;k++)sum+=BB[k]*ks[k][j];
            next[j]=v[j]+h*sum;
        }
        rhs(a+h,next,ks[12]);
        double e5sq=0,e3sq=0;
        for(int j=0;j<12;j++) {
            double e3=0,e5=0;
            for(int k=0;k<13;k++){e3+=EE3[k]*ks[k][j];e5+=EE5[k]*ks[k][j];}
            double sc=atol+rtol*fmax(fabs(v[j]),fabs(next[j]));
            e5sq+=pow(e5/sc,2);e3sq+=pow(e3/sc,2);
        }
        double err=(e5sq==0&&e3sq==0)?0:fabs(h)*e5sq/sqrt((e5sq+.01*e3sq)*12);
        if(!isfinite(err)){fprintf(stderr,"nonfinite error at %g\n",a);return 3;}
        double factor=err>0?.9*pow(err,-.125):5.;
        factor=fmin(5.,fmax(.2,factor));
        if(err<=1) {
            for(int j=0;j<NO;j++)old[j]=q[j];
            double aa=a;a+=h;
            for(int j=0;j<ND;j++){v[j]=next[j];ks[0][j]=ks[12][j];}
            quantities(a,v,q);steps++;
            maxG=fmax(maxG,fabs(q[17])/N0);maxW=fmax(maxW,fabs(q[18])/E0);
            if(wi<7) {
                span[wi]+=h;
                for(int j=0;j<NO;j++)sums[wi][j]+=.5*h*(old[j]+q[j]);
                int jj[3]={0,3,9};
                for(int k=0;k<3;k++)sq[wi][k]+=.5*h*(old[jj[k]]*old[jj[k]]+q[jj[k]]*q[jj[k]]);
                if(a>=edges[wi+1]*(1-1e-14))wi++;
            }
            if(a>=nextsample*(1-1e-14)||a>=end) {
                printf("S %.17g",a);for(int j=0;j<NO;j++)printf(" %.17g",q[j]);printf("\n");
                si++;nextsample=si<ns?exp(log(end)*(si+1)/ns):2*end;
            }
            (void)aa;
        } else rejects++;
        stepsize=fmin(.3,h*factor);
        if(steps>1500000000L||stepsize<1e-14){fprintf(stderr,"step limit at %g\n",a);return 4;}
    }
    for(int k=0;k<7;k++)if(span[k]>0){
        printf("W %.17g %.17g %.17g",edges[k],fmin(edges[k+1],end),span[k]);
        for(int j=0;j<NO;j++)printf(" %.17g",sums[k][j]/span[k]);
        for(int j=0;j<3;j++)printf(" %.17g",sq[k][j]/span[k]);
        printf("\n");
    }
    printf("F");for(int j=0;j<12;j++)printf(" %.17g",v[j]);printf("\n");
    printf("END %ld %ld %.17g %.17g\n",steps,rejects,maxG,maxW);
    return 0;
}
