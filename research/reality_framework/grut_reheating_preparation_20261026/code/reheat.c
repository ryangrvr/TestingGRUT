/* Jordan-frame homogeneous nonminimal X/Y/Sigma + quadratic inflaton + RD.
 * Author DOP853 implementation with SciPy/Hairer coefficient tables.
 * No spatial fluctuation, hidden decay or quantum backreaction certification.
 */
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "dop853.h"
#define ND 19
#define NO 29
static double beta,decay,xi,eta2;
static double lx=.001,ly=1,ls=1,mm=30,hh=40,mu,kap;

static void matter(const double *v,double complex *z,double complex *vel,double complex *force,double *pot) {
    for(int i=0;i<3;i++){z[i]=v[2+2*i]+I*v[3+2*i];vel[i]=v[8+2*i]+I*v[9+2*i];}
    double x2=pow(cabs(z[0]),2),y2=pow(cabs(z[1]),2),s2=pow(cabs(z[2]),2);
    force[0]=2*lx*x2*z[0]+3*kap*z[1]*conj(z[0])*conj(z[0]);
    force[1]=(mm*mm+2*ly*y2)*z[1]+kap*z[0]*z[0]*z[0]+2*mu*z[2]*conj(z[1]);
    force[2]=2*ls*(s2-hh*hh/2)*z[2]+mu*z[1]*z[1];
    *pot=lx*x2*x2+ly*y2*y2+mm*mm*y2+ls*pow(s2-hh*hh/2,2)
          +2*creal(kap*conj(z[1])*z[0]*z[0]*z[0]+mu*conj(z[2])*z[1]*z[1]);
}

static void geometry(const double *v,const double complex *z,const double complex *vel,const double complex *force,
                     double pot,double *ff,double *fd,double *ricci,double *hdot) {
    double K=0;for(int i=0;i<3;i++)K+=pow(cabs(vel[i]),2);
    double x2=pow(cabs(z[0]),2),vx2=pow(cabs(vel[0]),2);
    *ff=1-2*xi*eta2*x2;
    *fd=-4*xi*eta2*creal(conj(z[0])*vel[0]);
    double numerator=2*beta*beta*v[0]*v[0]-v[1]*v[1]
      +eta2*(4*pot-2*K+12*xi*vx2-12*xi*creal(conj(z[0])*force[0]));
    *ricci=numerator/(*ff+12*xi*xi*eta2*x2);
    double complex ax=-3*v[16]*vel[0]-force[0]-xi*(*ricci)*z[0];
    double fdd=-4*xi*eta2*(vx2+creal(conj(z[0])*ax));
    *hdot=-(v[1]*v[1]+2*eta2*K+4*v[14]/3+fdd-v[16]*(*fd))/(2*(*ff));
}

static void rhs(double t,const double *v,double *out) {
    (void)t;double complex z[3],vel[3],force[3];double pot,ff,fd,ricci,hdot;
    matter(v,z,vel,force,&pot);geometry(v,z,vel,force,pot,&ff,&fd,&ricci,&hdot);
    out[0]=v[1];out[1]=-(3*v[16]+decay)*v[1]-beta*beta*v[0];
    double K=0;
    for(int i=0;i<3;i++){
        double complex acc=-3*v[16]*vel[i]-force[i]-(i==0?xi*ricci*z[0]:0);
        out[2+2*i]=creal(vel[i]);out[3+2*i]=cimag(vel[i]);
        out[8+2*i]=creal(acc);out[9+2*i]=cimag(acc);K+=pow(cabs(vel[i]),2);
    }
    out[14]=decay*v[1]*v[1]-4*v[16]*v[14];out[15]=v[16];out[16]=hdot;
    out[17]=-3*v[16]*(v[1]*v[1]+2*eta2*K+4*v[14]/3)+.5*fd*ricci;
    out[18]=-6*v[16]*K-2*xi*ricci*creal(conj(z[0])*vel[0]);
}

static void quantities(double t,const double *v,double *q) {
    double complex z[3],vel[3],force[3];double pot,ff,fd,ricci,hdot;
    matter(v,z,vel,force,&pot);geometry(v,z,vel,force,pot,&ff,&fd,&ricci,&hdot);
    double a=exp(v[15]),a3=a*a*a,a4=a3*a,K=0;
    double n[3];for(int i=0;i<3;i++){n[i]=-2*cimag(conj(z[i])*vel[i])*a3;K+=pow(cabs(vel[i]),2);}
    double rho=.5*v[1]*v[1]+.5*beta*beta*v[0]*v[0]+eta2*(K+pot)+v[14];
    double ex=a4*(pow(cabs(vel[0]),2)+lx*pow(cabs(z[0]),4));
    double nmax=2*sqrt(2*lx)*pow(ex/(3*lx),.75);
    double A=mm*mm;double complex B=2*mu*z[2],J=kap*z[0]*z[0]*z[0];
    double D=A*A-pow(cabs(B),2),nh=-1,nl=-1;
    if(D>0 && cabs(z[2])>1e-20){
        double complex yp=(-A*J+B*conj(J))/D;
        double complex jd=3*kap*z[0]*z[0]*vel[0]+2*mu*vel[2]*conj(yp);
        double complex ypd=(-A*jd+B*conj(jd))/D;
        double complex rot=cexp(-.5*I*carg(z[2]));
        double theta=cimag(conj(z[2])*vel[2])/pow(cabs(z[2]),2);
        double complex dy=(z[1]-yp)*rot,dv=(vel[1]-ypd-.5*I*theta*(z[1]-yp))*rot;
        double wh=sqrt(A+cabs(B)),wl=sqrt(A-cabs(B));
        nh=a3*(pow(creal(dv),2)+wh*wh*pow(creal(dy),2))/wh;
        nl=a3*(pow(cimag(dv),2)+wl*wl*pow(cimag(dy),2))/wl;
    }
    double sr=cabs(z[2]);double sd=sr>1e-20?creal(conj(z[2])*vel[2])/sr:0;
    double ns=a3*(sd*sd+ls*pow(sr*sr-hh*hh/2,2))/(sqrt(2*ls)*hh);
    double vi=.5*beta*beta*v[0]*v[0],ki=.5*v[1]*v[1];
    q[0]=t;q[1]=a;q[2]=v[16];q[3]=ricci;q[4]=(ki-vi)/fmax(ki+vi,1e-300);
    q[5]=v[14]/fmax(rho,1e-300);q[6]=n[0];q[7]=n[0]+3*n[1];q[8]=n[0]+3*n[1]+6*n[2];
    q[9]=ex;q[10]=n[0]/nmax;q[11]=nh;q[12]=nl;q[13]=ns;q[14]=sr/(hh/sqrt(2));
    q[15]=(3*ff*v[16]*v[16]+3*v[16]*fd-rho)/fmax(rho,1e-300);
    q[16]=rho-v[17];q[17]=K+pot-v[18];q[18]=ff;
    q[19]=ricci-6*(2*v[16]*v[16]+hdot);
    q[20]=xi*ricci;
    q[21]=4*lx*pow(cabs(z[0]),2)+xi*ricci-cabs(2*lx*z[0]*z[0]+6*kap*z[1]*conj(z[0]));
    q[22]=v[14];q[23]=rho;q[24]=2*creal(kap*conj(z[1])*z[0]*z[0]*z[0]);
    q[25]=2*creal(mu*conj(z[2])*z[1]*z[1]);q[26]=cabs(z[0]);q[27]=cabs(z[1]);q[28]=K+pot;
}

int main(int argc,char **argv){
    if(argc!=19){fprintf(stderr,"reheat beta Gamma xi eta2 mu kappa end rtol atol Xre Xim Yre Yim Sre Sim I0 rad0 samples\n");return 1;}
    beta=atof(argv[1]);decay=atof(argv[2]);xi=atof(argv[3]);eta2=atof(argv[4]);mu=atof(argv[5]);kap=atof(argv[6]);
    double end=atof(argv[7]),rtol=atof(argv[8]),atol=atof(argv[9]);
    double v[ND]={0};for(int i=0;i<6;i++)v[2+i]=atof(argv[10+i]);
    v[0]=atof(argv[16]);v[14]=atof(argv[17]);v[16]=1.;int samples=atoi(argv[18]);
    double q[NO],E0,Eh0,Qscale=1;quantities(0,v,q);E0=q[23];Eh0=q[28];
    Qscale=fmax(1.,2*sqrt(2*lx)*pow(q[9]/(3*lx),.75));
    double t=0,step=.0001;long nstep=0,nreject=0;
    double ks[13][ND],tmp[ND],next[ND];rhs(0,v,ks[0]);
    double maxf=0,maxg=0,maxt=0,maxe=0,maxh=0,minff=q[18];
    printf("S");for(int j=0;j<NO;j++)printf(" %.17g",q[j]);for(int j=0;j<ND;j++)printf(" %.17g",v[j]);printf("\n");
    for(int sample=1;sample<=samples;sample++){
        double target=end*sample/samples;
        while(t<target){
            double dt=fmin(step,target-t);if(dt<1e-14){fprintf(stderr,"step too small %g\n",t);return 2;}
            for(int st=1;st<12;st++){
                for(int j=0;j<ND;j++){double sum=0;for(int k=0;k<st;k++)sum+=AA[st][k]*ks[k][j];tmp[j]=v[j]+dt*sum;}
                rhs(t+CC[st]*dt,tmp,ks[st]);
            }
            for(int j=0;j<ND;j++){double sum=0;for(int k=0;k<12;k++)sum+=BB[k]*ks[k][j];next[j]=v[j]+dt*sum;}
            rhs(t+dt,next,ks[12]);double e3sq=0,e5sq=0;
            for(int j=0;j<ND;j++){
                double e3=0,e5=0;for(int k=0;k<13;k++){e3+=EE3[k]*ks[k][j];e5+=EE5[k]*ks[k][j];}
                double sc=atol+rtol*fmax(fabs(v[j]),fabs(next[j]));e3sq+=pow(e3/sc,2);e5sq+=pow(e5/sc,2);
            }
            double err=(e3sq==0&&e5sq==0)?0:fabs(dt)*e5sq/sqrt((e5sq+.01*e3sq)*ND);
            if(!isfinite(err)){fprintf(stderr,"nonfinite at %g\n",t);return 3;}
            double factor=err>0?.9*pow(err,-.125):5.;factor=fmax(.2,fmin(5.,factor));
            if(err<=1){
                t+=dt;for(int j=0;j<ND;j++){v[j]=next[j];ks[0][j]=ks[12][j];}nstep++;
                quantities(t,v,q);maxf=fmax(maxf,fabs(q[15]));maxg=fmax(maxg,fabs(q[8])/Qscale);
                maxt=fmax(maxt,fabs(q[19])/fmax(1.,fabs(q[3])));
                maxe=fmax(maxe,fabs(q[16]/E0-1));maxh=fmax(maxh,fabs(q[17]/Eh0-1));minff=fmin(minff,q[18]);
            }else nreject++;
            step=fmin(.08/fmax(1.,beta),dt*factor);
            if(nstep>50000000L){fprintf(stderr,"step limit\n");return 4;}
        }
        quantities(t,v,q);printf("S");for(int j=0;j<NO;j++)printf(" %.17g",q[j]);for(int j=0;j<ND;j++)printf(" %.17g",v[j]);printf("\n");
    }
    printf("CHECK %ld %ld %.17g %.17g %.17g %.17g %.17g %.17g\n",nstep,nreject,maxf,maxg,maxt,maxe,maxh,minff);
    return 0;
}
