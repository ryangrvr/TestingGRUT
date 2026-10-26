/* Original coupled common-mode extension of own DOP853 engine. */
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

#define NC 3
#define ND (13+24*NC)
static double gauge_e,wavenumber[NC];
static double complex initial_ratio[NC];
#define NO 27
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
    /* Additional original diagnostics, following the mass-eigenmode
     * dictionary credited to Claude f1aa7df. Evolution equations unchanged. */
    o[20]=-1.;o[21]=-1.;
    double sr=cabs(u[2]);
    double phaseprime=sr>1e-30?cimag(conj(u[2])*p[2])/(sr*sr):0.;
    if(D>0 && sr>1e-30) {
        double complex J=kap*u[0]*u[0]*u[0];
        double complex yp=(-A*J+B*conj(J))/D;
        double complex JD=3*kap*u[0]*u[0]*p[0]+2*mm*mm*a*yp+2*mu*(u[2]+a*p[2])*conj(yp);
        double complex ypd=(-A*JD+B*conj(JD))/D;
        double complex rotation=cexp(-.5*I*carg(u[2]));
        double complex d=(u[1]-yp)*rotation;
        double complex dd=(p[1]-ypd-(u[1]-yp)/a-.5*I*phaseprime*(u[1]-yp))*rotation;
        double wh=sqrt(A+cabs(B)),wl=sqrt(A-cabs(B));
        o[20]=(pow(creal(dd),2)+wh*wh*pow(creal(d),2))/wh;
        o[21]=(pow(cimag(dd),2)+wl*wl*pow(cimag(d),2))/wl;
    }
    double rd=sr>1e-30?creal(conj(u[2])*p[2])/sr-sr/a:0.;
    double Erad=rd*rd+ls*pow(sr*sr-a*a*hh*hh/2,2);
    o[22]=Erad/(a*sqrt(2*ls)*hh);o[23]=Erad;
    o[24]=sr>1e-30?pow(cimag(conj(u[2])*p[2]),2)/(sr*sr):0.;
    o[25]=phaseprime;o[26]=2*ly*y2/A;

}


/* Coulomb-gauge scalar/longitudinal Hamiltonian, original derivation.
 * dq = g + w, canonical momentum = g' + r for one initial common mode.
 * H2=1/2 p^2+1/2 q^T(k^2+M6)q+e^2/(2k^2)(g.p-g'.q)^2.
 */
static void gauge_vectors(const double *v,double *g,double *gp) {
    const double charges[]={1,3,6};
    for(int i=0;i<3;i++) {
        g[2*i]=-sqrt(2.)*charges[i]*v[2*i+1];
        g[2*i+1]=sqrt(2.)*charges[i]*v[2*i];
        gp[2*i]=-sqrt(2.)*charges[i]*v[7+2*i];
        gp[2*i+1]=sqrt(2.)*charges[i]*v[6+2*i];
    }
}
static void projection(double a,const double *v) {
    double g[6],gp[6];gauge_vectors(v,g,gp);double gg=0,ggp=0;
    for(int i=0;i<6;i++){gg+=g[i]*g[i];ggp+=g[i]*gp[i];}
    printf("P %.17g",a);
    for(int j=0;j<NC;j++) {
        int off=13+24*j;double complex dw[6],dr[6],dtheta=0,jbar=0,numerator=0;
        for(int i=0;i<6;i++) {
            dw[i]=v[off+2*i]+I*v[off+2*i+1];
            dr[i]=v[off+12+2*i]+I*v[off+13+2*i];
            dtheta+=g[i]*dw[i];jbar+=g[i]*dr[i]-gp[i]*dw[i];
            numerator+=gp[i]*dw[i]+g[i]*dr[i];
        }
        dtheta/=gg;
        double coupling=gauge_e*gauge_e/(wavenumber[j]*wavenumber[j]);
        double complex theta=1+dtheta;
        double complex thetap=(numerator+coupling*gg*jbar-2*dtheta*ggp)/gg;
        printf(" %.17g %.17g %.17g %.17g",creal(theta),cimag(theta),creal(thetap),cimag(thetap));
    }
    printf("\n");
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

    double complex A[3][3]={{0}},B[3][3]={{0}};double M[6][6];
    A[0][0]=4*lx*x2;B[0][0]=2*lx*X*X+6*kap*Y*conj(X);
    A[0][1]=3*kap*conj(X)*conj(X);A[1][0]=conj(A[0][1]);
    A[1][1]=4*ly*y2+mm*mm*a*a;B[1][1]=2*ly*Y*Y+2*mu*a*S;
    A[1][2]=2*mu*a*conj(Y);A[2][1]=conj(A[1][2]);
    A[2][2]=4*ls*s2-ls*a*a*hh*hh;B[2][2]=2*ls*S*S;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
        M[2*i][2*j]=creal(A[i][j]+B[i][j]);
        M[2*i][2*j+1]=-cimag(A[i][j]-B[i][j]);
        M[2*i+1][2*j]=cimag(A[i][j]+B[i][j]);
        M[2*i+1][2*j+1]=creal(A[i][j]-B[i][j]);
    }
    double g[6],gp[6],tf[6],gg=0;gauge_vectors(v,g,gp);
    const double charge[]={1,3,6};
    for(int i=0;i<3;i++){tf[2*i]=-sqrt(2.)*charge[i]*cimag(F[i]);tf[2*i+1]=sqrt(2.)*charge[i]*creal(F[i]);}
    for(int i=0;i<6;i++)gg+=g[i]*g[i];
    for(int j=0;j<NC;j++) {
        int off=13+24*j;double complex w[6],r[6],jbar=0;
        double wave2=wavenumber[j]*wavenumber[j],coupling=gauge_e*gauge_e/wave2;
        for(int i=0;i<6;i++) {
            w[i]=v[off+2*i]+I*v[off+2*i+1];
            r[i]=v[off+12+2*i]+I*v[off+13+2*i];
            jbar+=g[i]*r[i]-gp[i]*w[i];
        }
        double complex common=0;
        for(int i=0;i<6;i++)common+=g[i]*w[i];
        common/=gg;
        for(int i=0;i<6;i++) {
            double complex wd=r[i]+coupling*g[i]*jbar;
            double complex rd=-wave2*(w[i]+g[i])+coupling*gp[i]*jbar-common*tf[i];
            /* M g = T force: use the exact gauge Ward identity for the
             * common component instead of subtracting large Hessian terms. */
            for(int q=0;q<6;q++)rd-=M[i][q]*(w[q]-common*g[q]);
            f[off+2*i]=creal(wd);f[off+2*i+1]=cimag(wd);
            f[off+12+2*i]=creal(rd);f[off+13+2*i]=cimag(rd);
        }
    }

}

int main(int argc,char **argv) {
    if(argc!=25){fprintf(stderr,"release M h mu kappa end rtol atol x_re x_im y_re y_im s_re s_im nsamp\n");return 1;}
    mm=atof(argv[1]);hh=atof(argv[2]);mu=atof(argv[3]);kap=atof(argv[4]);
    double end=atof(argv[5]),rtol=atof(argv[6]),atol=atof(argv[7]);int ns=atoi(argv[14]);
    double v[ND]={0};for(int j=0;j<6;j++)v[j]=v[j+6]=atof(argv[j+8]);
    gauge_e=atof(argv[15]);
    for(int j=0;j<NC;j++) {wavenumber[j]=atof(argv[16+3*j]);initial_ratio[j]=atof(argv[17+3*j])+I*atof(argv[18+3*j]);}
    double ginit[6],gpinit[6],gginit=0;gauge_vectors(v,ginit,gpinit);
    for(int i=0;i<6;i++)gginit+=ginit[i]*ginit[i];
    for(int j=0;j<NC;j++) {
        double wave2=wavenumber[j]*wavenumber[j];
        for(int i=0;i<6;i++) {
            double complex r=ginit[i]*wave2/(wave2+gauge_e*gauge_e*gginit)*initial_ratio[j];
            v[13+24*j+12+2*i]=creal(r);v[13+24*j+13+2*i]=cimag(r);
        }
    }
    projection(1.,v);
    double q[NO],old[NO];quantities(1.,v,q);v[12]=q[15];
    double E0=q[15],N0=q[10],maxG=0,maxW=0;
    double stepsize=.001,a=1.;long steps=0,rejects=0;
    double ks[13][ND],temp[ND],next[ND];rhs(a,v,ks[0]);
    const double edges[]={1,5,10,20,30,40,60,100};
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
        for(int j=0;j<ND;j++) {
            if(j==12)continue;
            double e3=0,e5=0;
            for(int k=0;k<13;k++){e3+=EE3[k]*ks[k][j];e5+=EE5[k]*ks[k][j];}
            double sc=atol+rtol*fmax(fabs(v[j]),fabs(next[j]));
            e5sq+=pow(e5/sc,2);e3sq+=pow(e3/sc,2);
        }
        double err=(e5sq==0&&e3sq==0)?0:fabs(h)*e5sq/sqrt((e5sq+.01*e3sq)*(ND-1));
        if(!isfinite(err)){fprintf(stderr,"nonfinite error at %g\n",a);return 3;}
        double factor=err>0?.9*pow(err,-.125):5.;
        factor=fmin(5.,fmax(.2,factor));
        if(err<=1) {
            for(int j=0;j<NO;j++)old[j]=q[j];
            double aa=a;a+=h;
            for(int j=0;j<ND;j++){v[j]=next[j];ks[0][j]=ks[12][j];}
            quantities(a,v,q);steps++;
            if(fabs(a-10)<1e-12||fabs(a-20)<1e-12||fabs(a-30)<1e-12||fabs(a-40)<1e-12)projection(a,v);
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
