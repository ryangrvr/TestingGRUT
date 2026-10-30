/* Six real scalars plus eliminated A0 on a prescribed inflaton/radiation
 * background. Metric, inflaton and radiation perturbations are excluded.
 * Canonical 12x12 transfer follows the inherited Coulomb-gauge Hamiltonian.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "values.h"
#include "dop853.h"
static double cr,ci,eg;static int nk,nv;static double *wave;
static void charge_rotate(const double *r,double *g){const double qs[]={1,3,6};for(int i=0;i<3;i++){g[2*i]=-qs[i]*r[2*i+1];g[2*i+1]=qs[i]*r[2*i];}}
static void matter(const double *r,double R,double *V,double *F,double *VR,double D[6][6]){
 double complex x=(r[0]+I*r[1])/sqrt(2.),y=(r[2]+I*r[3])/sqrt(2.),s=(r[4]+I*r[5])/sqrt(2.),lock=cr+I*ci,u=MU+lock*R;
 double xx=creal(x*conj(x)),yy=creal(y*conj(y)),ss=creal(s*conj(s));
 *V=LX*xx*xx+LY*yy*yy+MM*MM*yy+LS*pow(ss-HH*HH/2,2)+2*creal(KAP*conj(y)*x*x*x+u*conj(s)*y*y)+XI*R*xx;
 *VR=XI*xx+2*creal(lock*conj(s)*y*y);
 double complex ff[]={2*LX*xx*x+3*KAP*y*conj(x)*conj(x)+XI*R*x,(MM*MM+2*LY*yy)*y+KAP*x*x*x+2*conj(u)*s*conj(y),2*LS*(ss-HH*HH/2)*s+u*y*y};
 for(int i=0;i<3;i++){F[2*i]=sqrt(2.)*creal(ff[i]);F[2*i+1]=sqrt(2.)*cimag(ff[i]);}
 double complex A[3][3]={{0}},B[3][3]={{0}};
 A[0][0]=4*LX*xx+XI*R;B[0][0]=2*LX*x*x+6*KAP*y*conj(x);
 A[0][1]=3*KAP*conj(x)*conj(x);A[1][0]=conj(A[0][1]);
 A[1][1]=4*LY*yy+MM*MM;B[1][1]=2*LY*y*y+2*conj(u)*s;
 A[1][2]=2*conj(u)*conj(y);A[2][1]=conj(A[1][2]);
 A[2][2]=4*LS*ss-LS*HH*HH;B[2][2]=2*LS*s*s;
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){D[2*i][2*j]=creal(A[i][j]+B[i][j]);D[2*i][2*j+1]=-cimag(A[i][j]-B[i][j]);D[2*i+1][2*j]=cimag(A[i][j]+B[i][j]);D[2*i+1][2*j+1]=creal(A[i][j]-B[i][j]);}
}
static void rhs(const double *v,double *out){
 double infl=v[0],vel=v[1],rad=v[2],a=exp(v[3]);double H=sqrt((.5*vel*vel+4.5*infl*infl+rad)/3.),R=18*infl*infl-vel*vel,Rdot=54*infl*vel+2*(3*H+.5)*vel*vel;
 double V,F[6],VR,D[6][6];matter(v+4,R,&V,F,&VR,D);
 out[0]=vel;out[1]=-(3*H+.5)*vel-9*infl;out[2]=.5*vel*vel-4*H*rad;out[3]=H;out[16]=Rdot*VR;
 for(int i=0;i<6;i++){out[4+i]=v[10+i];out[10+i]=-3*H*v[10+i]-F[i];out[16]-=3*H*v[10+i]*v[10+i];}
 if(!nk)return;
 double M[6][6],chi[6],chip[6],g[6],gp[6],ward[6],TF[6],gg=0;
 charge_rotate(F,TF);
 for(int i=0;i<6;i++){chi[i]=a*v[4+i];chip[i]=a*a*(v[10+i]+H*v[4+i]);for(int j=0;j<6;j++)M[i][j]=a*a*(D[i][j]-(i==j?R/6:0));}
 charge_rotate(chi,g);charge_rotate(chip,gp);
 for(int i=0;i<6;i++){ward[i]=a*a*a*TF[i]-a*a*R*g[i]/6;gg+=g[i]*g[i];}
 for(int n=0;n<nk;n++){
  const double *U=v+17+144*n;double *du=out+17+144*n;double k2=wave[n]*wave[n],c=k2?eg*eg/k2:0;
  for(int col=0;col<12;col++){
   double current=0,common=0;for(int i=0;i<6;i++){current+=g[i]*U[(i+6)*12+col]-gp[i]*U[i*12+col];common+=g[i]*U[i*12+col];}common/=gg;
   for(int i=0;i<6;i++){
    double mq=common*ward[i];for(int j=0;j<6;j++)mq+=M[i][j]*(U[j*12+col]-common*g[j]);
    du[i*12+col]=(U[(i+6)*12+col]+c*g[i]*current)/a;
    du[(i+6)*12+col]=(-k2*U[i*12+col]-mq+c*gp[i]*current)/a;
   }
  }
 }
}
int main(int argc,char **argv){
 if(argc!=3)return 1;
 FILE *f=fopen(argv[1],"rb");if(!f)return 2;
 int32_t head[2];if(fread(head,4,2,f)!=2)return 3;nk=head[0];int nt=head[1];if(nk<0||nk>16||nt<2)return 4;
 double pars[6];if(fread(pars,8,6,f)!=6)return 5;double end=pars[0],rtol=pars[1],atol=rtol*.001;eg=pars[2];cr=pars[3];ci=pars[4];double maxstep=pars[5];
 nv=17+144*nk;double *v=calloc(nv,8),*next=calloc(nv,8),*tmp=calloc(nv,8),*st=calloc(13*nv,8),*times=calloc(nt,8);wave=calloc(nk,8);
 if(fread(v,8,17,f)!=17||fread(wave,8,nk,f)!=(size_t)nk||fread(times,8,nt,f)!=(size_t)nt)return 6;
 fclose(f);
 for(int n=0;n<nk;n++)for(int i=0;i<12;i++)v[17+144*n+12*i+i]=1;
 FILE *out=fopen(argv[2],"wb");if(!out)return 7;
 double t=0,step=.0001;long accepted=0,rejected=0;rhs(v,st);
 for(int sample=0;sample<nt;sample++){
  while(t<times[sample]){
   double dt=fmin(step,times[sample]-t);if(dt<1e-15)return 8;
   for(int stage=1;stage<12;stage++){
    for(int j=0;j<nv;j++){double sum=0;for(int k=0;k<stage;k++)sum+=AA[stage][k]*st[k*nv+j];tmp[j]=v[j]+dt*sum;}rhs(tmp,st+stage*nv);
   }
   for(int j=0;j<nv;j++){double sum=0;for(int k=0;k<12;k++)sum+=BB[k]*st[k*nv+j];next[j]=v[j]+dt*sum;}
   rhs(next,st+12*nv);double e3sq=0,e5sq=0;
   for(int j=0;j<nv;j++){double e3=0,e5=0;for(int k=0;k<13;k++){e3+=EE3[k]*st[k*nv+j];e5+=EE5[k]*st[k*nv+j];}double sc=atol+rtol*fmax(fabs(v[j]),fabs(next[j]));e3sq+=pow(e3/sc,2);e5sq+=pow(e5/sc,2);}
   double err=(e3sq==0&&e5sq==0)?0:fabs(dt)*e5sq/sqrt((e5sq+.01*e3sq)*nv);
   if(!isfinite(err)){fprintf(stderr,"nonfinite at %.17g\n",t);return 9;}
   double fac=err?fmax(.2,fmin(5.,.9*pow(err,-.125))):5;
   if(err<=1){t+=dt;memcpy(v,next,nv*8);memcpy(st,st+12*nv,nv*8);accepted++;}else rejected++;
   step=fmin(maxstep,dt*fac);if(accepted>100000000)return 10;
  }
  if(fwrite(&t,8,1,out)!=1||fwrite(v,8,nv,out)!=(size_t)nv)return 11;
 }
 fclose(out);printf("COMPLETE end=%.17g accepted=%ld rejected=%ld\n",end,accepted,rejected);free(v);free(next);free(tmp);free(st);free(times);free(wave);return 0;
}
