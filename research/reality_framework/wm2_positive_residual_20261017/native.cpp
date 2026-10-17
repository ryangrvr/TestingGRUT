#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cfenv>
#include <limits>

struct Z { double r,i; };
static_assert(sizeof(Z)==2*sizeof(double),"complex buffer layout");
inline Z add(Z a,Z b){return {a.r+b.r,a.i+b.i};}
inline Z mul(Z a,Z b){return {a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r};}
inline Z conjz(Z a){return {a.r,-a.i};}
inline Z scale(Z a,double b){return {a.r*b,a.i*b};}

extern "C" {
int round_to_nearest(){return std::fegetround()==FE_TONEAREST;}
int long_double_digits(){return std::numeric_limits<long double>::digits;}

double positive_norm(int64_t n,const Z*x){
    long double value=0;
    for(int64_t j=0;j<n;j++){
        long double a=x[j].r,b=x[j].i; value+=a*a+b*b;
    }
    return static_cast<double>(value);
}

void action(int64_t total,int words,const int32_t*lengths,const int64_t*offset,
            const int32_t*parent,const double*lambda,const int64_t*nocc,
            const int64_t*annoff,const int32_t*modes,const int32_t*lower,
            const double*factor,const Z*contact,const Z*x,Z*y){
    std::memset(y,0,sizeof(Z)*total);
    for(int v=1;v<words;v++){
        int k=lengths[v]; int64_t base=annoff[k],count=nocc[k];
        int64_t dst0=offset[v],src0=offset[parent[v]];
        for(int64_t j=0;j<count;j++){
            int64_t dst=dst0+j; Z value=x[dst];
            for(int p=0;p<k;p++){
                int64_t a=base+j*k+p;int m=modes[a];if(m<0)continue;
                int64_t src=src0+lower[a];
                Z coefficient=scale(contact[m],-lambda[v]*factor[a]);
                y[dst]=add(y[dst],mul(coefficient,x[src]));
                y[src]=add(y[src],mul(conjz(coefficient),value));
            }
        }
    }
}

void fill_density(int words,int rank,const int32_t*lengths,const int64_t*offset,
                  const int64_t*auxoff,const int64_t*nocc,const int64_t*annoff,
                  const int32_t*modes,const int32_t*lower,const double*factor,
                  const Z*x,Z*d){
    for(int v=1;v<words;v++){
        int k=lengths[v];int64_t start=annoff[k];
        for(int64_t j=0;j<nocc[k];j++){
            Z value=x[offset[v]+j];
            for(int p=0;p<k;p++){
                int64_t a=start+j*k+p;int m=modes[a];if(m<0)continue;
                d[(auxoff[v]+lower[a])*rank+m]=scale(value,factor[a]);
            }
        }
    }
}

void vacancy_weights(int words,int rank,const int32_t*lengths,const int64_t*offset,
                     const int64_t*nocc,const int64_t*annoff,const int32_t*modes,
                     const Z*x,double*q){
    long double values[64]={0};
    for(int v=1;v<words;v++){
        int k=lengths[v];int64_t start=annoff[k];
        for(int64_t j=0;j<nocc[k];j++){
            Z z=x[offset[v]+j];long double r=z.r,i=z.i,value=r*r+i*i;
            for(int p=0;p<k;p++){
                int m=modes[start+j*k+p];if(m>=0)values[m]+=value;
            }
        }
    }
    for(int m=0;m<rank;m++)q[m]=static_cast<double>(values[m]);
}

void beam(int64_t rows,int64_t columns,int block,int64_t groups,
          const int32_t*index,const Z*matrix,Z*x){
    for(int64_t row=0;row<rows;row++){
        Z*data=x+row*columns;
        for(int64_t g=0;g<groups;g++){
            Z old[8];
            for(int j=0;j<block;j++)old[j]=data[index[g*block+j]];
            for(int j=0;j<block;j++){
                Z value{0,0};
                for(int l=0;l<block;l++)value=add(value,mul(matrix[j*block+l],old[l]));
                data[index[g*block+j]]=value;
            }
        }
    }
}
}
