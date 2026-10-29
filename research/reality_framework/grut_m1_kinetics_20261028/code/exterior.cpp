#include <complex>
#include <cmath>
extern "C" void ratios(double *k2, int ell, double h, double R, double Rfar, double *out) {
    int N=std::lround(R/h), far=std::max(N+1,(int)std::lround(Rfar/h));
    for(int c=0;c<3;++c){
        std::complex<double> q;
        if(k2[c]>0) q=std::exp(std::complex<double>(0,2*std::asin(.5*h*std::sqrt(k2[c]))));
        else q=std::exp(-2*std::asinh(.5*h*std::sqrt(-k2[c])));
        for(int i=far;i>N;--i) q=1.0/(2+h*h*(ell*(ell+1)/(i*h*i*h)-k2[c])-q);
        out[2*c]=q.real();out[2*c+1]=q.imag();
    }
}
