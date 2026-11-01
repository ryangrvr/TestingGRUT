import sys, math, numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp
TRACK = Path('/home/user/TestingGRUT/research/reality_framework/grut_formation_review_20261025/tracks/m1_transport')
sys.path.insert(0, str(TRACK/'code'))
import mixture as mx, transport_check as tc
bg=mx.bg; a_i=mx.a_i; f=mx.fgas; cs2=mx.v_rms**2/3
def trace(k, x_dec):
    N_end=math.log(1/a_i); N_dec=math.log(x_dec)
    aHi=mx.aH_Mpc(a_i); a_k=a_i*aHi/k if k>aHi else a_i; slope=1/(max(math.log(a_i/a_k),0)+1)
    def coef(N):
        a=a_i*math.exp(N); aH=mx.aH_Mpc(a); return 2+mx.dlnH_dlna(a), tc.Om/a**3/float(bg.E2(a)), (k*math.sqrt(cs2)*a_i/a/aH)**2
    def rhs(N,y,mode):
        fric,Om,kc=coef(N); dc,dcp,dg,dgp=y
        if mode=='cold': return [dcp,-fric*dcp+1.5*Om*dc,0,0]
        if mode=='tight': s=1.5*Om*dc; return [dcp,-fric*dcp+s-f*kc*dc,dcp,-fric*dcp+s-f*kc*dc]
        s=1.5*Om*((1-f)*dc+f*dg); return [dcp,-fric*dcp+s,dgp,-fric*dgp+s-kc*dg]
    s1=solve_ivp(rhs,(0,N_dec),[1,slope,1,slope],args=('tight',),method='Radau',rtol=1e-10,atol=1e-14,dense_output=True)
    y1=s1.y[:,-1]; y1=[y1[0],y1[1],y1[0],y1[1]]
    s2=solve_ivp(rhs,(N_dec,N_end),y1,args=('two',),method='Radau',rtol=1e-10,atol=1e-14,dense_output=True)
    sc=solve_ivp(rhs,(0,N_end),[1,slope,0,0],args=('cold',),method='Radau',rtol=1e-10,atol=1e-14,dense_output=True)
    N_eq=math.log(tc.a_eq/a_i)
    print(f'k={k:g} x_dec={x_dec:g} slope={slope:.3f} kc_s/aH(2MeV)={math.sqrt(coef(0)[2]):.3f}')
    for N in [N_dec, N_dec+2, N_dec+4, N_eq-1, N_eq, N_eq+1, N_eq+2, N_eq+4, N_eq+6, N_end]:
        if N<N_dec: continue
        y=s2.sol(min(N,N_end)); c=sc.sol(min(N,N_end))[0]; a=a_i*math.exp(N)
        print(f'  a={a:.2e} dc={y[0]:9.3f} dcp={y[1]:8.3f} dg={y[2]:9.3f} dgp={y[3]:8.3f} cold={c:9.3f} dm/cold={((1-f)*y[0]+f*y[2])/c:.4f} kcs/aH={math.sqrt(coef(N)[2]):.3f}')
trace(3e5, 2073.)
trace(3e5, 1.0001)
trace(1e6, 6.55)
