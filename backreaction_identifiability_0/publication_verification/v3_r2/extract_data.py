"""Extract complete authoritative data at reference resolution nx=240, dt=5e-4."""
import numpy as np, json
T_STAR=0.5; T_DIAG=[0.25,0.75,1.0]; NB_GRID=[4,8,16,32,64,128]
L=6.0; NX=240; DT=5e-4
def s(u): return 10*u**3-15*u**4+6*u**5
def q_p1(t): return s(t/np.pi) if t<=np.pi else 1.0
xg,wg = np.polynomial.legendre.leggauss(NX)
xs = L*xg; wx = L*wg
X0 = np.broadcast_to(xs[:,None],(NX,NX)).copy()
P0 = np.broadcast_to(xs[None,:],(NX,NX)).copy()
rho = np.exp(-0.5*P0**2-0.5*X0**2-0.25*X0**4)
W = wx[:,None]*wx[None,:]
Wn_raw = W*rho; norm_raw = Wn_raw.sum()
Wn = Wn_raw/norm_raw
def integrate(X0,P0,T,eps,dt):
    steps=max(1,int(round(T/dt))); h=T/steps
    x=X0.copy(); v=P0.copy()
    for k in range(steps):
        q=q_p1(k*h)
        k1v=-x-x**3+eps*q; k1x=v
        x2=x+0.5*h*k1x; v2=v+0.5*h*k1v
        qm=q_p1(k*h+0.5*h)
        k2v=-x2-x2**3+eps*qm; k2x=v2
        x3=x+0.5*h*k2x; v3=v+0.5*h*k2v
        k3v=-x3-x3**3+eps*qm; k3x=v3
        x4=x+h*k3x; v4=v+h*k3v
        q1=q_p1(k*h+h)
        k4v=-x4-x4**3+eps*q1; k4x=v4
        x=x+(h/6)*(k1x+2*k2x+2*k3x+k4x)
        v=v+(h/6)*(k1v+2*k2v+2*k3v+k4v)
    return x
def moments(f,Wn):
    m1=(f*Wn).sum(); m2=(f**2*Wn).sum(); m3=(f**3*Wn).sum(); m4=(f**4*Wn).sum()
    var=m2-m1*m1; k3=m3-3*m1*m2+2*m1*m1*m1
    return m1,m2,var,k3,m4
out={"settings":{"nx":NX,"dt":DT,"norm_raw":float(norm_raw)}}
# P0 at all times
p0={}
for T in [0.0,0.25,0.5,0.75,1.0]:
    xT=integrate(X0,P0,T,0.0,DT)
    m1,m2,var,k3,m4=moments(xT,Wn)
    p0[str(T)]={"mean":float(m1),"m2":float(m2),"var":float(var),"k3":float(k3),"m4":float(m4)}
out["p0"]=p0
out["gibbs_identity_residual"]=float(p0["0.0"]["var"]+p0["0.0"]["m4"]-1.0)
# Independent 1-D check
xm,wxm=np.polynomial.legendre.leggauss(4000)
xm=L*xm; wxm=L*wxm
rx=np.exp(-0.5*xm**2-0.25*xm**4); wxn=wxm*rx; wxn/=wxn.sum()
out["independent_1d"]={"m2":float((xm**2*wxn).sum()),"m4":float((xm**4*wxn).sum()),"m6":float((xm**6*wxn).sum())}
# P1 all times
p1={}
for T in [T_STAR]+T_DIAG:
    row={}
    for nb in NB_GRID:
        eps=1.0/np.sqrt(nb)
        xT=integrate(X0,P0,T,eps,DT)
        m1,m2,var,k3,m4=moments(xT,Wn)
        k3F=k3/np.sqrt(nb); g1F=k3F/var**1.5
        row[str(nb)]={"epsilon":float(eps),"mean_X":float(m1),"var_X":float(var),
                      "kappa3_X":float(k3),"kappa3_F":float(k3F),
                      "gamma1_F":float(g1F),"NB_gamma1":float(nb*g1F)}
        # assertions
        assert abs(k3F-k3/np.sqrt(nb))<1e-30
        assert abs(g1F-k3F/var**1.5)<1e-15
    p1[str(T)]=row
out["p1"]=p1
# Odd-epsilon
odd={}
for nb in [4,16]:
    eps=1.0/np.sqrt(nb)
    xTp=integrate(X0,P0,T_STAR,eps,DT)
    xTm=integrate(X0,P0,T_STAR,-eps,DT)
    _,_,_,k3p,_=moments(xTp,Wn)
    _,_,_,k3m,_=moments(xTm,Wn)
    odd[str(nb)]={"k3_plus":float(k3p),"k3_minus":float(k3m),
                  "sum":float(k3p+k3m),"ratio":float(k3m/k3p)}
out["odd_epsilon"]=odd
# Scaling fit
gam=np.array([p1[str(T_STAR)][str(nb)]["gamma1_F"] for nb in NB_GRID])
nb_arr=np.array(NB_GRID,dtype=float)
p_g,a_g=np.polyfit(np.log(nb_arr),np.log(np.abs(gam)),1)
local_p=[]
for i in range(len(nb_arr)-1):
    lp=-np.log(abs(gam[i+1])/abs(gam[i]))/np.log(nb_arr[i+1]/nb_arr[i])
    local_p.append(float(lp))
out["scaling_fit"]={"global_p":float(-p_g),"intercept":float(a_g),
                    "local_p":local_p,
                    "NB_gamma1":[float(nb*g) for nb,g in zip(nb_arr,gam)]}
json.dump(out,open("v3_r2_results.json","w"),indent=2)
print("SAVED v3_r2_results.json")
print(f"Gibbs identity residual: {out['gibbs_identity_residual']:.3e}")
print(f"Independent 1D: m2={out['independent_1d']['m2']:.12f} m4={out['independent_1d']['m4']:.12f} m6={out['independent_1d']['m6']:.12f}")
print(f"Odd-eps nb=4: k3+={odd['4']['k3_plus']:.6e} k3-={odd['4']['k3_minus']:.6e} sum={odd['4']['sum']:.3e}")
print(f"Odd-eps nb=16: k3+={odd['16']['k3_plus']:.6e} k3-={odd['16']['k3_minus']:.6e} sum={odd['16']['sum']:.3e}")
print(f"Global p={float(-p_g):.6f}, local_p={local_p}")
print(f"NB_gamma1: {[f'{v:.6e}' for v in out['scaling_fit']['NB_gamma1']]}")
print("P1 t*=0.5:")
for nb in NB_GRID:
    r=p1[str(T_STAR)][str(nb)]
    print(f"  {nb}: eps={r['epsilon']:.4f} m1={r['mean_X']:.6e} var={r['var_X']:.8f} k3X={r['kappa3_X']:.6e} k3F={r['kappa3_F']:.6e} g1F={r['gamma1_F']:.6e} NB*g={r['NB_gamma1']:.6e}")
