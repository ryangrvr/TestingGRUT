"""en_rerun INDEPENDENT check: linear tangent growth on the reduced homogeneous X orbit, tau 20 -> 60,
for k inside and below the seed band (1.8 <= k <= 8). Own code; reads only the matched mean field and clock
from Codex's MATCH_L12.npz. Purpose: is the omitted IR band (k < 1.8) as unstable as the seeded band edge?
usage: python3 -I ir_tangent_scan.py <packet_dir> <out_json>"""
import sys,json
import numpy as np
from scipy.integrate import solve_ivp
pkt,outp=sys.argv[1],sys.argv[2];LX=.001
m=np.load(pkt+'/MATCH_L12.npz');a_m=float(m['a_match']);C=float(m['C_match'])
q0=m['light_mean_q'][:2];p0=m['light_mean_p'][:2]
def s_of_tau(t):return 2*(t-20)/(a_m+np.sqrt(a_m*a_m+2*C*(t-20)))
ks=np.r_[np.linspace(.2,1.8,17),np.linspace(1.9,4,22),5.,6.,8.]
# lattice momenta below 1.8 at L=12 and L=16
lat12=sorted({2*np.pi/12*np.sqrt(n) for n in range(1,40) if 2*np.pi/12*np.sqrt(n)<1.8 and n%8!=7 and n not in (28,)})
ks=np.unique(np.r_[ks,lat12])
nk=len(ks)
def rhs(s,y):
    x=y[:2];v=y[2:4];r2=x@x;T=y[4:].reshape(nk,4,4);out=np.empty_like(y);out[:2]=v;out[2:4]=-LX*r2*x
    M=LX*r2*np.eye(2)+2*LX*np.outer(x,x)
    dT=np.empty_like(T);dT[:,:2,:]=T[:,2:,:]
    dT[:,2:,:]=-(ks[:,None,None]**2*T[:,:2,:]+np.einsum('ij,kjl->kil',M,T[:,:2,:]))
    out[4:]=dT.ravel();return out
taus=np.array([30.,40.,50.,55.,60.])
y0=np.r_[q0,p0,np.tile(np.eye(4).ravel(),nk)]
sol=solve_ivp(rhs,(0,s_of_tau(60.)),y0,t_eval=s_of_tau(taus),method='DOP853',rtol=1e-10,atol=1e-12)
assert sol.success
growth={}
for j,t in enumerate(taus):
    T=sol.y[4:,j].reshape(nk,4,4);growth[str(t)]=[float(np.linalg.norm(Ti,2)) for Ti in T]
g55=np.array(growth['55.0'])
res=dict(ks=ks.tolist(),lattice_L12_k_below_1p8=lat12,max_singular_value_of_tangent_map=growth,
    k_of_max_growth_tau55=float(ks[np.argmax(g55)]),max_growth_tau55=float(g55.max()),
    growth_tau55_at_k_1p737=float(g55[np.argmin(abs(ks-1.737))]),growth_tau55_at_k_1p8=float(g55[np.argmin(abs(ks-1.8))]),
    growth_tau55_at_k_2p0=float(g55[np.argmin(abs(ks-2.0))]),
    note='Reduced quartic homogeneous orbit (no backreaction); maximum singular value of the 4x4 X tangent map from tau20. '
         'A comparable growth below k=1.8 means the omitted IR seed is a missing sensitivity control (onset clock), not an error.')
json.dump(res,open(outp,'w'),indent=1)
for k,g in zip(ks,g55):print(f'k={k:.3f} growth55={g:.3e} growth60={growth["60.0"][list(ks).index(k)]:.3e}')
