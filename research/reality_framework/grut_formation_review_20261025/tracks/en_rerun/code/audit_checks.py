"""en_rerun audit checks (INDEPENDENT code by Claude; reads Codex outputs/initial states, never imports Codex modules).

usage: python3 -I audit_checks.py <packet_dir> <early_workdir> <out_json>

A. radiation clock tau(s)=20+a_m s+C s^2/2 versus the full background a(tau) in dense_b3.npz (20<=tau<=240)
B. own real-space FD4 energy (np.roll stencil) versus the C++ trace rho at tau 20 and 800
C. whitened empirical covariance of the generated N48 initial Fourier coefficients versus the specified
   light marginal C_light/V from MATCH_L12.npz (per-mode and pooled); pseudo-covariance E[c c^T]
D. symplectic eigenvalues of the full 12x12 tau-20 covariance (pure state -> all 1/2) and of the light X marginal
E. reference projections (massless / isotropic Hartree) at tau 20: seeded zero-point content per net charge
F. own n_pm implementation (full complex FFT, no half-space weights) at tau 800, plus mass2 variants:
   Hartree with/without the zero mode, radial 3 lam <r2>, tangential lam <r2>, scan
"""
import sys,json,struct
import numpy as np

pkt,work,outp=sys.argv[1],sys.argv[2],sys.argv[3]
LX,G2,LP=.001,.01,.03
res={}

def fd4_symbol(N,L):
    """own FD4 symbol from the stencil (-1,16,-30,16,-1)/12 per axis, full complex FFT grid"""
    k=2*np.pi*np.fft.fftfreq(N)
    s=(30-32*np.cos(k)+2*np.cos(2*k))/12/(L/N)**2
    return s[:,None,None]+s[None,:,None]+s[None,None,:]

def lap_real(f,dx):
    out=-7.5*f
    for ax in range(3):
        out=out+(4/3)*(np.roll(f,1,ax)+np.roll(f,-1,ax))-(1/12)*(np.roll(f,2,ax)+np.roll(f,-2,ax))
    return out/dx**2

def energy_real(q,p,L):
    dx=L/q.shape[-1];r2=q[0]**2+q[1]**2
    grad=sum(-.5*np.mean(q[c]*lap_real(q[c],dx)) for c in range(3))
    return float(.5*np.mean(np.sum(p*p,axis=0))+grad+np.mean(LX/4*r2*r2+G2/2*q[2]**2*r2+LP/4*q[2]**4))

def read_initial(path):
    with open(path,'rb') as f:
        magic=f.read(8);N,=struct.unpack('<i',f.read(4));L,=struct.unpack('<d',f.read(8));par=struct.unpack('<4d',f.read(32))
        a=np.fromfile(f,dtype='<f8').reshape(6,N,N,N)
    return N,L,par,a[:3].copy(),a[3:].copy()

def read_snapshot(path):
    with open(path,'rb') as f:
        N,=struct.unpack('<i',f.read(4));L,=struct.unpack('<d',f.read(8));s,=struct.unpack('<d',f.read(8));bits,=struct.unpack('<i',f.read(4))
        a=np.fromfile(f,dtype='<f8').reshape(6,N,N,N)
    return N,L,s,a[:3].copy(),a[3:].copy()

def trace(tag):
    return np.genfromtxt(f'{work}/{tag}/state.tsv',names=True,delimiter='\t')

m=np.load(pkt+'/MATCH_L12.npz');a_m=float(m['a_match']);C=float(m['C_match'])
# ---------- A. clock
d=np.load(pkt+'/dense_b3.npz');t=d['times'];a=d['a'];sel=(t>=20)&(t<=240)
aclock=np.sqrt(a_m**2+2*C*(t[sel]-20))
res['A_clock']=dict(a_match=a_m,C_match=C,max_rel_a_deviation_20_240=float(np.max(abs(aclock/a[sel]-1))),
    a_full_240=float(a[sel][-1]),a_clock_240=float(aclock[-1]),
    note='a(s)=a_m+C s with d tau=a d eta gives a(tau)=sqrt(a_m^2+2C(tau-20)); compared with the full quiet(3,.5) background')

# ---------- B. energies
N,L,par,q,p=read_initial(f'{work}/N48/initial.bin');tr=trace('N48')
e0=energy_real(q,p,L)
snaps=sorted(__import__('glob').glob(f'{work}/N48/state_s*.bin'));Nf,Lf,sf,qf,pf=read_snapshot(snaps[-1])
ef=energy_real(qf,pf,Lf)
res['B_energy']=dict(own_FD4_energy_tau20=e0,cpp_trace_rho_tau20=float(tr['rho'][0]),rel20=abs(e0/tr['rho'][0]-1),
    own_FD4_energy_tau800=ef,cpp_trace_rho_tau800=float(tr['rho'][-1]),rel800=abs(ef/tr['rho'][-1]-1),
    snapshot_eta=sf,trace_eta_last=float(tr['eta'][-1]),file_params_nu_g2_lx_lp=list(par))

# ---------- C. initial covariance test (N48, L=12, cutoff 8)
V=L**3;keys=m['keys'];lookup={int(v):i for i,v in enumerate(keys)}
ES=m['evolved_scalar'];EP=m['evolved_phi']
def light_sqrt(ix):
    """rows of the 6 x 14 map z -> (qX1,qX2,qphi,pX1,pX2,pphi)"""
    A=np.zeros((6,14));A[[0,1,3,4],:12]=ES[ix][[0,1,6,7],:];A[[2,5],12:]=EP[ix][[0,1],:]
    return A
ch=np.fft.fftn(q,axes=(1,2,3),norm='forward');cp=np.fft.fftn(p,axes=(1,2,3),norm='forward')
c=np.concatenate([ch,cp],axis=0)  # order qX1,qX2,qphi,pX1,pX2,pphi
n=np.rint(np.fft.fftfreq(N)*N).astype(int)
I,J,K=np.meshgrid(n,n,n,indexing='ij');nn=(I**2+J**2+K**2).astype(int);kk=2*np.pi/L*np.sqrt(nn)
# one representative of each +-n pair: first nonzero component positive
first=np.where(I!=0,I,np.where(J!=0,J,K))
band=(kk>=1.8)&(kk<=8.)&(first>0)
idx=np.argwhere(band)
W=np.zeros((6,6),complex);Wpseudo=np.zeros((6,6),complex);count=0;chi2=[];per_shell={}
outmask=(~((kk>=1.8)&(kk<=8.)))&(nn>0);outside=float(np.max(abs(c[:,outmask])))  # excluding k=0 (mean)
for (i,j,k) in idx:
    ix=lookup[int(nn[i,j,k])];A=light_sqrt(ix)
    Cl=A@A.T/V  # specified light covariance of one complex coefficient: E[c c^dagger]
    Lc=np.linalg.cholesky(Cl);w=np.linalg.solve(Lc,c[:,i,j,k])
    W+=np.outer(w,w.conj());Wpseudo+=np.outer(w,w);count+=1
    per_shell.setdefault(ix,[]).append(np.real(w.conj()@w))
W/=count;Wpseudo/=count
res['C_initial_covariance']=dict(independent_pairs=count,
    whitened_cov_max_abs_dev_from_identity=float(np.max(abs(W-np.eye(6)))),
    whitened_cov_diag=np.real(np.diag(W)).tolist(),
    whitened_pseudo_cov_max_abs=float(np.max(abs(Wpseudo))),
    expected_statistical_sigma=float(1/np.sqrt(count)),
    mean_whitened_norm2_per_component=float(np.mean([np.mean(v) for v in per_shell.values()])/6),
    max_abs_coefficient_outside_band_excluding_k0=outside,
    note='E[c c^dag]=C_light/V and E[c c^T]=0 are the intended complex-Gaussian moments for one +-n pair; '
         'whitened coefficients should give identity and zero pseudo-covariance within ~1/sqrt(pairs)')
# variance factor check in real space: sum over all modes of specified C_qq/V vs empirical <dq^2>
spec_qX=0.;spec_pX=0.
for (i,j,k) in idx:
    A=light_sqrt(lookup[int(nn[i,j,k])]);Cl=A@A.T/V;spec_qX+=2*(Cl[0,0]+Cl[1,1]);spec_pX+=2*(Cl[3,3]+Cl[4,4])
dq=q[:2]-q[:2].mean(axis=(1,2,3),keepdims=True);dp=p[:2]-p[:2].mean(axis=(1,2,3),keepdims=True)
res['C_initial_covariance'].update(specified_X_fluct_variance=spec_qX,empirical_X_fluct_variance=float(np.mean(np.sum(dq*dq,axis=0))),
    specified_piX_fluct_variance=spec_pX,empirical_piX_fluct_variance=float(np.mean(np.sum(dp*dp,axis=0))))

# ---------- D. symplectic eigenvalues
def symp_eigs(Cm):
    n2=Cm.shape[0]//2;Jm=np.block([[np.zeros((n2,n2)),np.eye(n2)],[-np.eye(n2),np.zeros((n2,n2))]])
    ev=np.linalg.eigvals(1j*Jm@Cm);return np.sort(abs(ev.real))[::2]
full=[];light=[];ks=m['ks']
for ix in range(len(ks)):
    Cf=ES[ix]@ES[ix].T;full.append(symp_eigs(Cf));CX=Cf[np.ix_([0,1,6,7],[0,1,6,7])];light.append(symp_eigs(CX))
full=np.array(full);light=np.array(light);inband=(ks>=1.8)&(ks<=8)
res['D_symplectic']=dict(full_12x12_symplectic_eigs_max_abs_dev_from_half=float(np.max(abs(full[inband]-.5))),
    light_X_marginal_symplectic_eigs_min=float(light[inband].min()),light_X_marginal_symplectic_eigs_max=float(light[inband].max()),
    light_X_marginal_mean_occupation_excess=float(np.mean(light[inband]-.5)),
    k_of_max_mixedness=float(ks[inband][np.argmax(light[inband].max(axis=1))]),
    note='Full scalar covariance at tau20 remains pure (1/2). X-marginal symplectic eigenvalue nu>1/2 measures heavy-field entanglement retained in the light marginal (nu-1/2 = thermal-like excess per mode)')

# ---------- E/F. projections
def projections(q,p,L,mass2_list):
    N=q.shape[-1];k2=fd4_symbol(N,L);qh=np.fft.fftn(q,axes=(1,2,3),norm='forward');ph=np.fft.fftn(p,axes=(1,2,3),norm='forward')
    nz=k2>0;Qt=float(np.mean(q[0]*p[1]-q[1]*p[0]))
    Qm=np.real(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0])
    out={}
    for name,m2 in mass2_list:
        w=np.sqrt(np.where(nz,k2+m2,1.))
        nr=(abs(ph[0])**2+abs(ph[1])**2+w**2*(abs(qh[0])**2+abs(qh[1])**2))/(2*w)
        nm=float(np.sum(((nr-Qm)/2)[nz]));npl=float(np.sum(((nr+Qm)/2)[nz]))
        out[name]=dict(mass2=m2,nminus_over_net=nm/Qt,nplus_over_net=npl/Qt,charge_identity=(npl-nm)/Qt)
    return out,Qt
# tau20 (initial) for N48
r2=q[0]**2+q[1]**2;hf0=2*LX*np.mean(r2)+G2*np.mean(q[2]**2)
p20,Q20=projections(q,p,L,[('massless',0.),('isotropic_Hartree',hf0)])
modes=2*count;Qtot=Q20*V
res['E_zero_point']=dict(tau20_projections=p20,net_current_density=Q20,net_charge_total=Qtot,seeded_modes=modes,
    half_quantum_per_seeded_mode_nminus_per_net=modes*.5/Qtot,
    full_lattice_half_quantum_nminus_per_net=(N**3-1)*.5/Qtot,
    note='Classical seeding gives n=1/2 per seeded mode per charge sign (no subtraction).')
# tau800 variants
r2f=qf[0]**2+qf[1]**2;mean2=float(qf[0].mean()**2+qf[1].mean()**2);X2=float(np.mean(r2f));phi2=float(np.mean(qf[2]**2))
variants=[('massless',0.),('isotropic_Hartree_with_zero_mode(Codex)',2*LX*X2+G2*phi2),
          ('Hartree_fluctuations_only',2*LX*(X2-mean2)+G2*phi2),
          ('tangential_lambda_r2',LX*X2+G2*phi2),('radial_3lambda_r2',3*LX*X2+G2*phi2)]
pf800,Q800=projections(qf,pf,Lf,variants)
scan=projections(qf,pf,Lf,[(f'm2={x:g}',x) for x in [0,1,2,5,10,19.26,30,50,100]])[0]
res['F_tau800']=dict(own_projection=pf800,mass2_scan={k:v['nminus_over_net'] for k,v in scan.items()},
    X2_total=X2,X2_zero_mode=mean2,zero_mode_fraction_of_X2=mean2/X2,phi2=phi2,net_current_density=Q800,
    zero_point_seeded_share_of_massless=(modes*.5/(Q800*V))/pf800['massless']['nminus_over_net'])
# n_minus spectrum in physical-momentum bins
k2=fd4_symbol(Nf,Lf);qh=np.fft.fftn(qf,axes=(1,2,3),norm='forward');ph=np.fft.fftn(pf,axes=(1,2,3),norm='forward')
Qm=np.real(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0]);nzm=k2>0
kphys=2*np.pi/Lf*np.sqrt(nn.astype(float)) if Nf==N else None
edges=np.array([0,1,1.8,3,4.4,6,8,10,1e9]);spec={}
for name,m2 in [('massless',0.),('Hartree',variants[1][1])]:
    w=np.sqrt(np.where(nzm,k2+m2,1.));nr=(abs(ph[0])**2+abs(ph[1])**2+w**2*(abs(qh[0])**2+abs(qh[1])**2))/(2*w)
    nm=np.where(nzm,(nr-Qm)/2,0.)
    spec[name]=[float(nm[(kphys>=lo)&(kphys<hi)].sum()/Q800) for lo,hi in zip(edges[:-1],edges[1:])]
res['F_tau800']['nminus_per_net_by_k_bin']=dict(edges=edges.tolist(),**spec)
json.dump(res,open(outp,'w'),indent=1,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else float(x))
print(json.dumps(res,indent=1,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else float(x)))
