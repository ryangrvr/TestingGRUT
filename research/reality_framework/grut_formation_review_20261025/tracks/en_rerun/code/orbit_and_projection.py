"""en_rerun INDEPENDENT analysis of the REPRODUCED traces/snapshots (own code; Codex modules not imported).
For each history: orbit (pericentre) phase of the coherent mean, cycle-averaged fluctuation diagnostics,
late-window averages, and reference n_pm projections with mass^2 variants and IR/UV split.
usage: python3 -I orbit_and_projection.py <early_workdir> <out_json> tag [tag...]"""
import sys,json,struct,glob,os
import numpy as np
from scipy.signal import argrelextrema
work,outp=sys.argv[1],sys.argv[2];tags=sys.argv[3:]
LX,G2,LP=.001,.01,.03

def fd4_symbol(N,L):
    k=2*np.pi*np.fft.fftfreq(N);s=(30-32*np.cos(k)+2*np.cos(2*k))/12/(L/N)**2
    return s[:,None,None]+s[None,:,None]+s[None,None,:]

def read_snapshot(path):
    with open(path,'rb') as f:
        N,=struct.unpack('<i',f.read(4));L,=struct.unpack('<d',f.read(8));s,=struct.unpack('<d',f.read(8));bits,=struct.unpack('<i',f.read(4))
        a=np.fromfile(f,dtype='<f8').reshape(6,N,N,N)
    return N,L,s,a[:3].copy(),a[3:].copy()

out={}
for tag in tags:
    d=os.path.realpath(f'{work}/{tag}');t=np.genfromtxt(d+'/state.tsv',names=True,delimiter='\t')
    tau=t['tau'];rho=t['rho'];mr2=t['q0']**2+t['q1']**2;X2=t['X2'];dX=X2-mr2
    EX0=.5*(t['p0']**2+t['p1']**2)+LX/4*mr2**2+G2/2*t['qphi']**2*mr2;EP0=.5*t['pphi']**2+LP/4*t['qphi']**4
    outside=1-(EX0+EP0)/rho;pm=(t['Kphi']+t['Gphi']+t['Uphi']+t['Umix'])/rho
    mq=(t['q0']*t['p1']-t['q1']*t['p0'])/t['nQ']
    con=np.sqrt(np.maximum(t['charge_RMS']**2-t['nQ']**2,0))/abs(t['nQ'])
    peri=argrelextrema(mr2,np.less,order=50)[0]
    cyc=[]
    for i in range(len(peri)-1):
        s=slice(peri[i],peri[i+1])
        cyc.append(dict(tau0=float(tau[peri[i]]),tau1=float(tau[peri[i+1]]),fluct_fraction=float(dX[s].mean()/X2[s].mean()),
                        contrast_mean=float(con[s].mean()),contrast_max=float(con[s].max()),mean_r2_peri=float(mr2[peri[i]])))
    fc=[c for c in cyc if c['contrast_mean']>=1.]
    last=tau>=tau[peri[-1]]
    w=(tau>=700)&(tau<=800)
    row=dict(pericentre_taus_first8=[float(x) for x in tau[peri[:8]]],
        mean_r2_pericentre_first=float(mr2[peri[0]]),mean_r2_apocentre_first=float(mr2[argrelextrema(mr2,np.greater,order=50)[0][0]]),
        instantaneous_X_fluct_fraction_first_half=float(tau[np.argmax(dX/X2>=.5)]),
        cycle_containing_first_half=[c for c in cyc if c['tau0']<=tau[np.argmax(dX/X2>=.5)]<c['tau1']][:1],
        first_cycle_with_mean_contrast_ge_1=fc[0] if fc else None,
        max_cycle_avg_fluct_fraction=float(max(c['fluct_fraction'] for c in cyc)),
        cycles_with_cycle_avg_fluct_fraction_ge_half=[c for c in cyc if c['fluct_fraction']>=.5],
        final_values=dict(outside=float(outside[-1]),phi_plus_mix=float(pm[-1]),meanQ_over_Q=float(mq[-1]),fluct_fraction=float(dX[-1]/X2[-1])),
        last_partial_orbit=dict(tau0=float(tau[peri[-1]]),outside_min=float(outside[last].min()),outside_max=float(outside[last].max()),
            phi_plus_mix_min=float(pm[last].min()),phi_plus_mix_max=float(pm[last].max()),meanQ_min=float(mq[last].min()),meanQ_max=float(mq[last].max())),
        window_700_800_mean=dict(outside=float(outside[w].mean()),phi_plus_mix=float(pm[w].mean()),meanQ_over_Q=float(mq[w].mean()),fluct_fraction=float(dX[w].mean()/X2[w].mean())),
        cycles=cyc)
    snap=sorted(glob.glob(d+'/*_s*.bin'))[-1];N,L,s,q,p=read_snapshot(snap);V=L**3
    k2=fd4_symbol(N,L);qh=np.fft.fftn(q,axes=(1,2,3),norm='forward');ph=np.fft.fftn(p,axes=(1,2,3),norm='forward')
    n=np.rint(np.fft.fftfreq(N)*N);kphys=2*np.pi/L*np.sqrt(n[:,None,None]**2+n[None,:,None]**2+n[None,None,:]**2)
    nz=k2>0;Qt=float(np.mean(q[0]*p[1]-q[1]*p[0]));Qm=np.real(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0])
    r2=q[0]**2+q[1]**2;Xm=float(q[0].mean()**2+q[1].mean()**2);Xa=float(r2.mean());ph2=float(np.mean(q[2]**2))
    proj={}
    for name,m2 in [('massless',0.),('isotropic_Hartree_Codex',2*LX*Xa+G2*ph2),('Hartree_without_zero_mode',2*LX*(Xa-Xm)+G2*ph2),
                    ('tangential',LX*Xa+G2*ph2),('radial',3*LX*Xa+G2*ph2)]:
        om=np.sqrt(np.where(nz,k2+m2,1.));nr=(abs(ph[0])**2+abs(ph[1])**2+om**2*(abs(qh[0])**2+abs(qh[1])**2))/(2*om)
        nm=np.where(nz,(nr-Qm)/2,0.)
        proj[name]=dict(mass2=m2,nminus_over_net=float(nm.sum()/Qt),k_below_1p8=float(nm[kphys<1.8].sum()/Qt),k_below_1=float(nm[kphys<1.].sum()/Qt),
                        modes_k_below_1=int(((kphys<1.)&nz).sum()))
    seeded=None
    pj=f'{d}/chunk_state.json'
    if os.path.exists(pj):seeded=json.load(open(pj))['params']['nonzero_Fourier_modes']
    elif os.path.exists(d+'/parameters.json'):seeded=json.load(open(d+'/parameters.json'))['nonzero_Fourier_modes']
    row.update(projections=proj,X2_zero_mode_fraction=Xm/Xa,net_current_density=Qt,
        seeded_half_quanta_nminus_per_net=(.5*seeded/(Qt*V)) if seeded else None)
    out[tag]=row
    print(tag,json.dumps({k:v for k,v in row.items() if k not in('cycles',)},default=float)[:1800],'\n')
json.dump(out,open(outp,'w'),indent=1,default=float)

# ---- appended: equal-time 'virial' projection floor (own code). For each mode the quadratic form
# (|p|^2+w^2|q|^2)/(2w) is minimised at w^2=|p|^2/|q|^2; per-shell common w(k) minimises the shell sum.
def floors(tag):
    d=os.path.realpath(f'{work}/{tag}');snap=sorted(glob.glob(d+'/*_s*.bin'))[-1];N,L,s,q,p=read_snapshot(snap)
    k2=fd4_symbol(N,L);qh=np.fft.fftn(q,axes=(1,2,3),norm='forward');ph=np.fft.fftn(p,axes=(1,2,3),norm='forward')
    nz=k2>0;Qt=float(np.mean(q[0]*p[1]-q[1]*p[0]));Qm=np.real(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0])
    A=abs(ph[0])**2+abs(ph[1])**2;B=abs(qh[0])**2+abs(qh[1])**2
    permode=float(np.sum(((np.sqrt(A*B)-Qm)/2)[nz])/Qt)
    key=np.round(k2,9);u,inv=np.unique(key[nz],return_inverse=True)
    As=np.bincount(inv,weights=A[nz]);Bs=np.bincount(inv,weights=B[nz]);Qs=np.bincount(inv,weights=Qm[nz])
    pershell=float(np.sum((np.sqrt(As*Bs)-Qs)/2)/Qt)
    weff2=As/Bs;keff=np.sqrt(u)
    sel=[np.argmin(abs(keff-k)) for k in [.52,1.05,2.,3.,4.4,6.,8.]]
    return dict(per_mode_floor_nminus_over_net=permode,per_FD4_shell_floor_nminus_over_net=pershell,
        effective_omega2_minus_k2_at=dict((f'{keff[i]:.2f}',float(weff2[i]-u[i])) for i in sel))
fl={t:floors(t) for t in tags}
res=json.load(open(outp))
for t in tags:res[t]['virial_projection_floor']=fl[t];print(t,fl[t])
json.dump(res,open(outp,'w'),indent=1,default=float)
