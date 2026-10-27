"""en_rerun INDEPENDENT: reference projections (massless, isotropic Hartree, per-shell virial floor) on every
retained N48 snapshot (Delta tau = 20) of the reproduced main history. Own code.
usage: python3 -I projection_history.py <N48_dir> <out_json>"""
import sys,glob,json,struct
import numpy as np
LX,G2=.001,.01
d,outp=sys.argv[1],sys.argv[2]
def sym(N,L):
    k=2*np.pi*np.fft.fftfreq(N);s=(30-32*np.cos(k)+2*np.cos(2*k))/12/(L/N)**2
    return s[:,None,None]+s[None,:,None]+s[None,None,:]
rows=[]
for f in sorted(glob.glob(d+'/state_s*.bin')):
    with open(f,'rb') as h:
        N,=struct.unpack('<i',h.read(4));L,=struct.unpack('<d',h.read(8));s,=struct.unpack('<d',h.read(8));h.read(4)
        a=np.fromfile(h,dtype='<f8').reshape(6,N,N,N)
    q,p=a[:3],a[3:];k2=sym(N,L);qh=np.fft.fftn(q,axes=(1,2,3),norm='forward');ph=np.fft.fftn(p,axes=(1,2,3),norm='forward')
    nz=k2>0;Qt=float(np.mean(q[0]*p[1]-q[1]*p[0]));Qm=np.real(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0])
    A=abs(ph[0])**2+abs(ph[1])**2;B=abs(qh[0])**2+abs(qh[1])**2
    n=np.rint(np.fft.fftfreq(N)*N);kp=2*np.pi/L*np.sqrt(n[:,None,None]**2+n[None,:,None]**2+n[None,None,:]**2)
    mh=2*LX*float(np.mean(q[0]**2+q[1]**2))+G2*float(np.mean(q[2]**2))
    r=dict(eta=s,tau=20+8.869218738806*s+.5*1.9857962750634526*s*s,hartree_mass2=mh)
    for name,m2 in [('massless',0.),('hartree',mh)]:
        w=np.sqrt(np.where(nz,k2+m2,1.));nm=np.where(nz,((A+w*w*B)/(2*w)-Qm)/2,0.)
        r[name]=float(nm.sum()/Qt);r[name+'_k_ge_1p8']=float(nm[kp>=1.8].sum()/Qt)
    u,inv=np.unique(np.round(k2[nz],9),return_inverse=True)
    As=np.bincount(inv,weights=A[nz]);Bs=np.bincount(inv,weights=B[nz]);Qs=np.bincount(inv,weights=Qm[nz])
    r['shell_floor']=float(np.sum((np.sqrt(As*Bs)-Qs)/2)/Qt);rows.append(r)
    print(f"tau={r['tau']:7.2f} massless={r['massless']:.3f} hartree={r['hartree']:.3f} floor={r['shell_floor']:.3f} m2H={mh:.2f} (k>=1.8: {r['massless_k_ge_1p8']:.3f}/{r['hartree_k_ge_1p8']:.3f})")
late=[r for r in rows if r['tau']>=599]
summ={k:dict(min=min(r[k] for r in late),max=max(r[k] for r in late),mean=float(np.mean([r[k] for r in late]))) for k in ['massless','hartree','shell_floor']}
json.dump(dict(rows=rows,late_600_800=summ),open(outp,'w'),indent=1);print(summ)
