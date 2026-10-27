"""en_neutral: OWN X-fluctuation spectrum of the reproduction N48 snapshots: where <|delta chi_X|^2> sits in lattice k.
usage: python3 -I xspec.py <rundir> <out_json>"""
import sys, struct, json, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import en_lib as E
from pathlib import Path
rd=Path(sys.argv[1]); N=48; L=12.0
nf=np.fft.fftfreq(N)*N; NX,NY,NZ=np.meshgrid(nf,nf,nf,indexing='ij'); klat=np.sqrt(E.lattice_k2(np.stack([NX,NY,NZ],-1),N,L))
a0=8.869218738806; C=1.9857962750634526
out={}
for path in sorted(rd.glob('state_s*.bin')):
    with open(path,'rb') as f:
        n_,L_,s_,bits=struct.unpack('<iddi',f.read(24)); a=np.fromfile(f,dtype='<f4').reshape(6,N,N,N).astype(float)
    tau=20+a0*s_+.5*C*s_*s_
    if min(abs(tau-np.array([60,70,80,100,150,200,400,800])))>1e-6: continue
    S=0; parts={}
    pw=sum(abs(np.fft.fftn(a[c])/N**3)**2 for c in (0,1)); pw[0,0,0]=0
    tot=pw.sum()
    out[f'{tau:.0f}']=dict(S=float(tot), frac_k_lt_1p8=float(pw[klat<1.8].sum()/tot), frac_1p8_4=float(pw[(klat>=1.8)&(klat<4)].sum()/tot),
        frac_4_8=float(pw[(klat>=4)&(klat<8)].sum()/tot), frac_gt_8=float(pw[klat>=8].sum()/tot))
    print(f'{tau:.0f}', out[f'{tau:.0f}'])
json.dump(out, open(sys.argv[2],'w'), indent=1)
