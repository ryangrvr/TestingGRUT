"""A reproducible classical Fourier realization of the matched light marginal."""
from en_common import *

def initialize(N=48,L=12.,seed=1,cutoff=8.,noise_scale=1.):
    data=np.load(ROOT/f'MATCH_L{L:g}.npz');keys=data['keys'];lookup={int(v):i for i,v in enumerate(keys)}
    coeff=np.zeros((6,N,N,N),complex);nmax=int(np.ceil(cutoff*L/(2*np.pi)));modes=0
    for i in range(-nmax,nmax+1):
      for j in range(-nmax,nmax+1):
       for k in range(-nmax,nmax+1):
        n=(i,j,k);first=next((x for x in n if x),0)
        if first<=0:continue
        nn=i*i+j*j+k*k;mom=2*np.pi/L*np.sqrt(nn)
        if not 1.8<=mom<=cutoff:continue
        assert max(abs(i),abs(j),abs(k))<N/2
        ix=lookup[nn]
        rng=np.random.default_rng(np.random.SeedSequence([seed,int(L*100),i+128,j+128,k+128]))
        z=(rng.normal(size=14)+1j*rng.normal(size=14))/np.sqrt(2)
        full=data['evolved_scalar'][ix]@z[:12];ph=data['evolved_phi'][ix]@z[12:]
        light=np.r_[full[:2],ph[0],full[6:8],ph[1]]*noise_scale/np.sqrt(L**3)
        coeff[(slice(None),i%N,j%N,k%N)]=light
        coeff[(slice(None),(-i)%N,(-j)%N,(-k)%N)]=light.conjugate();modes+=2
    coeff[:,0,0,0]=np.r_[data['light_mean_q'],data['light_mean_p']]
    fields=np.fft.ifftn(coeff,axes=(-3,-2,-1),norm='forward')
    assert np.max(abs(fields.imag))<1e-10
    q,p=fields.real[:3].copy(),fields.real[3:].copy()
    meanq=data['light_mean_q'];meanp=data['light_mean_p'];bgQ=meanq[0]*meanp[1]-meanq[1]*meanp[0]
    Q=float(np.mean(q[0]*p[1]-q[1]*p[0]));E=energy(q,p,L)
    params=dict(N=N,L=L,seed=seed,upper_seed_cutoff=cutoff,lower_seed_cutoff=1.8,noise_scale=noise_scale,
        nonzero_Fourier_modes=modes,tau_match=20.,a_match=float(data['a_match']),C_match=float(data['C_match']),
        nQ=Q,rho=E,homogeneous_nQ=float(bgQ),sample_charge_offset=float(Q/bgQ-1),
        full_fields_seeded_linearly=['X','Y','Sigma','phi'],nonlinear_fields=['X','phi'],
        scope='Finite-band classical Gaussian realization; no charge/energy retuning. Omitted IR/UV quantum state, counterterms and later heavy/gauge feedback remain uncomputed.')
    return q,p,params

if __name__=='__main__':
    q,p,r=initialize();print(json.dumps(r,indent=2))
