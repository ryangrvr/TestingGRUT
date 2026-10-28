"""Nonlinear incoming monopole packets, evolved with the unchanged reviewed engine."""
from common import *
from scipy.fft import dst,idst

CASES={
 'baseline01':dict(h=.1,dt=.0125,k=1.,Nin=0.,core=True),
 'baseline005':dict(h=.05,dt=.00625,k=1.,Nin=0.,core=True),
 'packet04':dict(h=.1,dt=.0125,k=.4,Nin=100.,core=True),
 'packet10':dict(h=.1,dt=.0125,k=1.,Nin=100.,core=True),
 'packet15':dict(h=.1,dt=.0125,k=1.5,Nin=100.,core=True),
 'packet10_small':dict(h=.1,dt=.0125,k=1.,Nin=10.,core=True),
 'packet10_phase':dict(h=.1,dt=.0125,k=1.,Nin=100.,core=True,phase=np.pi/2),
 'packet10_fine':dict(h=.05,dt=.00625,k=1.,Nin=100.,core=True),
 'free10':dict(h=.1,dt=.0125,k=1.,Nin=100.,core=False),
}

def run_packet(tag,par):
    src=DYN/'code/radial.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('radial_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off',str(src),'-o',str(exe)],check=True)
    h,dt,k,Nin,core=[par[n] for n in ('h','dt','k','Nin','core')]
    phase=par.get('phase',0.);L=300.;T=1000.;r0=80.;width=12.;n=round(L/h);r=(np.arange(n)+.5)*h
    vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    F,S=profile(r)
    if not core:F=np.zeros_like(F);S=np.ones_like(S)
    z=F.astype(complex);p=-1j*OMEGA*z
    if core:p*=2000/np.sum(-np.imag(np.conj(z)*p)*vol)
    # Odd radial extension; positive lab-time frequency gives negative X charge.
    g=lambda x:np.exp(-.5*((x-r0)/width)**2+1j*k*(x-r0)+1j*phase)
    u=g(r)-g(-r);modes=dst(u,type=2,norm='ortho')
    km=(np.arange(n)+1)*np.pi/L;klat=2*np.sin(km*h/2)/h
    energies=np.sqrt(1+klat**2)
    pu=1j*idst(energies*modes,type=2,norm='ortho')
    zw=u/r;pw=pu/r;charge=float(np.sum(-np.imag(np.conj(zw)*pw)*vol)/G2)
    scale=np.sqrt(Nin/(-charge)) if Nin else 0.
    zw*=scale;pw*=scale;modes*=scale
    z+=zw;p+=pw
    actualQ=float(np.sum(-np.imag(np.conj(z)*p)*vol)/G2)
    d=WORK/tag;d.mkdir(exist_ok=True);initial=d/'initial.csv'
    np.savetxt(initial,np.c_[r,z.real,z.imag,S,p.real,p.imag,np.zeros(n)],delimiter=',',
               header='rho,x,y,s,px,py,ps',comments='')
    input_hash=hashlib.sha256(initial.read_bytes()).hexdigest();record=d/'result.json'
    if record.exists():
        previous=json.loads(record.read_text())
        if previous['initial_sha256']==input_hash and previous['engine_sha256']==hashlib.sha256(src.read_bytes()).hexdigest():return previous
        raise RuntimeError('Stale packet result '+tag)
    specweight=4*np.pi*h*energies*abs(modes)**2/G2
    np.savez_compressed(d/'initial_spectrum.npz',k_effective=klat,energy=energies,negative_number_weight=specweight)
    start=time.time()
    with open(d/'run.log','w') as log:
        subprocess.run([str(exe),str(initial),str(d/'state'),str(h),str(dt),str(L),str(T),'40','180','.2'],check=True,stdout=log,stderr=subprocess.STDOUT)
    hist=np.genfromtxt(d/'state_history.csv',names=True,delimiter=',')
    e=hist['E']+hist['E_abs'];q=hist['Q']+hist['Q_abs']
    out={'tag':tag,'parameters':par,'L':L,'T':T,'r0':r0,'width':width,
         'engine_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'initial_sha256':input_hash,
         'actual_initial_Q':actualQ,'initial_packet_Q':float(np.sum(-np.imag(np.conj(zw)*pw)*vol)/G2),
         'initial_spectral_number':float(specweight.sum()),'seconds':time.time()-start,
         'energy_receiver_error':float(np.max(abs(e/e[0]-1))),
         'charge_receiver_absolute_error_in_quanta':float(np.max(abs(q-q[0]))/G2)}
    record.write_text(json.dumps(out,indent=2)+'\n');return out

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('cases',nargs='*');args=parser.parse_args()
    for tag in args.cases or list(CASES):print(json.dumps(run_packet(tag,CASES[tag])),flush=True)
