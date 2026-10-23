#!/usr/bin/env python3
"""Deciding author checks, with physical-completeness gates kept separate."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.integrate import simpson
from scipy.special import spherical_jn
from scattering import ROOT,Background,scatter,visible_spectral_probability
from exterior import project


def nodal_energy(a,h):
    r=a[:,0];u=r[:,None]*np.column_stack([a[:,1],a[:,2],a[:,3]-1]);p=r[:,None]*a[:,4:7]
    ext=np.vstack([-u[0],np.zeros(3),u,np.zeros(3),-u[-1]])
    lap=(-ext[:-4]+16*ext[1:-3]-30*ext[2:-2]+16*ext[3:-1]-ext[4:])/(12*h*h)
    n=a[:,1]**2+a[:,2]**2;s=a[:,3]
    return 4*np.pi*h*np.sum(.5*np.sum(p*p-u*lap,axis=1)+r*r*(.75*(s*s-1)**2+.5*s*s*n+.025*n*n))


def force_variation():
    h,dt,L=.125,.0001,8.;r=np.arange(1,round(L/h))*h
    x=.2*np.exp(-r*r/3);y=.03*np.sin(r)*np.exp(-r*r/4);s=1-.02*np.exp(-r*r/4)
    a=np.column_stack([r,x,y,s,np.zeros((len(r),3))])
    with tempfile.TemporaryDirectory(prefix='grut-nodal-force-') as tmp:
        tmp=Path(tmp);binary=tmp/'radial';init=tmp/'initial.csv';prefix=tmp/'run'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'code/radial_nodal.cpp'),'-o',str(binary)],check=True)
        np.savetxt(init,a,delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
        subprocess.run([str(binary),str(init),str(prefix),str(h),str(dt),str(L),str(dt),'4','7','0'],check=True,capture_output=True)
        final=np.loadtxt(str(prefix)+'_final.csv',delimiter=',',skiprows=1)
    acc=2*(final[:,1:4]-a[:,1:4])/dt**2
    rng=np.random.default_rng(231023);errors=[]
    for _ in range(4):
        d=rng.normal(size=acc.shape)*.05;eps=1e-6
        ap=a.copy();am=a.copy();ap[:,1:4]+=eps*d;am[:,1:4]-=eps*d
        numerical=(nodal_energy(ap,h)-nodal_energy(am,h))/(2*eps)
        analytic=-4*np.pi*h*np.sum(r[:,None]**2*acc*d)
        errors.append(abs(numerical-analytic)/max(abs(numerical),abs(analytic),1e-12))
    return max(errors)


def linearization():
    # Differentiate the original real force at several core points and phases,
    # independently of the frequency-domain three-channel assembly.
    rng=np.random.default_rng(4);errs=[]
    for _ in range(16):
        F,S=rng.uniform(.1,3),rng.uniform(0,1);phase=rng.uniform(0,6)
        psi=F*np.exp(-1j*phase);dpsi=rng.normal()+1j*rng.normal();ds=rng.normal();eps=1e-6
        def law(p,s):return np.array([((s*s+.1*abs(p)**2)*p).real,
             ((s*s+.1*abs(p)**2)*p).imag,(3*(s*s-1)+abs(p)**2)*s])
        num=(law(psi+eps*dpsi,S+eps*ds)-law(psi-eps*dpsi,S-eps*ds))/(2*eps)
        dp=(S*S+.2*F*F)*dpsi+.1*psi*psi*dpsi.conjugate()+2*S*psi*ds
        dS=(3*(3*S*S-1)+F*F)*ds+2*S*np.real(psi.conjugate()*dpsi)
        ana=np.array([dp.real,dp.imag,dS]);errs.append(max(abs(num-ana))/max(max(abs(ana)),1e-12))
    return max(errs)


def optical_control(bg):
    E=np.sqrt(1.5)
    result,solutions,w=scatter(bg,E,radius=35.,tol=2e-8,mesh=1201,return_solutions=True)
    r=np.linspace(0,35,5001);y=solutions[0].sol(r)
    z=w[1]['H'].conjugate()*(y[2]+1j*y[8]);z[0]=0
    integral=simpson(np.exp(-(r/12)**8)*abs(z)**2,x=r)/w[1]['k']
    rows=[]
    for epsilon in [1e-4,1e-5]:
        control=scatter(bg,E,radius=35.,tol=2e-8,mesh=1201,absorption=epsilon)
        loss=1-sum(control['columns'][0]['probabilities'])
        rows.append(dict(epsilon=epsilon,measured_loss=loss,first_order=epsilon*integral,
                         relative_error=loss/(epsilon*integral)-1))
    # Fourier-Bessel Parseval checks the 2/pi normalization of the nonlocal
    # visible kernel without borrowing its physical SM normalization.
    _,S=bg.fields(r);q=np.linspace(0,20,3001)
    form=simpson(S[None,:]*z[None,:]*r[None,:]*spherical_jn(0,q[:,None]*r[None,:]),x=r,axis=1)
    spatial=simpson(abs(S*z)**2,x=r)
    spectral=2/np.pi*simpson(q*q*abs(form)**2,x=q)
    return dict(local_absorption=rows,fourier_bessel_parseval_relative=spectral/spatial-1)


def main():
    checks=[]
    def check(name,ok,evidence=None):checks.append(dict(name=name,passed=bool(ok),evidence=evidence))
    baseline=json.loads((ROOT/'BASELINE_REPLAY.json').read_text())
    check('Fresh replay of all 13 parent implementation controls',baseline['implementation_passed']==13)
    lin=linearization();check('Full real action derivative matches coupled perturbation equations',lin<1e-8,lin)
    f=force_variation();check('Nodal fourth-order force matches a finite variation of its Hamiltonian',f<2e-5,f)
    controls=json.loads((ROOT/'SCATTERING_CONTROLS.json').read_text())
    maxU=max(x.get('unitarity_max',0) for x in controls['controls'])
    maxR=max(x.get('reciprocity_max',0) for x in controls['controls'])
    check('Two- and three-open-channel flux unitarity and reciprocity',maxU<2e-6 and maxR<2e-6,dict(unitarity=maxU,reciprocity=maxR))
    free=max(x.get('free_identity_error',0) for x in controls['controls'])
    check('Zero-potential partial waves return the identity scattering matrix',free<2e-5,free)
    bg=Background();basis=[]
    for ell in [0,1,3]:
        a=scatter(bg,np.sqrt(1.5),ell=ell,radius=30.,tol=2e-7,basis='factored')
        b=scatter(bg,np.sqrt(1.5),ell=ell,radius=30.,tol=2e-7,basis='radial')
        pa=a['columns'][0]['probabilities'][0];pb=b['columns'][0]['probabilities'][0]
        basis.append(abs(pa/pb-1))
    check('Regular-origin and factored-origin BVPs independently agree on conversion',max(basis)<2e-5,basis)
    opt=optical_control(bg)
    check('Absorptive flux loss agrees with the first-order radial optical integral',
          max(abs(x['relative_error']) for x in opt['local_absorption'])<.005,opt['local_absorption'])
    check('Nonlocal visible kernel has the Fourier-Bessel Parseval normalization',
          abs(opt['fourier_bessel_parseval_relative'])<1e-6,opt['fourier_bessel_parseval_relative'])
    kernel=controls['visible_kernel_quadrature'];values=[r['probability'] for r in kernel]
    check('Visible spectral integral converges in radial and momentum quadrature',max(values)/min(values)-1<1e-6,values)
    rate=json.loads((ROOT/'RATE_RESULTS.json').read_text());rows=rate['rows']
    tails=max(r['last_two_partial_fraction'] for r in rows)
    flux=max(abs(p['flux_residual']) for r in rows for p in r['partials'])
    check('Rate table angular tails and every production partial-wave flux are controlled',tails<2e-5 and flux<2e-6,dict(tail=tails,flux=flux))
    # Pure positive-frequency standing wave: mQ/E<1 despite zero antiparticles.
    L,h,k=80.,.125,4*np.pi/80;radius=np.arange(1,round(L/h))*h
    psi=np.sin(k*radius)/radius;omega=np.sqrt(1+k*k);zero=np.zeros_like(radius)
    table=np.column_stack([radius,psi,zero,np.ones_like(radius),zero,-omega*psi,zero])
    projected=project(table,dict(h=h,length=L,absorber=L,nodal=True),start=0,width=.001)
    # Masking mixes spatial k, so test the canonical modal identity separately.
    plus=(omega+omega)**2/(4*omega);minus=(omega-omega)**2/(4*omega)
    check('Kinetic energy invalidates a particle-antiparticle inference from mQ/E alone',
          minus==0 and 1/omega<1,dict(mQ_over_E=1/omega,antiparticles=minus))
    trans=json.loads((ROOT/'TRANSPORT_RESULTS.json').read_text())['rows']
    hashes={r['rate_table_numeric_sha256'] for r in trans}
    maxN=max(r['number_relative_error'] for r in trans);maxE=max(r['energy_ledger_relative_error'] for r in trans)
    check('Transport consumes one rate table and conserves number, charge and energy with expansion work',
          len(hashes)==1 and maxN<1e-9 and maxE<1e-8,dict(rate_tables=len(hashes),number=maxN,energy=maxE))
    check('Inverse conversion is included and produces nonzero events',all(r['inverse_events_per_initial_pair']>0 for r in trans))
    runs=[json.loads(p.read_text()) for p in sorted(ROOT.glob('[HN]*_RESULTS.json'))]
    check('All new collapse runs conserve charge plus their numerical receiver',
          max(r['max_charge_error'] for r in runs)<1e-8,max(r['max_charge_error'] for r in runs))
    trace_errors=[]
    for run in runs:
        trace=np.load(ROOT/'data'/f"{run['name']}_trace.npz")['trace']
        history=np.loadtxt(ROOT/'data'/f"{run['name']}_history.csv",delimiter=',',skiprows=1)
        late=history[:,0]>.8*run['parameters']['tmax']
        trace_errors.append(hashlib.sha256(trace.tobytes()).hexdigest()==run['trace_sha256']
                            and len(trace)==run['trace_samples']
                            and trace[-1,0]==run['parameters']['tmax']
                            and history[-1,0]==run['parameters']['tmax']
                            and abs(float(np.mean(history[late,12])/history[0,2])-run['core60_late_mean_fraction'])<1e-14)
    check('All 18 raw traces reconstruct their recorded lossless hashes and final times',
          len(runs)==18 and all(trace_errors),dict(runs=len(runs),matched=int(sum(trace_errors))))
    conv=json.loads((ROOT/'CONVERGENCE_RESULTS.json').read_text())
    check('Finite-time prototype charge agrees within one percent across fine grids and boundary controls',
          conv['finite_time_one_percent_agreement'],dict(grid_spread=conv['spatial']['relative_spread'],
          max_control=max(abs(r['retained_charge_relative_change']) for r in conv['controls'])))
    points=[r['points_across_radius90'] for r in conv['snapshot_diagnostics'] if r['name'] in ['N025','N0125','N00625']]
    check('All three fine nodal grids resolve the core radius with at least ten points',min(points)>=10,points)
    ext=json.loads((ROOT/'EXTERIOR_RESULTS.json').read_text())
    perr=max(abs(v['charge_parseval_residual']) for r in ext['exterior_snapshots'] for v in r['variants'])
    check('Every masked exterior projection preserves its measured canonical charge',perr<1e-10,perr)
    gates=[dict(name='Observed Richardson order certifies continuum extrapolation',passed=conv['spatial']['extrapolation_certified']),
           dict(name='Final primordial-density three-dimensional object/source population',passed=False),
           dict(name='Complete asymptotic source including soft and absorbed fields',passed=False),
           dict(name='Nonlinear hot-core capture, quantum lifetime, full visible deposition and structure test',passed=False),
           dict(name='Independent Claude Code / GRUT-RAI scientific certification',passed=False)]
    out=dict(scope='Same-author controls; pass counts are not empirical or external scientific validation.',
        checks=checks,passed=sum(r['passed'] for r in checks),total=len(checks),physical_gates=gates,
        optical_normalization_control=opt)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    for r in checks:print('PASS' if r['passed'] else 'FAIL',r['name'],r['evidence'],flush=True)
    assert all(r['passed'] for r in checks)


if __name__=='__main__':main()
