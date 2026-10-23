#!/usr/bin/env python3
"""Separate implementation checks from unresolved physical population gates."""
import json
import hashlib
from pathlib import Path
import subprocess
import tempfile
import numpy as np

from daughters import ROOT,MX,MS,MD,G,V,KAPPA,MH,VH,WIDTHS,CM3S,ann,decay_grid,stream,OR,OM,H0,C_MPC_S
from run_core import load_trace
from spectra import field_spectrum


def discrete_energy(a,h):
    n=len(a);edges=np.arange(n+1)*h
    vol=4*np.pi/3*np.diff(edges**3);face=4*np.pi*edges[1:]**2/h
    face[-1]*=2
    fields=a[:,1:4];vel=a[:,4:7]
    nextfields=np.vstack([fields[1:],[0,0,1]])
    grad=.5*np.sum(face[:,None]*(nextfields-fields)**2)
    x2=fields[:,0]**2+fields[:,1]**2;s=fields[:,2]
    pot=.75*(s*s-1)**2+.5*s*s*x2+.025*x2*x2
    return float(grad+np.sum(vol*(.5*np.sum(vel*vel,axis=1)+pot)))


def force_control():
    """C++ first drift acceleration against a finite energy variation."""
    h,dt,L=.25,.002,10.
    r=(np.arange(round(L/h))+.5)*h
    x=.1*np.exp(-r*r/8)*(1+.2*np.sin(r));y=.03*np.exp(-r*r/8)*np.sin(.7*r)
    s=1+.01*np.exp(-r*r/8)
    a=np.column_stack([r,x,y,s,np.zeros_like(r),-x,np.zeros_like(r)])
    with tempfile.TemporaryDirectory(prefix='grut-force-control-') as tmp:
        tmp=Path(tmp);binary=tmp/'radial';inp=tmp/'initial.csv';prefix=tmp/'run'
        np.savetxt(inp,a,delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/'code/radial.cpp'),'-o',str(binary)],check=True)
        subprocess.run([str(binary),str(inp),str(prefix),str(h),str(dt),str(L),str(dt),'5','8','0'],check=True,capture_output=True)
        final=np.loadtxt(str(prefix)+'_final.csv',delimiter=',',skiprows=1)
    force=2*(final[:,1:4]-a[:,1:4]-dt*a[:,4:7])/dt**2
    edges=np.arange(len(r)+1)*h;vol=4*np.pi/3*np.diff(edges**3)
    rng=np.random.default_rng(231023);errs=[]
    for _ in range(4):
        direction=rng.normal(size=(len(r),3))*.1
        plus=a.copy();minus=a.copy();eps=1e-6
        plus[:,1:4]+=eps*direction;minus[:,1:4]-=eps*direction
        variation=(discrete_energy(plus,h)-discrete_energy(minus,h))/(2*eps)
        from_force=-np.sum(vol[:,None]*force*direction)
        errs.append(abs(variation-from_force)/max(abs(variation),abs(from_force),1e-12))
    return max(errs)


def main():
    data=json.loads((ROOT/'DAUGHTER_RESULTS.json').read_text())
    core=json.loads((ROOT/'CORE_RESULTS.json').read_text())
    spectra=json.loads((ROOT/'SPECTRAL_RESULTS.json').read_text())
    checks=[]
    def check(name,passed,evidence=None):checks.append(dict(name=name,passed=bool(passed),evidence=evidence))
    # Direct spin-summed ee amplitude cross section, independently of the
    # spectral-width formula used by ann(). Include the same mixed poles.
    W=2*MX;me=.00051099895
    theta=KAPPA*V*VH/(MH*MH-MS*MS)
    point=data['point'];Gamma=point['Gamma_s_GeV']
    prop=theta*(1/(W*W-MS*MS+1j*MS*Gamma)-1/(W*W-MH*MH+1j*MH*.0041))
    direct=(2*G*G*V)**2*abs(prop)**2*(me/VH)**2/(4*np.pi)*(1-4*me*me/W**2)**1.5
    check('Cold ee annihilation agrees with direct invariant spin-sum formula',abs(ann(W,0)/direct-1)<2e-12,
          ann(W,0)/direct-1)
    check('Light-pole free-pair envelope remains below named Planck screen even at unit free fraction',
          max(r['ratio_to_named_Planck_screen_if_ALL_DM_free'] for r in data['annihilation']['light_resonance_envelope'])<1e-10)
    capture=data['capture']
    check('Integer-charge removal/addition energies agree with omega and block cold one-free-s emission',
          max(abs(capture[k]/capture['chemical_potential_from_omega_GeV']-1) for k in ['removal_energy_GeV','addition_energy_GeV'])<2e-6
          and not capture['cold_single_free_s_channel_open'] and capture['deficit_GeV']>.05)
    boost_errors=[]
    for gamma in [1.,10.,1000.]:
        for mass in [MX,MD]:
            E,u,w=decay_grid(gamma,mass)
            boost_errors+=[abs(sum(w)-1),abs(np.dot(w,E)/(gamma*MS/2)-1),
                           float(max(abs(E*E-mass*mass*(1+u*u)))/max(E*E))]
    check('Isotropic boosted daughters preserve normalization, mass shell and energy',max(boost_errors)<2e-12,max(boost_errors))
    check('Free-streaming numerical integral agrees with radiation-era closed form',
          abs(data['free_streaming']['radiation_analytic_control_relative'])<2e-9,
          data['free_streaming']['radiation_analytic_control_relative'])
    cold=[r for r in data['free_streaming']['rows'] if r['parent_gamma']==1]
    check('Later formation increases the calculated collisionless free-streaming scale',
          all(a['mean_free_stream_Mpc_today']<b['mean_free_stream_Mpc_today'] for a,b in zip(cold,cold[1:])))
    force_error=force_control()
    check('Full-field finite-volume forces equal the variation of the discrete energy',force_error<2e-5,force_error)
    runs=core['runs']
    check('Noether charge plus explicit absorber receiver conserves in every production run',
          max(r['max_total_charge_relative_error'] for r in runs)<5e-10,
          max(r['max_total_charge_relative_error'] for r in runs))
    check('Hamiltonian-plus-receiver energy error decreases with time/space refinement',
          next(r for r in runs if r['name']=='A_finer')['max_total_energy_relative_error']<4e-5
          and next(r for r in runs if r['name']=='A_time')['max_total_energy_relative_error']<1.5e-4)
    parseval=[];fold=[];qfold=[]
    for row in spectra['rows']:
        for key in ['X','scalar']:
            r=row[key]
            parseval.append(abs(r['total_windowed_net_energy']-r['time_domain_windowed_net_energy'])/max(abs(r['total_windowed_net_energy']),1e-20))
        r=row['conditional_particle_fold']
        fold.append(abs(r['scalar_to_daughters_energy_relative_residual']))
        qfold.append(abs(r['returned_net_charge']-r['direct_free_X_net_charge'])/max(abs(r['direct_free_X_net_charge']),1.))
    check('Complex/real Fourier flux agrees with direct time-domain Parseval accounting',max(parseval)<2e-12,max(parseval))
    check('Restricted free-s decay folding conserves energy and creates no net charge',max(fold)<2e-12 and max(qfold)<2e-12,
          dict(energy=max(fold),charge=max(qfold)))
    trace_errors=[]
    for r in runs:
        name=r['name'];trace=load_trace(r)
        assert len(trace)==r['trace_samples'] and np.isfinite(trace).all()
        assert hashlib.sha256(trace.tobytes()).hexdigest()==r['trace_float64_sha256']
        if name not in ['A_finer','A_remote']:continue
        fresh,_=field_spectrum(trace,r,'X')
        old=next(x for x in spectra['rows'] if x['run']==name)['X']
        trace_errors.append(abs(fresh['total_windowed_net_energy']/old['total_windowed_net_energy']-1))
    check('Lossless trace shards reconstruct exact arrays and regenerate selected fluxes',max(trace_errors)<2e-13,max(trace_errors))
    # Free outgoing wave normalization, including both complex charge signs
    # and the real-scalar factor relative to a complex field amplitude.
    test_errors=[]
    n=4096;dt=.125;t=np.arange(n)*dt;T=n*dt;test_run=dict(parameters=dict(probe=10.,h=.25))
    r=(round(10/.25-.5)+.5)*.25
    for kind,mass,mode in [('X',1.,1),('X',1.,-1),('s',np.sqrt(6),1)]:
        trace=np.zeros((n,10));trace[:,0]=t
        j=100 if kind=='X' else 256;omega=2*np.pi*j/T;k=np.sqrt(omega*omega-mass*mass)
        phase=-mode*omega*t;amplitude=.01/r
        if kind=='X':
            z=amplitude*np.exp(1j*phase);der=(1j*mode*k-1/r)*z;p=-1j*mode*omega*z
            trace[:,1]=z.real;trace[:,2]=z.imag;trace[:,4]=der.real;trace[:,5]=der.imag;trace[:,7]=p.real;trace[:,8]=p.imag
            expected_E=4*np.pi*r*r*amplitude**2*omega*k
            expected_Q=4*np.pi*r*r*amplitude**2*mode*k
        else:
            trace[:,3]=amplitude*np.cos(phase);trace[:,6]=-k*amplitude*np.sin(phase)-trace[:,3]/r
            trace[:,9]=omega*amplitude*np.sin(phase)
            expected_E=2*np.pi*r*r*amplitude**2*omega*k;expected_Q=0.
        result,_=field_spectrum(trace,test_run,kind,alpha=0)
        test_errors.append(abs(result['total_windowed_net_energy']/(expected_E*T)-1))
        if kind=='X':test_errors.append(abs(result['total_windowed_net_charge']/(expected_Q*T)-1))
    check('Free positive/negative-charge and real-field waves have the correct flux normalization',max(test_errors)<2e-12,max(test_errors))
    # These are physical completeness gates, deliberately not asserted green.
    a=[r for r in runs if r['name'] in ['A','A_fine','A_finer']]
    values=[r['late_mean_core_charge_fractions']['60'] for r in a]
    spread=(max(values)-min(values))/np.mean(values)
    coverage=max(r['X']['controlled_fraction_of_outgoing_linear_energy'] for r in spectra['rows'])
    gates=[dict(name='Late retained charge converged within one percent across the tested spatial grids',
                passed=bool(spread<.01),relative_spread=spread,values=values),
           dict(name='At least 95 percent of outgoing charged energy covered by the controlled linear particle band',
                passed=bool(coverage>=.95),best_measured_fraction=coverage),
           dict(name='Remaining capture/internal-excitation portal leakage rate computed',passed=False),
           dict(name='Nonthermal mixed-population Boltzmann and Lyman-alpha comparison computed',passed=False),
           dict(name='Actual primordial-density 3D final population produced',passed=False)]
    output=dict(scope='Implementation and scoped rate checks are separate from physical completeness; same author.',
                implementation_controls=checks,implementation_passed=sum(r['passed'] for r in checks),
                implementation_total=len(checks),physical_population_gates=gates)
    (ROOT/'VERIFICATION_RESULTS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    for r in checks:print('PASS' if r['passed'] else 'FAIL',r['name'],r['evidence'])
    for r in gates:print('MET' if r['passed'] else 'OPEN',r['name'])
    assert all(r['passed'] for r in checks),'Implementation control failed; do not promote the physical results.'


if __name__=='__main__':main()
