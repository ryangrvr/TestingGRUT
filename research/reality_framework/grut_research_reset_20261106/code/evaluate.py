#!/usr/bin/env python3
"""Independent, small checks for the research reset."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_legendre

ROOT=Path(__file__).resolve().parents[1]
TRACKS=ROOT.parent/'grut_working_theory_v2_20261101/tracks'
BOUNDARIES={'at035':'early_boundary','at055':'reference','late':'late'}
MASS=1.
FREQUENCY=2.4
S0=2891.2
CRITICAL_DENSITY=1.87834e-29/1.7826619216e-24


def load_json(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_slice(root, cases, output):
    arrays={}
    provenance={}
    for tag,name in BOUNDARIES.items():
        source=root/'production'/name
        kinetic=root/'kinetic'/tag
        release=root/'abundance'/tag/'release_work/engine/release/release_x11_h55.json'
        history=source/'history.txt'
        spectrum=source/'spectra.txt'
        reference=cases['kinetic'][tag]
        found={k:sha(p) for k,p in (('history',history),('spectra',spectrum))}
        if any(reference[f'source_{k}_sha256']!=v for k,v in found.items()):
            raise ValueError(f'Raw {tag} source differs from pinned evidence')
        h=np.genfromtxt(history,names=True)
        zero=int(np.searchsorted(h['t'],0))
        indices=np.unique(np.r_[zero-1,zero,np.flatnonzero(h['t']>=h['t'][-1]-.02)])
        arrays[f'{tag}_source']=h[indices]
        arrays[f'{tag}_spectra']=np.genfromtxt(spectrum,names=True)
        ledger=kinetic/'native_ledger.txt'
        arrays[f'{tag}_kinetic']=np.genfromtxt(ledger,names=True)[[0,-1]]
        arrays[f'{tag}_charge']=np.array(load_json(release)['result']['N_X'])
        provenance[tag]=dict(source_sha256=found,kinetic_ledger_sha256=sha(ledger),
                             release_sha256=sha(release),full_source_hashes_matched=True)
    np.savez_compressed(output/'raw_slices.npz',**arrays)
    (output/'raw_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


def audit_boundary(data, tag, masses, abundance):
    source=data[f'{tag}_source']
    history=source[-1]
    spectra=data[f'{tag}_spectra']
    kinetic=data[f'{tag}_kinetic']
    a=float(history['a'])
    if (not math.isclose(float(kinetic['t'][0]),float(history['t']),abs_tol=1e-13) or
        not math.isclose(float(kinetic['a'][0]),a,abs_tol=1e-13)):
        raise ValueError(f'Kinetic boundary {tag} does not match its source')
    p=spectra['q']/a
    dq=2*min(spectra['q'])
    lo,hi=spectra['q']-dq/2,spectra['q']+dq/2
    w=dq*(lo*lo+lo*hi+hi*hi)/(6*math.pi**2)
    measure_error=float(np.max(abs(w/spectra['w']-1)))
    f=spectra['f_free']
    moments=[]
    moment_error=0.
    for i,M in enumerate(masses,1):
        select=spectra['species']==i
        hminus=select&(spectra['h']==-1)
        hplus=select&(spectra['h']==1)
        if not np.array_equal(spectra['q'][hminus],spectra['q'][hplus]):
            raise ValueError('Helicity quadratures do not match')
        P=p[select]
        E=np.sqrt(M*M+P*P)
        population=w[select]*f[select]/a**3
        vals=np.array([sum(population),sum(E*population),sum(P*P*population/(3*E))])
        expected=np.array([history[f'n_{i}'],history[f'rho_{i}'],history[f'p_free_{i}']])
        moment_error=max(moment_error,float(np.max(abs(vals/expected-1))))
        moments.append(vals)
    rows=np.array(moments)
    pressure=sum(rows[:,2])+(history['rho_dr']+history['rho_sm'])/3
    work=sum(source[f'work_{i}'] for i in range(1,5))/source['rho']
    work=work[source['t']>=history['t']-.02]
    total_pressure=float((pressure-history['P'])/history['rho'])
    energy_jump=float(kinetic['rho'][0]/history['rho']-1)
    # Same scale-factor normalization as the release. Recompute visible
    # entropy from the raw kinetic endpoint, without counting dark radiation.
    a0=np.interp(0,source['t'],source['a'])
    T=(30*kinetic['rho_sm'][-1]/(math.pi**2*106.75))**.25
    entropy=2*math.pi**2/45*106.75*(kinetic['a'][-1]*T/a0)**3
    charge=abs(float(data[f'{tag}_charge']))
    Y=charge/entropy
    omega=Y*abundance['population']['M_central']*abundance['X_mass_GeV']*S0/CRITICAL_DENSITY
    return dict(time=float(history['t']),helicity_count=2,mode_range=[float(min(f)),float(max(f))],
                shell_measure_error=measure_error,moment_error=moment_error,
                energy_jump=energy_jump,
                fermion_pressure_jump=float((sum(rows[:,2])-sum(history[f'p_{i}'] for i in range(1,5)))/
                                            history['rho']),
                complete_pressure_jump=total_pressure,
                kinetic_pressure_readout_error=float((pressure-kinetic['pressure'][0])/history['rho']),
                largest_sampled_work_last_window=float(max(abs(work))),
                number_by_species=rows[:,0].tolist(),energy_by_species=rows[:,1].tolist(),
                pressure_by_species=rows[:,2].tolist(),
                Y_charge=float(Y),conditional_omega_h2=float(omega),
                conditional_omega_replay_error=float(omega/abundance['conditional_Omega_h2']-1),
                energy_only_gate=abs(energy_jump)<1e-6,
                full_handoff_gate=bool(max(abs(energy_jump),abs(total_pressure),max(abs(work)))<1e-6))


def rotating_mass():
    M,mu=.7,.8
    results={}
    X=np.array([[0.,1.],[1.,0.]])
    Z=np.diag([1.,-1.])
    times=np.linspace(0,4,161)
    for p in (-.9,.9):
        energy,axes=np.linalg.eigh(M*X+p*Z)
        E=energy[1]

        def rhs(t,psi):
            vertex=M*np.exp(2j*mu*t)
            return -1j*np.array([[p,vertex],[vertex.conjugate(),-p]])@psi

        sol=solve_ivp(rhs,(0,4),axes[:,0].astype(complex),method='DOP853',
                      t_eval=times,rtol=2e-13,atol=2e-15)
        if not sol.success:
            raise ValueError('Direct rotating-mass reference did not converge')
        # Extract the instantaneous physical positive level in the lab frame.
        unrotate=sol.y*np.exp(np.array([-1j,1j])[:,None]*mu*times)
        measured=abs((axes.conj().T@unrotate)[1])**2
        effective=math.hypot(M,p+mu)
        exact=(M*mu/(E*effective))**2*np.sin(effective*times)**2
        r_base=-.5*np.array([M,0.,p])/E
        trans=np.array([0.,math.sqrt(.75),0.])
        states=np.stack((r_base+trans,r_base-trans))
        n=.5*(1+states@np.array([M,0.,p])/E)
        pressure=p*(states[:,2]+p/E)/6
        work=-M*mu*states[:,1]
        results[str(p)]=dict(max_probability_error=float(max(abs(exact-measured))),
                             norm_error=float(max(abs(np.sum(abs(sol.y)**2,axis=0)-1))),
                             max_probability=float(max(exact)),
                             static_probability=0.,
                             matched_number=n.tolist(),matched_energy=(E*n).tolist(),
                             matched_pressure=pressure.tolist(),opposite_work=work.tolist())
    return results


def pulse_transform(nu,tau):
    return (tau*math.sqrt(2*math.pi)/2 *
            (np.exp(-.5*(tau*(nu-FREQUENCY))**2)+np.exp(-.5*(tau*(nu+FREQUENCY))**2)))


def invariant_work(k,tau):
    threshold=math.sqrt(k*k+4*MASS*MASS)
    end=max(threshold,FREQUENCY)+14/tau
    def integral(u):
        nu=threshold+u*u
        # Use the difference of squares to resolve the physical threshold.
        excess=u*u*(nu+threshold)
        s=4*MASS*MASS+excess
        return 2*u*nu*pulse_transform(nu,tau)**2*math.sqrt(s*excess)
    value=quad(integral,0,math.sqrt(end-threshold),epsabs=0,epsrel=2e-11)[0]
    # The spatial mean-square of cos(kz) is 1/2. The local control
    # uses this same factor even when evaluating the k=0 response.
    return MASS*MASS*value/(32*math.pi**2)


def laboratory_work(k,tau,nodes):
    x,w=roots_legendre(nodes)
    p,wp=4*(x[:,None]+1),4*w[:,None]
    u,wu=x[None,:],w[None,:]
    E1=np.sqrt(MASS*MASS+p*p)
    E2=np.sqrt(MASS*MASS+p*p+k*k-2*k*p*u)
    energy=E1+E2
    invariant=energy*energy-k*k
    integral=np.sum(wp*wu*p*p*energy*invariant*pulse_transform(energy,tau)**2/(E1*E2))
    return float(MASS*MASS*integral/(32*math.pi**2))


def spatial_response():
    rows={}
    for tau in (12,24):
        local=invariant_work(0.,tau)
        for k in (1.,1.6):
            invariant=invariant_work(k,tau)
            first=laboratory_work(k,tau,256)
            second=laboratory_work(k,tau,512)
            rows[f'{tau}:{k}']=dict(k=k,tau=tau,threshold=math.hypot(k,2*MASS),
                                    invariant_energy=invariant,laboratory_energy=second,
                                    laboratory_refinement_error=abs(first/second-1),
                                    relative_error=abs(second/invariant-1),
                                    local_energy=local,physical_over_local=invariant/local)
    # Per squared coupling g^2=M^2/f^2 for a parent at rest. The
    # spatial amplitude is held at the same RMS value as the waves.
    widths={}
    expected=FREQUENCY*math.sqrt(1-4*MASS*MASS/FREQUENCY**2)/(16*math.pi)
    for tau in (24,48):
        integrated_parent_energy=FREQUENCY**2*math.sqrt(math.pi)*tau/4
        rate=invariant_work(0.,tau)/(MASS*MASS*integrated_parent_energy)
        widths[str(tau)]=dict(rate_per_g2=rate,expected=expected,relative_error=abs(rate/expected-1),
                             broken_double_counting_relative_error=abs(2*rate/expected-1))
    return rows,widths


def pressure_controls():
    a=1.7
    dx=2e-5
    q=np.array([.4,1.2,2.8])
    wf=np.array([.1,.3,.2])
    M=.7

    def comoving(loga,m):
        return float(np.sum(wf*np.sqrt(m*m+q*q*np.exp(-2*loga))))

    def pressure(m):
        p=q/a
        return float(np.sum(wf*p*p/np.sqrt(m*m+p*p))/(3*a**3))

    x=math.log(a)
    fd=-(comoving(x+dx,M)-comoving(x-dx,M))/(6*dx*a**3)
    exact=pressure(M)
    E0=comoving(x,0.)
    rad_ratio=pressure(0.)/(E0/a**3)
    rad_const=abs(2*a*comoving(math.log(2*a),0.)/(a*E0)-1)
    g=.4
    v=g/a
    w=math.hypot(M,v)

    def coupled(loga,rx,rz):
        v=g*math.exp(-loga)
        return .5*(v*rx+M*rz+math.hypot(M,v))

    fd_coupled=-(coupled(x+dx,1.,0.)-coupled(x-dx,1.,0.))/(6*dx*a**3)
    exact_coupled=(v+v*v/w)/(6*a**3)
    omitted=v*v/(6*a**3*w)
    vacuum_e=coupled(x,-v/w,-M/w)
    vacuum_p=- (coupled(x+dx,-v/w,-M/w)-coupled(x-dx,-v/w,-M/w))/(6*dx*a**3)
    return dict(volume_derivative_error=abs(fd/exact-1),massless_pressure_ratio=rad_ratio,
                massless_redshift_error=rad_const,interaction_pressure_error=abs(fd_coupled/exact_coupled-1),
                broken_pressure_error=abs(omitted/exact_coupled-1),vacuum_energy=vacuum_e,vacuum_pressure=vacuum_p)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    parser.add_argument('--raw-root',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    cases=load_json(TRACKS/'x11_h55/results/continuation.json')['cases']
    if args.raw_root:
        raw_slice(args.raw_root,cases,args.output)
    data_path=args.output/'raw_slices.npz' if args.raw_root else ROOT/'results/raw_slices.npz'
    masses=[r['M'] for r in load_json(TRACKS/'x2_strip/results/lower_H_lattice.json')['run']['channels']]
    with np.load(data_path,allow_pickle=False) as d:
        boundaries={tag:audit_boundary(d,tag,masses,cases['abundance'][tag]) for tag in BOUNDARIES}
    ref=boundaries['at055']['conditional_omega_h2']
    relative={k:v['conditional_omega_h2']/ref-1 for k,v in boundaries.items()}
    rotating=rotating_mass()
    spatial,widths=spatial_response()
    pressure=pressure_controls()
    tests={
        'rotating_solution':all(r['max_probability_error']<1e-9 and r['norm_error']<1e-9 for r in rotating.values()),
        'moments_leave_work_underdetermined':all(abs(r['matched_number'][1]-r['matched_number'][0])<1e-12 and
                                              abs(r['matched_pressure'][1]-r['matched_pressure'][0])<1e-12 and
                                              r['opposite_work'][0]==-r['opposite_work'][1]!=0 for r in rotating.values()),
        'independent_spatial_quadrature':all(max(r['relative_error'],r['laboratory_refinement_error'])<1e-6
                                            for r in spatial.values()),
        'rest_width_and_wrong_identical_factor':widths['48']['relative_error']<.005 and
            widths['48']['relative_error']<widths['24']['relative_error'] and
            widths['48']['broken_double_counting_relative_error']>.9,
        'resolve_spatial_threshold':spatial['24:1.6']['physical_over_local']<1e-6 and
                                     spatial['24:1.0']['physical_over_local']>.1,
        'pressure_and_expansion':max(pressure['volume_derivative_error'],pressure['interaction_pressure_error'],
                                     pressure['massless_redshift_error'],abs(pressure['massless_pressure_ratio']-1/3),
                                     abs(pressure['vacuum_energy']),abs(pressure['vacuum_pressure']))<1e-8 and
                                 pressure['broken_pressure_error']>.1,
        'raw_moments_and_abundance_arithmetic':all(max(r['shell_measure_error'],r['moment_error'],
                                                      abs(r['kinetic_pressure_readout_error']),
                                                      abs(r['conditional_omega_replay_error']))<1e-10
                                                for r in boundaries.values()),
        'reject_energy_only_splice':all(r['energy_only_gate'] and not r['full_handoff_gate']
                                        for r in boundaries.values()) and abs(relative['late'])>.01,
    }
    result=dict(rotating=rotating,spatial=spatial,rest_width=widths,pressure=pressure,
                boundaries=boundaries,moved_boundary_relative_abundance=relative,
                raw_slice_sha256=sha(data_path),tests=tests,observational_evidence=False)
    (args.output/'checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(tests,indent=2))
    raise SystemExit(not all(tests.values()))


if __name__=='__main__':
    main()
