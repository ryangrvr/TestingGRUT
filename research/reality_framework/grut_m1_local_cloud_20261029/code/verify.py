"""Independent algebraic and same-author numerical controls for local plume."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import numpy as np
from scipy.integrate import simpson,quad
from scipy.optimize import brentq
from scipy.signal import correlate,correlation_lags
from scipy.special import roots_legendre
from neutral_return import return_probability, MH
from local_rates import finite_width_control
from daughter_cloud import root_rates
from cross_source import kernels

ROOT=Path(__file__).resolve().parents[1]

def read(tag):return json.loads((ROOT/(tag+'.json')).read_text())


def synthetic_flux():
    t=np.arange(0,1000.0001,.25);rad=(2400+.5)*.025;E=2.;k=np.sqrt(E*E-1);amp=.01
    z=amp/rad*np.exp(1j*(k*rad-E*t));dz=(1j*k-1/rad)*z;p=-1j*E*z
    trace=np.column_stack([t,z.real,z.imag,0*t,dz.real,dz.imag,0*t,p.real,p.imag,0*t])
    path=ROOT/'.work'/'synthetic_trace.csv';path.parent.mkdir(exist_ok=True)
    np.savetxt(path,trace,delimiter=',',header='tau,x,y,s_minus1,dx,dy,ds,px,py,ps',comments='')
    cmd=[sys.executable,str(ROOT/'code'/'surface_source.py'),'--trace',str(path),'--width','400','--tag','synthetic_source']
    subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
    d=read('synthetic_source');expected=4*np.pi/.01*k*amp**2*len(t)*.25
    return {'expected_charge':expected,'actual_charge':d['raw_flux']['charge'],
            'relative_charge_error':d['raw_flux']['charge']/expected-1,
            'plus_signed_relative_error':d['species']['plus']['signed_number']/expected-1,
            'anti_fraction':d['species']['minus']['positive_projection_number']/expected,
            'apparent_incoming_fraction_k_gt_0_1':d['species']['plus']['split_k_gt_0.1']['in']/expected}


def analytic_overlap():
    # Constant beams of duration T: their time correlation is T-|delta|.
    T=1000.;dt=1.;n=int(T/dt);F=np.ones(n);A=.75;r0=60.;R=800.
    c=correlate(F,F,mode='full',method='fft')*dt;lags=correlation_lags(n,n,'full')*dt
    u=np.linspace(0,np.log(R/r0),2001);r=r0*np.exp(u)
    num=simpson(np.interp(A*(r-r0),lags,c,left=0,right=0)/r,x=u)
    stop=min(R,r0+T/A)
    exact=(T+A*r0)*(1/r0-1/stop)-A*np.log(stop/r0)
    return {'analytic':exact,'numeric':float(num),'relative_error':float(num/exact-1)}


def root_control():
    errors=[];rows=[]
    for Ed,mu in ((1.1,-1.),(1.4,-.5),(1.4,0.),(1.4,.5),(2.5,.9),(2.5,-.8)):
        pd=np.sqrt(Ed*Ed-1);B=pd*mu;D=Ed*Ed-B*B;disc=4-D
        derived=[]
        if disc>=0:
            for sign in (1,-1):
                E=(2*Ed+sign*B*np.sqrt(disc))/D;p=np.sqrt(max(E*E-1,0))
                if 1<E<8 and abs(Ed*E-B*p-2)<1e-8 and not any(abs(E-x)<1e-8 for x in derived):derived.append(E)
        f=lambda E:Ed*E-B*np.sqrt(max(E*E-1,0))-2
        grid=np.linspace(1,8,10001);values=np.array([f(x) for x in grid]);numeric=[]
        for i in range(len(grid)-1):
            if values[i]*values[i+1]<0:numeric.append(brentq(f,grid[i],grid[i+1]))
        derived=sorted(derived);numeric=sorted(numeric)
        assert len(derived)==len(numeric),(Ed,mu,derived,numeric)
        err=max([abs(a-b) for a,b in zip(derived,numeric)] or [0.]);errors.append(err)
        rows.append({'Ed':Ed,'mu':mu,'roots':numeric,'absolute_error':err})
    return {'max_absolute_error':max(errors),'cases':rows}


def return_control():
    z=np.load(ROOT/'source_t400.npz');n=z['neutral_positive'].sum(axis=0);e=z['energy_edges'];R=float(z['probe_radius_mX_inverse'])
    x,w=roots_legendre(32);x=(x+1)/2;w=w/2;expected=0.
    for i,N in enumerate(n):
        lo=max(MH,e[i]);hi=e[i+1]
        if N<=1e-12 or hi<=lo:continue
        expected+=N*sum(v*return_probability(lo+(hi-lo)*u,R,R) for u,v in zip(x,w))
    rows=[]
    for tag in ('daughter_262144_s19','daughter_262144_s37','daughter_262144_s59'):
        d=read(tag)['result'];rows.append({'tag':tag,'hits':[r['direct_matching_surface_hits_per_core'] for r in d['rows']]})
    return {'analytic_all_time_crossings_per_species':float(expected),'finite_time_QMC':rows}


def daughter_root_rate_control():
    rows=read('GENERAL_WIDTH_CONTROL')['rows'];tgrid=np.linspace(0,10000,101)
    egrid=np.linspace(1.0125,7.9875,280);flux=np.ones((len(tgrid),len(egrid)))
    errors=[]
    for case in rows:
        Ed=case['Ed'];mu=case['relative_cosine'];pd=np.sqrt(Ed*Ed-1);radius=100.
        got,_=root_rates(np.array([Ed]),np.array([pd]),np.array([[mu]]),np.array([[5000.]]),
                         np.array([[radius]]),60.,tgrid,egrid,flux)
        exact=0.
        for Ep in case['roots']:
            pp=np.sqrt(Ep*Ep-1);vp=pp/Ep;deriv=2*(Ed-pd*mu*Ep/pp)
            exact+=np.pi*.2**2/(2*Ed*Ep*abs(deriv))/(4*np.pi*radius**2*vp)
        errors.append(abs(float(got[0,0])/exact-1))
    return {'maximum_relative_error':max(errors),'cases':len(errors)}


def angular_kernel_control():
    # Independent integration in the angle, with explicit pole breakpoints.
    energy=np.array([1.1,1.4,1.8,2.5]);_,full=kernels(energy,128)
    rows=[]
    for i,j in ((0,0),(0,1),(1,2),(2,3)):
        a=energy[i];b=energy[j];pa=np.sqrt(a*a-1);pb=np.sqrt(b*b-1)
        f=lambda mu:.2**4*np.sqrt(max(0.,1-4/(2+2*(a*b-pa*pb*mu))))/(64*np.pi*a*b*((2+2*(a*b-pa*pb*mu)-6)**2+6*.00018756589919939715**2))
        mu0=(2+2*a*b-6)/(2*pa*pb);w=np.sqrt(6)*.00018756589919939715/(2*pa*pb)
        splits=[-1.,1.]
        splits.extend(float(mu0+v*w) for v in (-100,-10,-1,0,1,10,100) if -1<mu0+v*w<1)
        splits=sorted(set(splits))
        direct=sum(quad(f,x,y,epsabs=1e-12,epsrel=1e-9,limit=200)[0] for x,y in zip(splits[:-1],splits[1:]))
        rows.append({'energies':[float(a),float(b)],'angular_quad':direct,'transformed_quad':float(full[i,j]),
                     'relative_difference':float(full[i,j]/direct-1)})
    return {'max_relative_error':max(abs(r['relative_difference']) for r in rows),'rows':rows}


def main():
    checks=[]
    def add(name,passed,**details):checks.append({'name':name,'pass':bool(passed),**details})
    synthetic=synthetic_flux();add('monochromatic_outgoing_flux',abs(synthetic['relative_charge_error'])<1e-12 and
        abs(synthetic['plus_signed_relative_error'])<.002 and synthetic['anti_fraction']<.002,**synthetic)
    overlap=analytic_overlap();add('ballistic_rectangular_pulse_analytic_overlap',abs(overlap['relative_error'])<1e-4,**overlap)
    roots=root_control();add('general_resonance_roots_against_brent',roots['max_absolute_error']<1e-9,**roots)
    daughter_roots=daughter_root_rate_control();add('daughter_root_rate_against_independent_roots',daughter_roots['maximum_relative_error']<1e-9,**daughter_roots)
    for tag in ('source_t200','source_t400','source_t800','source_h005_r40','source_h005_r60','source_h005_r80'):
        d=read(tag);err=max(abs(x) for x in d['parseval_error'].values())
        a=max(d['species'][s]['negative_bin_weight']/d['species'][s]['positive_projection_number'] for s in ('plus','minus','neutral'))
        add('time_frequency_partition_'+tag,err<1e-7 and a<.001,max_absolute_parseval_error=err,
            max_signed_negative_fraction=a,interpretation='Small signed-bin weight does not prove an incoherent quantum source.')
    n=read('NEUTRAL_RETURN')['checks'];add('neutral_decay_kinematics',n['solid_angle_and_independent_root_max_absolute_error']<1e-11 and
        n['two_daughter_energy_and_radial_momentum_max_absolute_error_mX']<1e-11 and
        n['analytic_inward_energy_upper_bound_mX']<n['min_positive_reverse_exchange_threshold_mX'],**n)
    rc=return_control();mean=np.mean([v for r in rc['finite_time_QMC'] for v in r['hits']]);relative=abs(mean/rc['analytic_all_time_crossings_per_species']-1)
    add('QMC_return_vs_analytic_all_time',relative<.08,relative_difference=relative,**rc,
        caveat='Finite horizon versus all-time, windowed neutral source and uniform within-bin energy.')
    q64=read('continuous_e280_q64')['result']['expected_unattenuated_events'];q128=read('continuous_e280_q128')['result']['expected_unattenuated_events']
    add('continuous_energy_quadrature',abs(q64/q128-1)<.003,gauss64=q64,gauss128=q128)
    r32=read('continuous_e560_q32')['result']['expected_unattenuated_events'];r64=read('continuous_e560_q64')['result']['expected_unattenuated_events']
    add('fine_energy_quadrature',abs(r32/r64-1)<.003,gauss32=r32,gauss64=r64)
    g=read('ANGULAR_GEOMETRY');add('straight_ray_geometric_cone',g['tests']['max_sampled_sin_angle_above_geometric_bound']<1e-12 and
        g['tests']['support_fraction_monotonically_nonincreasing'],**g['tests'])
    ck=angular_kernel_control();add('cross_source_angular_kernel_against_direct_quad',ck['max_relative_error']<.002,**ck)
    fine=read('shot_noise_dr2p5')['rows'];coarse=read('shot_noise_dr5')['rows']
    numbererr=max(abs(row[key]) for row in fine for key in ('number_relative_error_plus','number_relative_error_minus'))
    momenterr=max(abs(a['poisson_mean_product_over_product_of_means']/b['poisson_mean_product_over_product_of_means']-1)
                  for a,b in zip(fine,coarse))
    add('poisson_profile_number_and_resolution',numbererr<.005 and momenterr<.01,
        max_number_relative_error=numbererr,max_cross_moment_relative_difference_dr2p5_vs_dr5=momenterr)
    c32=read('cross_e280_q32')['result']['rows'][1]['integrated'][-1]['events_per_reference_core']
    c64=read('cross_e280_q64')['result']['rows'][1]['integrated'][-1]['events_per_reference_core']
    add('poisson_cross_source_pole_quadrature',abs(c32/c64-1)<.001,gauss32=c32,gauss64=c64)
    closure=read('mixing_seed37')['regions'][0]['neutral_parents_and_daughters']['snapshots']
    enerr=max(x['per_parent_energy_closure_max_mX'] for x in closure)
    qerr=max(x['per_parent_charge_closure_max'] for x in closure)
    add('neutral_decay_sample_energy_and_charge',enerr<1e-12 and qerr==0,max_energy_error_mX=enerr,max_charge_error=qerr)
    f=finite_width_control();add('smooth_sample_nwa_vs_finite_width',max(abs(v['relative_difference']) for v in f)<.1,
        samples=f,interpretation='One sampled case differs by about 8%; no universal NWA bound follows.')
    out={'same_author':True,'independent_replication':False,'checks':checks,'all_engineering_checks_pass':all(x['pass'] for x in checks),
         'limitations':'These tests do not validate local angular distributions, wave-to-particle conversion, hidden channels, finite-time quantum coherence or the final cosmological mass budget.'}
    (ROOT/'VERIFICATION.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(out,indent=2))
    if not out['all_engineering_checks_pass']:raise SystemExit(1)

if __name__=='__main__':main()
