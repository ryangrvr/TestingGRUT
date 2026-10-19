#!/usr/bin/env python3
"""Matched first/second Born graphs for the original rotating scalar profile.

All results are dimensionless K1/(g*v). No soliton energy is inferred from these
graphs alone; the nonperturbative subtracted spectrum is a separate calculation.
"""
from __future__ import annotations
from pathlib import Path
import hashlib,json,platform,sys,time
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
REF=HERE.parent/'grut_candidate_1_20261018/code/results_profile.npz'
NAMES=['A','B','C','D','E']
LOOP=16*np.pi**2


def potentials(F,S):
    return np.column_stack([3*(S*S-1)+F*F,(S*S-1)+.2*F*F,.1*F*F,2*F*S,.01*(S*S-1)])


def weights(r):
    w=np.ones(len(r));w[1:-1:2]=4;w[2:-1:2]=2
    return w*(r[1]-r[0])/3


def transform(r,V,k):
    w=weights(r);out=np.empty((len(k),5))
    for start in range(0,len(k),256):
        wave=np.sinc(np.outer(k[start:start+256],r)/np.pi)
        out[start:start+256]=wave@(4*np.pi*r[:,None]**2*V*w[:,None])
    return out


def kernels(k,O,n,mu=1.):
    z,w=leggauss(n);z=(z+1)/2;w=w/2
    specs=[(2.,2.,0.),(1.,1.,0.),(1.,1.,4*O*O),(2.,1.,O*O),(1.01,1.01,0.)]
    out=np.empty((len(k),5));mins=[]
    for i,(a,b,freq2) in enumerate(specs):
        den=z*a+(1-z)*b-z*(1-z)*freq2+np.outer(k*k,z*(1-z))
        mins.append(float(np.min(den)))
        if np.any(den<=0):raise ValueError('A restoration bubble crossed its pair cut')
        out[:,i]=-np.log(den/(mu*mu))@w
    return out,mins


def first_graph(integrals,mu=1.):
    def A(t):return t*(np.log(t/(mu*mu))-1)/LOOP
    parts=np.array([.5*A(2)*integrals[0],A(1)*integrals[1],0.,0.,.5*A(1.01)*integrals[4]])
    return parts


def second_graph(k,Vhat,kernel):
    factors=np.array([1,2,2,2,1])
    return -simpson(k[:,None]**2*Vhat**2*kernel*factors/(2*np.pi**2),x=k,axis=0)/(64*np.pi**2)


def main():
    started=time.perf_counter();pb=(HERE/'graphs_protocol.json').read_bytes();protocol=json.loads(pb)
    raw=REF.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    if sha!='0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3':raise ValueError('Changed reference profile')
    p=np.load(REF);O=float(p['Omega']);spline=CubicSpline(p['rho'],p['y'],axis=1)
    k=np.linspace(0,80,3201);runs=[];saved={};last=None
    kers={n:kernels(k,O,n) for n in [32,64,128]}
    for grid in protocol['radial_grids']:
        r=np.linspace(0,grid['R'],grid['N']);F,_,S,_=spline(r);V=potentials(F,S)
        integrals=4*np.pi*simpson(r[:,None]**2*V,x=r,axis=0)
        square=4*np.pi*simpson(r[:,None]**2*V**2,x=r,axis=0)
        Vhat=transform(r,V,k);one=first_graph(integrals)
        for kg in protocol['momentum_grids']:
            count=kg['N'];assert k[count-1]==kg['K']
            for n in protocol['Feynman_quadratures']:
                parts=second_graph(k[:count],Vhat[:count],kers[n][0][:count])
                parseval=simpson(k[:count,None]**2*Vhat[:count]**2/(2*np.pi**2),x=k[:count],axis=0)
                runs.append({'radial_R':grid['R'],'radial_N':grid['N'],'Kmax':kg['K'],'momentum_N':count,'Feynman_n':n,
                             'first_parts':dict(zip(NAMES,map(float,one))),'first_total':float(np.sum(one)),
                             'second_parts':dict(zip(NAMES,map(float,parts))),'second_total':float(np.sum(parts)),
                             'restored_graph_total':float(np.sum(one)+np.sum(parts)),
                             'Parseval_relative_errors':dict(zip(NAMES,map(float,np.abs(parseval/square-1))))})
        saved['r_N'+str(grid['N'])]=r;saved['V_N'+str(grid['N'])]=V;saved['Vhat_N'+str(grid['N'])]=Vhat
        last=(r,V,integrals,square,Vhat)
    r,V,integrals,square,Vhat=last;A,B,C,D,E=V.T
    delta=np.zeros((len(r),4,4));delta[:,0,0]=A;delta[:,1,1]=B+C;delta[:,2,2]=B-C;delta[:,3,3]=E
    delta[:,0,1]=delta[:,1,0]=D
    eigen=np.linalg.eigvalsh(delta);norm=np.max(abs(eigen),axis=1)
    I2=float(square@[1,2,2,2,1]);Irad3=float(simpson(r**3*np.trace(delta@delta,axis1=1,axis2=2),x=r))
    matrices={'sampled_V_operator_norm_max':float(norm.max()),'rho_at_sampled_norm_max':float(r[np.argmax(norm)]),'I2_integral_full_matrix':I2,'I_rad3_full_matrix':Irad3,
              'sampled_profile_bound_certified':False,
              'trace_delta_squared_via_full_matrix':float(4*np.pi*simpson(r*r*np.trace(delta@delta,axis1=1,axis2=2),x=r))}
    mus=[]
    for mu in [.5,1.,2.]:
        kernel,_=kernels(k,O,128,mu);g1=first_graph(integrals,mu);g2=second_graph(k,Vhat,kernel)
        mus.append({'mu_over_gv':mu,'graph_total':float(np.sum(g1)+np.sum(g2))})
    derivative=(mus[2]['graph_total']-mus[0]['graph_total'])/np.log(4)
    mV=2*integrals[0]+2*integrals[1]+1.01*integrals[4]
    exact_derivative=-(2*mV+I2)/(32*np.pi**2)
    finest=runs[-1]
    by=lambda N,K,n:next(x for x in runs if x['radial_N']==N and x['Kmax']==K and x['Feynman_n']==n)
    radial_difference=abs(by(10001,80,128)['second_total']-by(5001,80,128)['second_total'])
    momentum_difference=abs(by(10001,80,128)['second_total']-by(10001,40,128)['second_total'])
    feynman_difference=abs(by(10001,80,128)['second_total']-by(10001,80,64)['second_total'])
    checks={'profile_input_matches':sha=='0417c30367b52d6fdac021b8fb113d8486547faddf06c332617805d584a461c3',
            'all_five_Parseval_errors_small':max(finest['Parseval_relative_errors'].values())<1e-7,
            'radial_refinement':radial_difference<1e-4,'momentum_refinement':momentum_difference<1e-4,
            'Feynman_refinement':feynman_difference<1e-5,'all_bubbles_below_cuts':min(kers[128][1])>0,
            'vacuum_potentials_vanish':np.max(abs(potentials(np.zeros(2),np.ones(2))))==0,
            'local_counterterm_scale_derivative':abs(derivative-exact_derivative)<1e-6,
            'trace_matrix_normalization':abs(matrices['trace_delta_squared_via_full_matrix']-I2)<1e-9}
    result={'protocol_sha256':hashlib.sha256(pb).hexdigest(),'profile_sha256':sha,'Omega':O,
            'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()},
            'runs':runs,'selected_finest':finest,'potential_diagnostics':matrices,'scale_runs':mus,
            'explicit_scale_derivative':float(derivative),'local_counterterm_expected_derivative':float(exact_derivative),
            'refinement_differences':{'radial':float(radial_difference),'momentum':float(momentum_difference),'Feynman':float(feynman_difference)},
            'checks':{k:bool(v) for k,v in checks.items()},'all_checks_pass':all(checks.values()),
            'elapsed_seconds':time.perf_counter()-started,'scope':'Only restored first/second Born graphs. Not the full renormalized soliton energy.'}
    (HERE/'graphs_results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(HERE/'graphs_fourier_data.npz',k=k,**saved)
    (HERE/'graphs_execution.jsonl').write_text(json.dumps({'command':'python graphs.py','checks':result['checks'],'elapsed_seconds':result['elapsed_seconds']})+'\n')
    print(json.dumps({'selected':finest,'potential_diagnostics':matrices,'checks':result['checks'],'refinements':result['refinement_differences'],'elapsed_seconds':result['elapsed_seconds']},indent=2))
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
