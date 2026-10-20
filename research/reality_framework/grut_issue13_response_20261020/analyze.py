#!/usr/bin/env python3
"""Retrospective threshold audit and field diagnostics, without importing old analysis.

Requires numpy and scipy. Reads original frozen data plus read-only reruns.
It does not fit a new benchmark, calibrate vacuum energy, or calculate a CMB
likelihood. Reproducible numerical checks are distinct from scientific evidence.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq

P = Path(__file__).resolve().parent
OLD = P.parent/'grut_3d_empirical_20261020'
AUDIT = P.parent/'grut_empirical_audit_20261020'
C = 299792.458
CHECKS = {}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def check(name, result):
    CHECKS[name] = bool(result)

def readcsv(path):
    return np.atleast_1d(np.genfromtxt(path,delimiter=',',names=True))

def cosmology():
    # Physical density held fixed; Hfloor is a derived unit, not an observed H0.
    ms = 1000 * min(.1,(.01*.001)**.25)
    wb, mb = .02233-3*.00036, .938*1.01
    wf = wb*(1+ms/mb)
    hf = 100*math.sqrt(wf)
    H100 = 100/3.0856775814913673e19*6.582119569e-25
    rho100 = 3*(2.435323e18)**2*H100**2
    table = json.loads((AUDIT/'maser_table1.json').read_text())['rows']
    quad_errors = []
    def da(z,L=0,K=0,order=None):
        # L and K are physical squared-rate terms in (km/s/Mpc)^2.
        def invH(x):
            return 1/np.sqrt(hf*hf*(1+x)**3 + L + K*(1+x)**2)
        if order:
            x,w = np.polynomial.legendre.leggauss(order)
            chi = C*z/2 * np.dot(w,invH((x+1)*z/2))
        else:
            val,err = quad(invH,0,z,epsabs=1e-14,epsrel=1e-12)
            quad_errors.append(C*err/(1+z))
            chi = C*val
        if K < 0:
            radius=C/math.sqrt(-K)
            # All roots explored here are safely before the sine maximum.
            if chi/radius >= math.pi/2:
                raise RuntimeError('Cannot use monotone closed-distance envelope here')
            dm=radius*math.sin(chi/radius)
        elif K > 0:
            radius=C/math.sqrt(K); dm=radius*math.sinh(chi/radius)
        else:
            dm=chi
        return dm/(1+z)
    result={'mstar_GeV':ms,'omega_floor':wf,'Hfloor_km_s_Mpc':hf,
            'rho_crit_100_GeV4':rho100,'corners':{}}
    for label,vpec in [('primary',-750.),('stress',-1500.)]:
        rows=[]
        for row in table:
            zobs=(row['v_km_s']+3*row['v_sigma_km_s'])/C
            beta=vpec/C
            z=(1+zobs)/math.sqrt((1+beta)/(1-beta))-1
            rows.append({'name':row['name'],'z':z,
                         'Dlow_Mpc':row['D_Mpc']-3*row['D_minus_Mpc']})
        def worst(L=0,K=0):
            return max(r['Dlow_Mpc']-da(r['z'],L,K) for r in rows)
        L=brentq(lambda x:worst(L=x),-.99*hf*hf,0,xtol=1e-9)
        K=brentq(lambda x:worst(K=x),-.99*hf*hf,0,xtol=1e-9)
        h0=math.sqrt(hf*hf+K)
        Om=hf*hf/(h0*h0)
        mattig=[]
        for r in rows:
            z=r['z']
            # Independent closed dust FRW expression, rather than quadrature.
            dmattig=2*C/(h0*Om**2*(1+z)**2)*(Om*z+(Om-2)*(math.sqrt(1+Om*z)-1))
            mattig.append(abs(dmattig-da(z,K=K)))
            r.update({'flat_DAmax_Mpc':da(z), 'DA_at_vacuum_threshold_Mpc':da(z,L=L),
                      'DA_at_curvature_threshold_Mpc':da(z,K=K)})
            check(f'{label}_{r["name"]}_quad_vs_Gauss',
                  abs(da(z,L=L,K=0)-da(z,L=L,K=0,order=32))<1e-9)
        critical=max(rows,key=lambda r:r['Dlow_Mpc']-r['DA_at_curvature_threshold_Mpc'])['name']
        result['corners'][label]={
            'vpec_km_s':vpec,'rows':rows,'critical_row':critical,
            'vacuum_threshold':{'L_km_s_Mpc_squared':L,'omega_L':L/10000,
                'rho_GeV4':L/10000*rho100,'abs_rho_fourth_root_meV':(-L/10000*rho100)**.25*1e12,
                'H0_floor_saturated_km_s_Mpc':math.sqrt(hf*hf+L)},
            'curvature_threshold':{'K_km_s_Mpc_squared':K,'omega_k':K/10000,
                'Omega_k_floor_saturated':K/(hf*hf+K),
                'H0_floor_saturated_km_s_Mpc':h0,'mattig_max_absolute_error_Mpc':max(mattig)}}
        check(label+'_closed_quadrature_vs_Mattig',max(mattig)<1e-8)
        check(label+'_vacuum_threshold_root',abs(worst(L=L))<1e-8)
        check(label+'_curvature_threshold_root',abs(worst(K=K))<1e-8)
    check('primary_curvature_review_reproduced',abs(result['corners']['primary']['curvature_threshold']['Omega_k_floor_saturated']+.5630192053576689)<1e-10)
    for z in [.001,.01,.04,1]:
        exact=2*C/hf/(1+z)*(1-1/math.sqrt(1+z))
        check(f'flat_limit_{z}',abs(da(z)-exact)<1e-9)
    result['maximum_QUADPACK_error_estimate_Mpc']=max(quad_errors)
    result['interpretation']='Necessary envelope sensitivity, separately varying K or constant vacuum. Equality is a boundary, not a strict mismatch. No viable full history or parameter calibration is established.'
    return result

def plateau_counts(F,S):
    result={}
    for a in [.5,.75,.9]:
        for x in [.01,.1]:
            result[f'absS_gt_{a}_F2_lt_{x}']={
                'negative':int(np.count_nonzero((S < -a)&(F<x))),
                'positive':int(np.count_nonzero((S > a)&(F<x)))}
    return result

def archived():
    out=[]
    raw=OLD/'raw'
    manifest=json.loads((OLD/'RAW_MANIFEST.json').read_text())['raw_files']
    check('all_190_original_raw_files_sha256',all(sha(raw/n)==v['sha256'] for n,v in manifest.items()))
    for path in sorted(raw.glob('*.csv')):
        if '.components.' in path.name:
            continue
        a=readcsv(path);idx=int(np.argmin(a['minS']))
        item={'name':path.stem,'samples':len(a),'final_t':float(a['t'][-1]),
              'minS_over_saved_full_grid_diagnostics':float(a['minS'][idx]),
              'time_of_saved_minS':float(a['t'][idx]),'endpoint_full_grid_minS':float(a['minS'][-1]),
              'negative_saved_times':int(np.count_nonzero(a['minS']<0))}
        meta=raw/(path.name+'.spatial.json')
        if meta.exists():
            m=json.loads(meta.read_text())
            A=np.fromfile(raw/(path.name+'.spatial.f32'),dtype='<f4').reshape(m['n'],m['n'],m['nz'],3)
            F,S=A[...,0],A[...,1]
            neg=S<0
            item['endpoint']={'source_grid':[m['N'],m['N'],m['Nz']],
                'retained_grid':list(S.shape),'sample_spacing':m['stride']*m['h'],
                'negative_cells':int(neg.sum()),'total_retained_cells':S.size,
                'minS':float(S.min()),'maxS':float(S.max()),
                'negative_F2_quantiles_0_25_50_75_100':np.quantile(F[neg],[0,.25,.5,.75,1]).tolist() if neg.any() else [],
                'plateaus':plateau_counts(F,S)}
            check(path.stem+'_subsample_min_consistent',float(S.min())>=float(a['minS'][-1])-2e-7)
        out.append(item)
    return out

def reruns():
    result=[]
    for name in json.loads((P/'PROTOCOL.json').read_text())['diagnostic_runs']:
        a=readcsv(P/'runs'/f'{name}.csv.walls.csv')
        old=readcsv(OLD/'raw'/f'{name}.csv')
        new=readcsv(P/'runs'/f'{name}.csv')
        check(name+'_scalar_columns_identical',old.dtype.names==new.dtype.names)
        check(name+'_scalar_times_identical',np.array_equal(old['t'],new['t']))
        errors={n:float(np.max(np.abs(old[n]-new[n])/(1+np.abs(old[n])))) for n in old.dtype.names}
        check(name+'_trajectory_replay',max(errors.values())<1e-9)
        times=np.searchsorted(a['t'],old['t'])
        check(name+'_every_step_times_match_old_output',np.max(np.abs(a['t'][times]-old['t']))<1e-10)
        check(name+'_minS_matches_old_saved_values',np.max(np.abs(a['minS'][times]-old['minS']))<1e-9)
        case=json.loads((P/'runs'/f'{name}.execution.json').read_text())['case']
        dt=.05;T=case.get('T',160)
        check(name+'_all_drift_steps_scanned',len(a)==round(T/dt)+1 and np.max(np.abs(np.diff(a['t'])-dt))<1e-10)
        count=int(a['negative'].max());i=int(a['minS'].argmin())
        maxima={n:int(a[n].max()) for n in ['n_lt_m10','n_lt_m25','n_lt_m50','n_lt_m75','n_lt_m90','negative_lowX']}
        f=P/'runs'/f'{name}.csv.spatial.f32';g=OLD/'raw'/f'{name}.csv.spatial.f32'
        x,y=np.fromfile(f,dtype='<f4'),np.fromfile(g,dtype='<f4')
        field_error=float(np.max(np.abs(x.astype(float)-y)))
        check(name+'_endpoint_replay',field_error<1e-5)
        result.append({'name':name,'h':case['h'],'dt':dt,'T':T,'scanned_states':len(a),
            'minS':float(a['minS'][i]),'time_of_minS':float(a['t'][i]),
            'maxS_over_run':float(a['maxS'].max()),'max_negative_cells':count,
            'maxima_of_threshold_counts':maxima,
            'first_negative_time':float(a['t'][np.flatnonzero(a['negative']>0)[0]]) if count else None,
            'minimum_meanF2_on_negative_cells':float(a['meanF2_negative'][a['negative']>0].min()) if count else None,
            'max_normalized_scalar_replay_error':max(errors.values()),
            'max_absolute_endpoint_float32_replay_error':field_error,
            'maximum_relative_E_drift':float(np.max(np.abs(new['E']/new['E'][0]-1))),
            'maximum_relative_Q_drift':float(np.max(np.abs(new['Q']/new['Q'][0]-1)))})
    return result

def wall_scales():
    v,lam,g=1000.,.01,.1
    analytic=2*math.sqrt(2*lam)/3*v**3
    val,err=quad(lambda p:math.sqrt(lam/2)*(v*v-p*p),-v,v,epsabs=1e-5,epsrel=1e-12)
    check('wall_tension_integral',abs(val/analytic-1)<1e-13)
    coeff=lam/4+g*g/6+.01**2/12
    check('thermal_trace_Hessian',(3*lam+2*g*g+.01**2)/12==coeff)
    # Rough ZKO-style dimensional diagnostic; the coefficient is not a C_l calculation.
    # Conventional non-reduced G=1/(8 pi Mred^2); H0=70 is an external illustration.
    M=2.435323e18
    H=70/3.0856775814913673e19*6.582119569e-25
    G=1/(8*math.pi*M*M)
    rough=1e-5*H/G
    return {'tree_planar_X_D_zero_tension_GeV3':analytic,'tension_cube_root_GeV':analytic**(1/3),
            'quadrature_error_estimate_GeV3':err,
            'tree_wall_width_mX_inverse_units':g*math.sqrt(2/lam),
            'zero_chemical_potential_highT_phi_mass_coefficient':coeff,
            'leading_curvature_zero_T_hidden_over_v':math.sqrt(lam/coeff),
            'rough_isotropy_sigma_GeV3':rough,'ratio_to_rough_isotropy_sigma':analytic/rough,
            'rough_isotropy_assumptions':'Persistent horizon-scale scaling network; G sigma/H0 of order DeltaT/T, coefficient set to one, H0=70 km/s/Mpc, DeltaT/T=1e-5. No Candidate-1 spectrum or likelihood. This diagnostic is not applicable to transient microscopic bubbles.',
            'thermal_assumptions':'Common equilibrated hidden scalar bath, high-temperature expansion, zero chemical potential; curvature zero is not a fully resummed transition temperature or a statement about the visible bath.'}

def controls():
    F=np.zeros((8,8,8));S=np.ones_like(F);S[:4]=-1
    p=plateau_counts(F,S)
    check('synthetic_two_vacua_counted',all(v=={'negative':256,'positive':256} for v in p.values()))
    S[:4]=-.2;F[:4]=2
    check('synthetic_shallow_charged_patch_not_negative_vacuum',all(v['negative']==0 for v in plateau_counts(F,S).values()))
    check('shallow_patch_still_has_sign_change',bool((S<0).any() and (S>0).any()))

def main():
    controls()
    output={'cosmology':cosmology(),'archive_diagnostics':archived(),
            'every_step_reruns':reruns(),'wall_scales':wall_scales(),
            'environment':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
                           'platform':platform.platform()},
            'source_sha256':sha(Path(__file__)), 'protocol_sha256':sha(P/'PROTOCOL.json')}
    output['checks']=CHECKS
    output['check_count']=len(CHECKS)
    output['all_checks_pass']=all(CHECKS.values())
    (P/'RESULTS.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    with (P/'ARCHIVED_MINIMA.csv').open('w') as f:
        keys=['name','samples','final_t','minS_over_saved_full_grid_diagnostics','time_of_saved_minS','negative_saved_times']
        writer=csv.DictWriter(f,keys,extrasaction='ignore',lineterminator='\n');writer.writeheader();writer.writerows(output['archive_diagnostics'])
    print(json.dumps({'all_checks_pass':output['all_checks_pass'],'checks':len(CHECKS),
                      'failed':[k for k,v in CHECKS.items() if not v],
                      'cosmology':output['cosmology']['corners'],
                      'reruns':output['every_step_reruns'],'wall_scales':output['wall_scales']},indent=2))
    if not output['all_checks_pass']:
        raise SystemExit(1)

if __name__=='__main__':
    main()
