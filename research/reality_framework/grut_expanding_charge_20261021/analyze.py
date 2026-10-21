"""Summarize every declared run, including negative formation outcomes."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
RAW=HERE/'raw'

def load(name):
    d=np.atleast_1d(np.genfromtxt(RAW/(name+'.csv'),delimiter=',',names=True))
    top=[json.loads(s) for s in (RAW/(name+'.topology.jsonl')).read_text().splitlines() if s.strip()]
    cfg=json.loads((RAW/(name+'.config.json')).read_text())
    meta=json.loads((RAW/(name+'.meta.json')).read_text())
    return d,top,cfg,meta

def summary(name):
    d,top,cfg,meta=load(name)
    formed=(d['fplus_05']>=.1)&(d['fminus_05']>=.1)
    ix=np.flatnonzero(formed)
    form=None
    if len(ix):
        i=int(ix[0]);den=d['fplus_05'][i]+d['fminus_05'][i]
        form={'a_interval':[float(d['a'][max(0,i-1)]),float(d['a'][i])],
              'a_first_sample':float(d['a'][i]),'pplus_among_raw_vacuum_cells':float(d['fplus_05'][i]/den)}
    after=[t for t in top if form is not None and t['a']>=form['a_first_sample']-1e-9]
    end=d[-1]
    residual=d['H_conf']-d['H_conf'][0]-d['work_conf']
    # A conservative sampled denominator: cumulative ABSOLUTE work is >= |net work|.
    # Passing with this smaller denominator suffices for the protocol's larger one.
    norm=np.maximum.reduce([np.full(len(d),abs(d['H_conf'][0])),np.abs(d['H_conf']),np.abs(d['work_conf']),np.full(len(d),1e-300)])
    q0=float(d['Q_com'][0])
    qdrift=float(np.max(np.abs(d['Q_com']/q0-1)))
    budget=float(np.max(np.abs(residual)/norm))
    last={x:float(end[x]) for x in ['a','t','fplus_04','fminus_04','fplus_05','fminus_05','fplus_06','fminus_06',
            'genuine_plus','genuine_minus','rho','p','area_core_com_over_V','area_smooth_com_over_V',
            'phi_core_energy_proxy','phi_smooth_energy_proxy','v2_core_proxy','v2_smooth_proxy','q_fraction_core']}
    last['w_all_fields']=last['p']/last['rho']
    initial=None
    N=int(cfg['arguments']['N']);dx=float(cfg['arguments']['dx']);L=N*dx
    if name.startswith('gaussian'):
        initial=cfg['noise']
    elif cfg.get('mode_sha256'):
        a=np.loadtxt(RAW/f'modes_L{L:.8g}_s{cfg["seed"]}.txt')
        c=a[:,3:].reshape(-1,3,2)
        cov=.5*np.einsum('mfc,mgc->fg',c,c)
        kval=np.linalg.norm(a[:,:3]*2*np.pi/L,axis=1)
        kap=float(cfg['arguments']['kappa'])
        initial={'sparse_modes':len(a),'minimum_nonzero_k':float(kval.min()),'covariance':cov.tolist(),
                 'leading_mean_chiS_acceleration_cross_term':float(-2*kap**(2/3)*.01*.001*cov[2,0]),
                 'cross_term_scope':'Leading restored-branch perturbative term only; not the full force or a demonstrated causal attribution.'}
    return {'name':name,'N':N,'dx':dx,'dt':float(cfg['arguments']['dt']),'L':L,'seed':cfg['seed'],
            'complete':abs(float(end['t'])-48)<1e-8,'formation':form,'final':last,'initial_noise':initial,
            'Q_initial':q0,'charge_relative_drift':qdrift,'conformal_budget_residual':budget,
            'numerical_gates_pass':qdrift<=1e-9 and budget<=.01,
            'topology_final':top[-1],
            'every_postformation_topology_has_interface':bool(after) and all(t['interface_links']>0 for t in after),
            'every_postformation_sample_both_vacua_wind_all_axes':bool(after) and all(t['plus']['any_winding']==7 and t['minus']['any_winding']==7 for t in after),
            'postformation_topology_samples':len(after),'wall_seconds':cfg['wall_seconds']}

def relative(a,b):return abs(a-b)/max(abs(b),1e-300)

def comparisons(runs):
    out={}
    for control in ['dt_half_s101','grid256_s101','grid384_s101']:
        if control not in runs:continue
        ref=runs['main_s101'];other=runs[control]
        a=ref['final'];b=other['final']
        out[control]={'reference':'main_s101','endpoint_absolute_fraction_differences':{k:abs(a[k]-b[k]) for k in ['fplus_05','fminus_05','genuine_plus','genuine_minus']},
                      'endpoint_relative_difference':{k:relative(a[k],b[k]) for k in ['rho','area_core_com_over_V','area_smooth_com_over_V']},
                      'basin_area_relative_difference':relative(ref['topology_final']['basin_link_area_comoving_over_V'],other['topology_final']['basin_link_area_comoving_over_V']),
                      'same_raw_formation_presence':(ref['formation'] is None)==(other['formation'] is None)}
    if 'grid256_s101' in runs and 'grid384_s101' in runs:
        a=runs['main_s101'];b=runs['grid256_s101'];c=runs['grid384_s101']
        keys=['rho','area_core_com_over_V','area_smooth_com_over_V']
        e1={k:relative(a['final'][k],b['final'][k]) for k in keys}
        e2={k:relative(b['final'][k],c['final'][k]) for k in keys}
        # None of these scalar/area tests substitutes for formation convergence.
        out['spatial_refinement']={'successive_relative_errors_192_256':e1,'successive_relative_errors_256_384':e2,
              'all_scalar_area_errors_below_10pct_and_decreasing':all(e2[k]<.1 and e2[k]<e1[k] for k in keys),
              'formation_convergence_certified':False if all(r['formation'] is None for r in [a,b,c]) else None,
              'scope':'A nonforming seed cannot certify convergence of two-phase formation or network persistence.'}
    return out

def main():
    status=json.loads((HERE/'RUN_STATUS.json').read_text())
    if (HERE/'GAUSSIAN_STATUS.json').exists():status+=json.loads((HERE/'GAUSSIAN_STATUS.json').read_text())
    runs={s['name']:summary(s['name']) for s in status if s['status']=='completed'}
    result={'reviewed_claude_commit':'c7b8480a49e3f0b8737e800396a30f33f6337981',
            'protocol_freeze':'188b92b9375df468802cb064c8a288bcc38f956f',
            'validated_source_freeze':'998ddd8b1252e3ecd7cab8be4dc99044574f5faa',
            'observational_scope':'No new data comparison; no fitted mass, yield or abundance prediction.',
            'status':status,'runs':runs,'comparisons':comparisons(runs),
            'all_completed_numerical_gates_pass':all(x['numerical_gates_pass'] for x in runs.values()),
            'negative_controls_and_limits':'Retain one-sign outcomes, threshold-dependent cores and every failed validation. No finite-time result establishes cosmic persistence or a bias law.'}
    (HERE/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    cols=['name','N','L','a_form','end_fplus','end_fminus','end_plus_winding','end_minus_winding','budget']
    lines=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
    for name,r in runs.items():
        vals=[name,str(r['N']),f'{r["L"]:.1f}',f'{r["formation"]["a_first_sample"]:.2f}' if r['formation'] else 'none',
              f'{r["final"]["fplus_05"]:.6g}',f'{r["final"]["fminus_05"]:.6g}',str(r['topology_final']['plus']['any_winding']),
              str(r['topology_final']['minus']['any_winding']),f'{r["conformal_budget_residual"]:.3g}']
        lines.append('| '+' | '.join(vals)+' |')
    notes='\nRaw vacuum fractions use S>0.5 / S<-0.5. Topology also requires F^2<0.5. Winding bit masks: 1=x, 2=y, 4=z, 7=all three, 0=none. The budget is a fraction (0.00118 is 0.118 percent), evaluated at saved scalar times. Formation requires each raw vacuum fraction >=10 percent. None means that criterion was never reached.\n'
    (HERE/'RUN_TABLE.md').write_text('\n'.join(lines)+'\n'+notes)
    print('\n'.join(lines))
    print(json.dumps({'completed':len(runs),'all_numerical_gates':result['all_completed_numerical_gates_pass'],'comparisons':result['comparisons']},indent=2))

if __name__=='__main__':main()
