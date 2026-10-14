"""Fresh small certificates plus explicit audits of saved large eigensolves."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import sys
import tempfile
from datetime import datetime,timezone
from fractions import Fraction as F

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
import numpy as np
import scipy
import representation_ground as rep
from cutoff_certificate import certify
from analyze import uniform_counts
from ground import topology
import frozen_model as frozen

def read(name):return json.loads((P/name).read_text())
def close(a,b,tol=1e-9):assert abs(a-b)<tol,(a,b,abs(a-b))
def main():
    groups=[]
    def passed(name,kind,details):
        groups.append(dict(name=name,kind=kind,details=details))
        print(name,': PASS',flush=True)
    assert hashlib.sha256((P/'frozen_model.py').read_bytes()).hexdigest()==\
        'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
    assert (frozen.MASS,frozen.HOP,frozen.COST,frozen.ELECTRIC,frozen.GAMMA)==(1.,1.,.5,.25,.1)
    passed('frozen_law_fidelity','source hash and declared constants','Frozen implementation hash matches.')
    small=[]
    with tempfile.TemporaryDirectory(prefix='cp12-verify-') as td:
        previous=rep.P
        try:
            rep.P=Path(td)
            for n,L in [(3,2),(4,1)]:
                rep.full_check(n,L)
                small.append(json.loads((Path(td)/f'representation_check_N{n}_L{L}.json').read_text()))
        finally:rep.P=previous
    passed('complete_small_spectra','fresh full-matrix eigensolves',[
        dict(N=r['N'],dimension=r['full_dimension'],max_error=r['max_full_spectrum_difference']) for r in small])
    certs=[]
    for saved in read('cutoff_certificates.json')['certificates']:
        fresh=certify(saved['N'],saved['L'],F(saved['lower_rational']))
        assert fresh['all_rational_pivots_positive']
        close(fresh['upper'],saved['upper'],1e-12)
        certs.append(fresh)
    passed('unbounded_rotor_certificates','fresh exact rational LDL and trial energies',certs)
    sector_summary=[]
    for n in [5,6]:
        r=read(f'representations_N{n}_L3.json')
        assert r['algorithm']=='stabilizer-multiplicity-v2'
        assert r['status']=='CP12_ALL_REPRESENTATIONS_GROUND_PASS'
        shapes={tuple(x['shape']) for x in r['irreps']}
        assert shapes==set(rep.partitions(n))
        accounted=sum(x['representation_dimension']*x['multiplicity_dimension'] for x in r['irreps'])
        assert accounted==r['full_dimension']==r['dimension_accounted_for']
        assert max(x['residual'] for x in r['irreps'])<2e-8
        minimum=min(r['irreps'],key=lambda x:x['lowest_energy'])
        assert minimum['shape']==[n]
        inv=read(f'ground_N{n}_L3_invariant.json')
        close(minimum['lowest_energy'],inv['energies'][0])
        close(minimum['observables']['connected'],inv['observables']['connected'])
        gap=min([x['lowest_energy'] for x in r['irreps'] if x['shape']!=[n]]
                +[inv['energies'][1]])-minimum['lowest_energy']
        if n==5:
            full=read('ground_N5_L3_full.json')
            close(full['energies'][0],minimum['lowest_energy'])
            close(full['gap'],gap)
            close(full['observables']['connected'],minimum['observables']['connected'])
            for key in ['connected_energy','disconnected_energy','cross_norm']:
                close(full['sector_comparison'][key],inv['sector_comparison'][key])
        sector_summary.append(dict(N=n,sectors=len(shapes),full_dimension=accounted,global_gap=gap))
    for n in [4,5]:
        full=read(f'ground_N{n}_L3_full.json');inv=read(f'ground_N{n}_L3_invariant.json')
        close(full['energies'][0],inv['energies'][0])
        assert max(full['eigen_residuals'])<1e-8
    passed('all_sector_large_ground','saved numerical records and independent-builder cross-agreement',sector_summary)
    refinements=[]
    for n in [4,5,6]:
        a=read(f'ground_N{n}_L3_invariant.json');b=read(f'ground_N{n}_L4_invariant.json')
        de=abs(a['energies'][0]-b['energies'][0])
        dp=abs(a['observables']['connected']-b['observables']['connected'])
        assert de<1e-10 and dp<1e-10
        assert b['cutoff_boundary_probability']<a['cutoff_boundary_probability']
        refinements.append(dict(N=n,energy_difference=de,connected_difference=dp))
    for n in [5,6]:
        c=read(f'cutoff_comparison_N{n}_L3.json')
        assert c['status']=='NUMERICAL_ALL_SECTOR_COMPARISON_PASS'
        assert {tuple(x['shape']) for x in c['sectors']}==set(rep.partitions(n))
        assert all(x['comparison_minimum']>100*x['residual'] for x in c['sectors'])
        assert not c['exact_rational_certificate']
    assert read('ground_N7_L2_invariant.json')['cutoff_boundary_probability']>1e-6
    passed('regulator_checks','saved refinements and all-sector floating comparison minima',refinements)
    analysis=read('analysis.json')
    for row in analysis['rows']:
        n=row['N'];assert uniform_counts(n)==row['uniform_labelled_graph_count']
        close(sum(row['connected_cycle_length_probability'].values()),1.)
        close(sum(row['largest_component_size_probability'].values()),1.)
        close(sum(int(k)*v for k,v in row['largest_component_size_probability'].items())/n,
              row['largest_component_fraction'])
        raw=read(f"ground_N{n}_L{row['L']}_invariant.json")
        es=list(__import__('itertools').combinations(range(n),2))
        ring=sum(r['probability'] for r in raw['graph_class_probability']
                 if topology(n,r['graph_mask'],es)[0] and topology(n,r['graph_mask'],es)[5])
        close(ring,row['single_ring_probability'])
        if n==7:assert 'provisional' in row['scope']
    passed('graph_counts_and_shapes','fresh uniform enumeration and saved-state projector accounting',
           'N=4..7 counts; N=7 remains explicitly provisional.')
    for row in analysis['rows']:
        if 'block_comparison' not in row:continue
        s=row['block_comparison'];p=row['connected_probability'];E=row['energy']
        ec=s['connected_energy'];ed=s['disconnected_energy'];b=s['cross_norm']
        dc=ec-E;dd=ed-E
        lower=dd**2/(dd**2+b*b);upper=b*b/(dc**2+b*b)
        assert lower-1e-9<=p<=upper+1e-9
        close(lower,s['connected_probability_lower']);close(upper,s['connected_probability_upper'])
        close(E,p*s['conditional_connected_energy']+(1-p)*s['conditional_disconnected_energy']+
              s['graph_mixing_energy'])
        assert s['connectedness_dephasing_energy_cost']>0
    passed('projector_and_coherence_bounds','recomputed algebra from saved eigensolve certificates',
           'Two projected eigen-equations and complete energy decomposition.')
    result=dict(status='CP12_GROUND_CONNECTIVITY_VERIFIED',groups=groups,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        large_eigensolves_rerun_by_this_verifier=False,
        scientific_limit='Finite Q=N ground diagnostics; no thermodynamic geometry or gravity claim.')
    (P/'VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],len(groups),'groups',flush=True)

if __name__=='__main__':main()
