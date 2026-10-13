"""Validate published records and optionally rerun the small independent checks.

Does not silently rerun larger invariant evolutions. Works with python -I.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
import math
import os
import subprocess
import sys

import numpy as np

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))


def read(name):return json.loads((P/name).read_text())


def run(execute=False):
    executed=[]
    if execute:
        for name in ['exact_destination.py','coherent_ring.py','entropy_control.py']:
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
            subprocess.run([sys.executable,'-I',str(P/name)],check=True,env=env)
            executed.append(name)
    groups=[]
    assert hashlib.sha256((P/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
    groups.append('byte_exact_pinned_frozen_model')
    exact=read('exact_results.json')
    assert exact['status']=='CP12_COHERENT_DESTINATION_EXACT_PASS'
    for case in exact['first_step_controls']:
        M=case['N']*(case['N']-1)//2
        expected=F(case['effective_removal_norm_square'])/(100*M*M)
        assert all(F(r['coefficient'])==expected for r in case['rows'])
    groups.append('coherent_charged_and_complex_phase_common_vacancy_coefficients')
    difference=exact['primary_distance_two_minus_distance_three']
    assert all(F(difference[str(k)])==0 for k in [2,4,6,8])
    assert F(difference['10'])==-F(1103,3240000000)
    assert exact['primary_coherent_taylor']['no_hopping_destination_equality_through_order']==10
    rows={tuple(r['edge']):r['coefficients'] for r in exact['primary_dephased_taylor']}
    assert F(rows[(0,4)]['10'])-F(rows[(1,4)]['10'])==-F(1103,6480000000)
    groups.append('unregulated_higher_order_destination_and_record_controls')
    from coherent_ring import case
    ring=read('ring_results.json')
    for row in ring['cases']:
        assert row==case(row['N'],row['R'],bool(row['locked_connections']))
    groups.append('fresh_exact_coherent_ring_and_locking_counts')
    for L in range(2,7):
        data=read(f'finite_L{L}.json')
        assert data['N']==5 and data['Q']==4 and data['status']=='CP12_COHERENT_DESTINATION_FINITE_PASS'
        assert abs(data['initial_energies'][-1]-2.)<1e-11
        for row in data['rows']:
            assert row['norm_error']<2e-9 and row['energy_error']<2e-8
            assert np.ptp(row['neutral_only_occupations'])<2e-12
            u=row['unordered_pregeometry_control']
            assert abs(u['all_labelled_edges_mean_occupation']-.4)<2e-12 and u['edge_mean_spread']<2e-12
            for r in row['record_visibility_controls']:
                eta=r['record_overlap']
                expected=eta*np.array(row['coherent_occupations'])+(1-eta)*np.array(row['dephased_occupations'])
                assert np.max(abs(expected-r['occupations']))<1e-13
    a=read('finite_L5.json');b=read('finite_L6.json')
    occupation_error=max(np.max(abs(np.array(x['coherent_occupations'])-y['coherent_occupations'])) for x,y in zip(a['rows'],b['rows']))
    assert occupation_error<2e-14
    assert b['rows'][4]['near_minus_far_coherent']<0 and b['rows'][5]['near_minus_far_coherent']>0
    groups.append('complete_five_event_evolution_cutoffs_records_and_sign_reversal')
    inter=read('invariant_intertwiner_results.json')
    for c in inter['controls']:
        assert c['isometry_residual']<1e-10 and c['Hamiltonian_intertwiner_residual']<1e-9
        assert max(r['state_error'] for r in c['dynamics_controls'])<1e-9
        assert max(r['shape_error'] for r in c['dynamics_controls'])<1e-9
    groups.append('full_builder_signed_invariant_intertwiners_and_dynamics')
    shapes=read('unordered_exact_results.json')['controls']+read('unordered_exact_N7.json')['controls']
    for row in shapes:
        n=row['N'];Q=n-1;M=n*(n-1)//2
        mean=F(n*math.comb(n-1,2)*Q*(Q-1),M*(M-1))
        assert mean==F(row['uniform_mean_degree_wedges'])
        factor=F(row['degree_bias_factor']);ell=F(M-1,10*M)
        assert factor==ell*(F(17,4)-ell)/90>0
        assert row['lower_probability_coefficients_vanish_through_order']==4
        assert sum(F(c['probability_T6_coefficient']) for c in row['graph_class_coefficients'])==0
        for c in row['graph_class_coefficients']:
            assert F(c['relative_per_graph_T6_coefficient'])==factor*(mean-c['degree_wedges'])
    groups.append('exact_unordered_degree_bias_for_every_graph_class_N4_to_N7')
    large=[]
    for n,Ls in [(6,range(2,7)),(7,range(2,6))]:
        for L in Ls:
            data=read(f'invariant_N{n}_L{L}.json')
            assert data['status']=='CP12_UNORDERED_INVARIANT_EVOLUTION_PASS'
            assert data['Q']==n-1
            assert abs(data['initial_energy']-(n-1)/2)<1e-10
            assert abs(data['initial_variance']-2*(n-1))<1e-9
            assert all(r['norm_error']<2e-9 and r['energy_error']<2e-8 for r in data['rows'])
            assert data['rows'][-1]['connected']>data['initial_shapes']['connected']
        aa=read(f'invariant_N{n}_L{max(Ls)-1}.json');bb=read(f'invariant_N{n}_L{max(Ls)}.json')
        keys=['connected','spanning_path','degree_wedges','maximum_degree','isolated_vertices']
        error=max(abs(x[k]-y[k]) for x,y in zip(aa['rows'],bb['rows']) for k in keys)
        assert error<2e-8
        large.append(dict(N=n,cutoff_pair=[max(Ls)-1,max(Ls)],maximum_shape_difference=error,
                          T16_connected=bb['rows'][-1]['connected'],T16_path=bb['rows'][-1]['spanning_path']))
    groups.append('complete_invariant_formation_N6_N7_and_observed_cutoff_stability')
    entropy=read('entropy_results.json')
    assert entropy['status']=='CP12_DIRECT_ENTROPY_AND_UNITARY_CONTROL_PASS'
    for c in entropy['finite_controls']:
        assert c['B_norm']<=c['B_norm_upper']+1e-10
        assert 0<=c['thermal_relative_entropy']+1e-12<=2*c['B_norm_upper']+1e-10
        for row in c['rows']:
            assert row['mixed_relative_entropy']<=c['exact_all_time_mixed_entropy_upper']+1e-10
            assert row['mixed_data_processing_gap']>-1e-10
        assert abs(c['rows'][0]['pure_relative_entropy']-c['rows'][0]['mixed_relative_entropy']-math.log(c['graph_count']))<1e-10
    counter=entropy['noncommuting_adversarial_control']
    assert counter['actual_probability']>1000*counter['invalid_multiplicative_bound']
    groups.append('finite_entropy_data_processing_unitary_rank_and_noncommuting_adversary')
    for c in entropy['conditional_gauge_sector_controls']:
        g=.1*c['Q']*(c['N']*(c['N']-1)/2-c['Q'])/(c['N']*(c['N']-1)/2)
        assert c['neutral_diagonal_max']<=c['Q']/2+g+1e-12
        assert c['neutral_diagonal_min']>=c['Q']/2-1e-12
    groups.append('separately_defined_invariant_Gibbs_ensemble_controls')
    out=dict(status='CP12_COHERENT_FORMATION_AND_ENTROPY_CHECKS_PASS',groups=groups,count=len(groups),
             small_scripts_freshly_executed=executed,five_event_cutoff_max_occupation_difference=occupation_error,
             large_invariant_controls=large,
             scope='Recorded larger evolutions plus fresh exact ring checks and optional small scripts; not a rerun of every large calculation or a certified continuum limit.')
    (P/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--execute-small',action='store_true')
    a=ap.parse_args();run(a.execute_small)
