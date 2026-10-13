"""Validate recorded independent controls and reproduce small exact certificates.

This is deliberately not a second large-H eigenvalue run. Individual
reproduction commands are in REPORT.md and MANIFEST.json.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal
import hashlib
import itertools
import json
import math

import numpy as np

P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())


def run():
    groups=[]
    sources=read('parent_sources.json')
    for source in sources['files']:
        data=(P.parents[2]/source['path']).read_bytes()
        assert hashlib.sha256(data).hexdigest()==source['sha256']
    groups.append('pinned_frozen_and_parent_inputs')
    exact=read('exact_short_time_results.json')
    assert exact['status']=='CP12_EXACT_SHORT_TIME_CHECKS_PASS'
    for case in exact['cases']:
        n=case['N'];q=case['Q'];M=n*(n-1)//2;d=case['source_degree'];total=Fraction(0)
        for edge in case['edges']:
            actual=Fraction(edge['exact_cubic_times_sqrt_degree']);total+=actual
            expected=Fraction(-(M-q),100*M*M) if edge['initial_present'] and 0 in edge['edge'] else (
                     Fraction(d,100*M*M) if not edge['initial_present'] else Fraction(0))
            assert actual==expected==Fraction(edge['derived_formula_times_sqrt_degree'])
        assert total==0
    groups.append('unregulated_exact_cubic_response_and_Q_sum_rule')
    for filename in ['local_response_path_L3.json','local_response_cycle_L3.json',
                     'local_response_star_L2.json','local_response_star_L3.json']:
        data=read(filename)
        for row in data['rows']:
            assert row['occupation_change_probability']<=row['probability_bound']+2e-10
            assert abs(row['response'])<=row['response_bound']+2e-10
            assert row['mean_outside_bonds']<=row['outside_mean_bound']+2e-10
            assert abs(row['fixed_block_response'])<2e-10
            if 'impulse_error' in row:assert row['impulse_error']<2e-7
        assert max(row['cubic_error'] for row in data['taylor'])<3e-11
    groups.append('complete_finite_local_connection_response_and_independent_pulses')
    neutral=read('neutral_spectrum_results.json')
    for case in neutral['independent_sparse_controls']:
        assert max(row['amplitude_error'] for row in case['checks'])<2e-11
    for row in neutral['complete_graph_survival_bounds']:
        n=row['N'];M=n*(n-1)//2;q=row['Q'];alpha=Fraction(q,M-q+1)
        assert Fraction(row['exact_all_time_departure_upper_bound'])==4*alpha*(1-alpha)
    groups.append('all_time_neutral_survival_and_flat_vacancy_spectrum')
    top=read('top_harmonic_dressing_results.json')
    for row in top['nested_fiber_controls']:
        assert abs(row['numerical_unlock_leakage_norm_square']-float(Fraction(row['exact_unlock_leakage_norm_square'])))<1e-11
        assert row['projector_identity_residual']<1e-11 and row['top_preserving_locking_residual']<1e-11
    full=top['full_control']; assert full['exact_V_top_identity_residual']<1e-11
    record=full['projected_preparation_record_control']
    assert abs(record['system_energy_change']-record['exact_change'])<1e-11 and record['system_energy_change']<0
    for row in full['projected_geometry_commutator_controls']:
        assert 0<row['projected_commutator_norm']<=row['exact_complement_product_bound']+1e-11
    from harmonic_ground_certificate import certificate
    for row in read('harmonic_ground_certificate_results.json')['rows']:
        assert row==certificate(row['N'])
        assert Decimal(row['top_harmonic_ground_probability_upper'])<=1
    groups.append('locking_unlocking_record_algebra_and_full_law_ground_certificates')
    tensor=read('tensor_channel_results.json')
    assert tensor['metric']['linear_connection_metric_rank']==4 and tensor['metric']['pair_tensor_rank']==6
    star=tensor['physical_star'];assert star['dimension']==126 and star['complete_star_basis_at_L1']
    assert star['total_E_spectral_weight']>1e-8
    for row in star['simultaneous_permutation_controls']:
        assert abs(row['E_character']-row['expected_E_character'])<1e-9
        assert abs(row['T_character']-row['expected_T_character'])<1e-9
    for row in star['responses']:assert row['impulse_error']<2e-8
    assert abs(star['identical_boundary_control']['difference']-2.065)<1e-12
    groups.append('independent_tensor_components_exact_star_response_and_energy_boundary_control')
    full2=read('full_tensor_L2.json');full3=read('full_tensor_L3.json')
    assert abs(full2['ground_energy']-full3['ground_energy'])<1e-9
    assert abs(full2['full_gap']-full3['full_gap'])<1e-9
    for data in [full2,full3]:
        assert max(data['eigen_residuals'])<1e-7
        for row in data['lowest_shell_permutation_representation']:
            assert abs(row['lowest_shell_character']-row['unordered_pair_minus_vertex_character'])<1e-7
        for channel in data['channels']:
            assert max(channel['linear_solve_residuals'])<1e-8
            assert channel['first_active_shell']['degeneracy']==5
            assert min(np.linalg.eigvalsh(channel['static_susceptibility']))>=-1e-9
    groups.append('complete_fluctuating_graph_tensor_spectra_cutoffs_and_static_matching')
    intertwiner=read('quotient_intertwiner_results.json')
    for row in intertwiner['controls']:
        assert row['isometry_residual']<1e-10 and row['full_H_intertwining_residual']<1e-9
        assert row['exact_integer_checks']['exact_H_C_commutator_nonzeros']==0
    groups.append('independent_full_builder_signed_quotient_intertwiners')
    quotient_files=sorted(P.glob('marked_pair_N*.json'))
    assert (P/'marked_pair_N7_Q6_L2.json') in quotient_files
    for filename in quotient_files:
        data=json.loads(filename.read_text());n=data['N']
        assert data['full_pair_irrep_dimension']==n*(n-3)//2
        for row in data['low_spectrum']:
            assert row['casimir_residual']<1e-7 and row['eigen_residual']<1e-7
        for row in data['invariant_activity_controls']:assert row['impulse_error']<2e-8
        if 'exact_integer_symmetry_controls' in data:
            assert data['exact_integer_symmetry_controls']['exact_C_cubic_polynomial_nonzeros']==0
    n6L2=read('marked_pair_N6_Q5_L2.json');n6L3=read('marked_pair_N6_Q5_L3.json')
    assert abs(n6L2['invariant_ground_energy']-n6L3['invariant_ground_energy'])<1e-9
    assert abs(n6L2['first_pair_shell']['gap']-n6L3['first_pair_shell']['gap'])<1e-9
    assert abs(read('marked_pair_N5_Q4_L2.json')['invariant_ground_energy']-full2['ground_energy'])<1e-10
    groups.append('larger_exact_symmetry_sectors_label_multiplicities_and_cutoff_controls')
    gr=read('spin2_consistency_results.json')
    for row in gr['canonical_constraint_controls']:
        assert row['constraint_rank']==4 and row['tt_projector_rank']==2 and row['physical_phase_dimension']==4
        assert row['first_class_bracket_norm']<1e-12
    for case in gr['universal_soft_coupling_controls']:
        for row in case['controls']:
            g1,g2=row['species_couplings']
            assert abs(row['gauge_polarization_change'])<1e-12 if g1==g2 else abs(row['gauge_polarization_change'])>1e-4
    groups.append('conventional_spin2_constraints_and_universal_soft_coupling_target')
    out=dict(status='CP12_GRAVITY_RESPONSE_CHECKS_PASS',groups=groups,count=len(groups),
             scope='Record validation plus recomputed outward certificates; individual independent simulation commands are recorded separately.')
    (P/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    return out


if __name__=='__main__':run()
