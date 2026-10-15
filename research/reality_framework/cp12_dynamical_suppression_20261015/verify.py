"""Verify the completed packet and regenerate its content manifest."""
from fractions import Fraction
from pathlib import Path
import hashlib,json,subprocess

P=Path(__file__).resolve().parent;RF=P.parent;ROOT=RF.parent.parent
BASE='768c46c96eae404a527d9cfee8f1b09abeb61a13'
PREFIX=P.relative_to(ROOT).as_posix()+'/'


def load(name):return json.loads((P/name).read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])


def main():
    frozen={
        'frozen_model.py':'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a',
        'exact_destination.py':'1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b'}
    for name,sha in frozen.items():assert digest(RF/'cp12_ground_connectivity'/name)==sha
    deps=[
        'cp12_ground_connectivity/frozen_model.py','cp12_ground_connectivity/exact_destination.py',
        'cp12_ground_connectivity/invariant_evolution.py',
        'charged_density_dynamics_20261015/fast_symmetry.py','charged_density_dynamics_20261015/symmetry.cpp',
        'charged_density_dynamics_20261015/stream.py',
        'cp12_shortcut_dynamics_20261015/operators.py','cp12_shortcut_dynamics_20261015/integrals.py',
        'cp12_shortcut_dynamics_20261015/recorded_transport.py']
    dep_hashes={}
    for name in deps:
        path=RF/name;rel=path.relative_to(ROOT).as_posix()
        before=subprocess.check_output(['git','show',f'{BASE}:{rel}'],cwd=ROOT)
        assert path.read_bytes()==before
        dep_hashes[rel]=digest(path)
    changed=subprocess.check_output(['git','diff','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    assert all(path.startswith(PREFIX) for path in changed)
    five=load('N5_B50.json');six=load('N6_SHELL2.json');checks=load('CHECKS.json')
    boundary=load('BOUNDARY_CHECKS.json');rare=load('N5_RARE_PAULI.json')
    pauli=load('PAULI_ALGEBRA.json');response=load('N4_CORE_RESPONSE.json')
    exterior=load('NEUTRAL_EXTERIOR.json')
    for d in [five,six,checks,boundary,rare,pauli,response,exterior]:assert d['status']=='PASS'
    for r in exterior['rows']:
        assert abs(r['neutral_exterior_action_square']-r['endpoint_vacancy_formula'])<1e-11
        assert r['neutral_exterior_other_outside_cross_abs']<1e-11
    assert checks['physical_hopping_Gauss_checks']==24636
    assert checks['full_physical_dephased_reference_error']<1e-12
    assert checks['orbit_diagonal_reference_would_be_wrong_by']>1
    assert checks['adaptive_rewiring_exit_transitions']>0
    assert checks['matrix_free_action_and_residual_difference']<1e-11
    assert checks['matrix_free_inherited_propagator_difference']<1e-11
    for c in checks['integer_hopping_bound_controls'].values():
        assert all(k>0 for k in c['positive_leading_minors'])
    for c in checks['exact_spectral_integral_controls']:
        assert c['quadrature_interval'][0]<=c['exact_spectral_integral']<=c['quadrature_interval'][1]
    assert boundary['maximum_full_column_entry_error']<1e-10
    assert boundary['sampled_columns']>=40
    assert boundary['dimension']==six['dimension']==2477611
    assert boundary['outgoing_dimension']==six['outside_dimension']==1658455
    eta_checks=0
    for data in [five,six]:
        assert data['hopping_norm_bound']==(14 if data['N']==5 else 17)
        for row in data['rows']:
            p=row['formation_probability'];pI=row['formation_interval'];assert pI[0]<=p<=pI[1]
            assert pI[0]>0
            assert row['record_entropy_interval_bits'][0]<=row['record_entropy_bits']<=row['record_entropy_interval_bits'][1]
            for r in row['probes']:
                ints=r['joint_integrals'];out=ints['outside']['full_law_interval'];inside=ints['inside']['full_law_interval']
                for z in ints.values():
                    assert z['retained_quadrature_interval'][0]<=z['value']<=z['retained_quadrature_interval'][1]
                    assert z['full_law_interval'][0]<=z['full_law_interval'][1]
                if inside[0]>0:
                    assert abs(r['eta_interval'][0]-out[0]/inside[1])<1e-12
                    assert abs(r['eta_interval'][1]-out[1]/inside[0])<1e-12
                else:assert r['eta_interval'][1] is None
                for key in ['outside','inside']:
                    expected=[ints[key]['full_law_interval'][0]/pI[1],
                        ints[key]['full_law_interval'][1]/pI[0]]
                    assert max(abs(a-b) for a,b in zip(expected,r['conditional_integrals'][key]))<1e-10
                eta_checks+=1
    assert all(r['eta_interval'][0]>2.53 for row in five['rows'] for r in row['probes'])
    r5=next(r for r in five['rows'] if r['wait']==2.)
    p5=next(r for r in r5['probes'] if r['lag']==.25)
    assert p5['interference_outside_interval'][0]>2.4
    assert p5['neutral_edge_vacancy_fraction_interval'][0]>.52
    old=json.loads((RF/'cp12_shortcut_dynamics_20261015/SPLIT_N6_B12.json').read_text())
    adaptive=load('N6_ADAPT225K_PILOT.json')
    improvement={};matched={}
    for row in six['rows']:
        wait=row['wait'];before=next(r for r in old['trajectory'] if r['time']==wait)
        middle=next(r for r in adaptive['rows'] if r['wait']==wait)
        assert row['waiting_state_error']<middle['waiting_state_error']<before['full_state_error']
        assert overlap(before['event_interval'],row['formation_interval'])
        assert overlap(middle['formation_interval'],row['formation_interval'])
        improvement[str(wait)]=dict(before=before['full_state_error'],after=row['waiting_state_error'],
            factor=before['full_state_error']/row['waiting_state_error'])
        five_row=next(r for r in five['rows'] if r['wait']==wait)
        a=next(r for r in five_row['probes'] if r['lag']==.25)['eta_interval']
        b=next(r for r in row['probes'] if r['lag']==.25)['eta_interval']
        conclusion=('resolved increase' if b[0]>a[1] else
            'resolved decrease' if b[1] is not None and b[1]<a[0] else 'unresolved comparison')
        matched[str(wait)]=dict(N5=a,N6=b,comparison=conclusion)
    for row in pauli['families']:
        n=row['N'];assert Fraction(row['filled_initial_inside_action_square'])==2*n
        assert Fraction(row['filled_outside_H_action_square'])==4*n
        assert Fraction(row['filled_eta_T2_coefficient'])==Fraction(2,3)
        assert Fraction(row['one_vacancy_eta_initial_limit'])==Fraction(1,n-1)
        assert row['direct_long_ring_distance_transfer_probability_T4_coefficient']=='1/4'
    assert rare['formation_interval'][0]>0 and rare['formation_interval'][1]<1e-6
    assert rare['initial_outside_action_square']==0
    assert abs(rare['eta_T2_coefficient']-2/3)<1e-12
    assert all(r['phase_pulse_fixed_receiver_response']==0 for r in rare['rows'])
    assert next(r for r in rare['rows'] if r['lag']==.25)['eta_interval'][1]<.04
    neighbor=next(r for r in response['rows'] if r['group']=='all' and r['lag']==.25)
    assert neighbor['responses']['core_neighbors_distance1']['zero_excluded']
    assert neighbor['responses']['core_neighbors_distance1_and_same_ring']['zero_excluded']
    for name in ['N6_ADAPT225K_BUILD.json','N6_SHELL1_BUILD.json','N6_SHELL2_BUILD.json']:
        e=load(name)['extension'];assert e['base_cross_H_error']<1e-10 and e['base_cross_V_error']<1e-10
        assert e['H_hermiticity']<1e-8 and e['V_hermiticity']<1e-8
    files=[p for p in sorted(P.iterdir()) if p.is_file() and p.name not in ['manifest.json','VERIFICATION.json']]
    assert all(p.suffix in ['.py','.json','.md','.svg'] for p in files)
    for name in ['REPORT.md','HANDOFF.md','README.md','DERIVATION.md']:
        s=(P/name).read_text();assert 'TO_BE_INSERTED' not in s and len(s)>1000
    assert (P/'dynamical_suppression.svg').stat().st_size>1000
    manifest=dict(base_commit=BASE,base_tree='f180c6ad1e3b1040aa899dd704a9ed7b309b5efd',
        packet=PREFIX.rstrip('/'),files={p.name:dict(sha256=digest(p),bytes=p.stat().st_size) for p in files},
        inherited_dependency_sha256=dep_hashes,
        exclusions='manifest.json and VERIFICATION.json excluded to avoid self-reference; large reproducible runtime caches excluded.',
        law='All inherited dependencies are byte-identical to the published parent; new numerical evolutions use full frozen H.')
    (P/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    verification=dict(status='PASS',base_commit=BASE,packet_files_hashed=len(files),
        frozen_sources_unchanged=frozen,inherited_dependencies_unchanged=len(deps),
        action_interval_consistency_checks=eta_checks,full_physical_hops_checked=24636,
        independent_large_boundary_columns=boundary['sampled_columns'],
        six_site_error_improvement=improvement,matched_early_time_eta_comparison=matched,
        rare_event_is_not_appreciable=True,rare_immediate_density_phase_response_exactly_zero=True,
        auxiliary_fixed_ring_neighbor_response_resolved=True,
        acceptance='The three-part physical locality gate is not passed by the tested experiments. Finite-size diagnostics do not decide every CP-12 phase.',
        manifest_sha256=digest(P/'manifest.json'),
        scope='Analytic truncation/quadrature bounds with ordinary floating-point reporting allowance, exact rational controls and exact integer norm certificates. No outward-rounded formal certificate.')
    (P/'VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification,indent=2))


if __name__=='__main__':main()
