"""Read-only audit of S6 and the separate G1 ledger; no new trajectories."""
import hashlib,json
from collections import Counter
from math import sqrt
from pathlib import Path
import campaign as P
import conditioned as C

BASE=Path(__file__).resolve().parent
PARENT_MANIFEST_SHA='9c5d0b7178d37348462d42da3972bd46b95f3eca713b06ad303b2dcabe20a096'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return json.loads((BASE/name).read_text())


def case(r):
    assert r['N']==sum(r['nonroot_type_counts'])+1
    assert r['receiver_modes']==64*r['N']**2
    assert r['observation_time']==4*r['N']**2
    assert r['maximum_word_rank']<=10
    assert r['maximum_coefficient_dimension']<=2_000_000
    assert r['maximum_cached_overlap_entries']<=8_000_000
    assert r['exact_ordered_histories'] and r['all_reverse_transitions'] and r['conditional_frames']
    assert r['history_count']==len(C.S5.type_words(r['nonroot_type_counts']))
    assert not r['extra_photon_cutoff'] and not r['roundoff_certified']
    assert abs(sum(r['populations'])-1)<3e-9
    assert abs(sum(r['error_components'].values())-r['analytic_state_error_account'])<3e-12
    assert abs(r['dt']*r['steps']-r['observation_time'])<1e-10
    assert abs(r['completed_steps_time']-r['dt']*r['completed_steps'])<1e-10
    assert (r['status']=='COMPLETE')==(r['completed_steps']==r['steps'])
    e=r['error_components']
    gates=(e['strang']<=.15+1e-10 and e['orbital_projection']<=.20+1e-9
           and e['projected_emission']<=.20+1e-10
           and sum(v for k,v in e.items() if k not in ('strang','orbital_projection','projected_emission'))<=.05
           and r['analytic_state_error_account']<=.60)
    if r['status']=='COMPLETE':
        def lower(p):
            p=min(1,max(0,p)); error=r['analytic_state_error_account']
            return max(0,max(0,sqrt(p)-error)**2,1-(sqrt(1-p)+error)**2)
        assert abs(r['formation_lower']-lower(r['macro_probability_estimate']))<2e-12
        assert abs(r['full_assembly_lower']-lower(r['full_assembly_probability_estimate']))<2e-12
        assert r['gate_pass']==bool(gates and r['formation_lower']>=.5)
    else:
        assert r['formation_lower'] is None and r['full_assembly_lower'] is None
        assert not r['gate_pass'] and r['incomplete_population_notice']


def check_families(raw):
    assert len(raw['attempts'])==10 and len(raw['families'])==6
    for r in raw['attempts']:case(r)
    for family,expected in zip(raw['families'],P.families()):
        for k,v in expected.items():assert family[k]==v
        lower=full=missing=0
        for root,c in enumerate(family['components']):
            target=P.component(expected,root)
            for k,v in target.items():assert c[k]==v
            r=raw['attempts'][c['index']]
            assert P.key(r['nonroot_type_counts'],r['detuning_values'],r['g']) in (
                P.key(c['counts'],c['values'],c['g']),P.mirror(c['counts'],c['values'],c['g']))
            if r['status']=='COMPLETE':
                lower+=c['weight']*r['formation_lower'];full+=c['weight']*r['full_assembly_lower']
            else:missing+=c['weight']
        assert abs(sum(c['weight'] for c in family['components'])-1)<1e-12
        assert abs(lower-family['formation_lower'])<2e-12
        assert abs(full-family['full_assembly_lower'])<2e-12
        assert abs(missing-family['unresolved_root_weight'])<1e-12
        assert (family['status']=='COMPLETE')==(missing==0)
        assert family['gate_pass']==bool(missing==0 and lower>=.5 and
                  all(raw['attempts'][c['index']]['gate_pass'] for c in family['components']))


def audit():
    original,extension=read('RESULTS.json'),read('GLOBAL_RESULTS.json')
    assert extension['original_S6_results_sha256']==sha(BASE/'RESULTS.json')
    assert not extension['physical_law_changed'] and not extension['physical_test_changed']
    assert not extension['global_error_limits_changed']
    check_families(original);check_families(extension)
    parent=BASE.parent/'wm2_history_formation_20261016'
    assert sha(parent/'manifest.json')==PARENT_MANIFEST_SHA
    for name,item in json.loads((parent/'manifest.json').read_text())['files'].items():
        assert sha(parent/name)==item['sha256'],'Prior S5 file changed: '+name
    controls=original['checks']
    assert len(controls['fock_checks'])==5 and len(controls['cross_history_checks'])==3
    assert len(controls['emission_step_checks'])==3 and len(controls['dynamics_checks'])==6
    assert max(abs(x['interference_change']) for x in controls['cross_history_checks'])>1e-4
    assert all(x['state_error']<=x['residual_integral_upper']+x['polynomial_tail']+3e-12
               for x in controls['emission_step_checks'])
    assert all(x['full_wavefunction_error']<=x['analytic_state_error_account']+3e-10
               for x in controls['dynamics_checks'])
    physical=controls['physical_receiver'];uniform=controls['uniform_receiver']
    case(physical['trajectory']);case(uniform['trajectory'])
    assert physical['complete_physical_dimension']==333505
    assert physical['full_wavefunction_error']<=physical['trajectory']['analytic_state_error_account']+3e-9
    assert uniform['quadratic_full_state_error']<=uniform['trajectory']['analytic_state_error_account']+3e-9
    extra=read('ADDITIONAL_CONTROLS.json')
    assert all(x['bit_identical'] for x in extra['rotation_cache_identity'])
    assert max(x['error'] for x in extra['matrix_free_stress'])<5e-11
    assert max(x['common_frame_S5_operator_error'] for x in extra['common_frame_limit'])<3e-11
    action=read('ACTION_CONTROLS.json')
    assert not action['used_in_primary_trajectories'] and not action['new_formation_result']
    assert len(action['coherent_density_controls'])==4
    assert max(x['physical_density_error'] for x in action['coherent_density_controls'])<5e-11
    assert action['exact_cancellation_density_norm']<5e-11
    assert max(x['missing_after_one_direction'] for x in action['one_photon_action_capture'])<5e-11
    assert max(max(x['all_path_support_errors']) for x in action['exact_Krylov_enclosure'])<5e-11
    attempts=original['attempts']+extension['attempts']
    return dict(numerical_consistency='PASS',
       scientific_status='Tested cross-history transition/error calculus; all-size formation OPEN',
       S6_families_complete=sum(f['status']=='COMPLETE' for f in original['families']),
       S6_families_pass=sum(f['gate_pass'] for f in original['families']),
       S6_stop_counts=dict(Counter(r['stop_reason'] for r in original['attempts'])),
       G1_families_complete=sum(f['status']=='COMPLETE' for f in extension['families']),
       G1_families_pass=sum(f['gate_pass'] for f in extension['families']),
       G1_stop_counts=dict(Counter(r['stop_reason'] for r in extension['attempts'])),
       original_unique_components=10,G1_unique_components=10,
       maximum_coefficient_dimension=max(r['maximum_coefficient_dimension'] for r in attempts),
       maximum_word_rank=max(r['maximum_word_rank'] for r in attempts),
       maximum_cached_overlap_entries=max(r['maximum_cached_overlap_entries'] for r in attempts),
       physical_reference_modes=576,physical_reference_dimension=333505,
       physical_reference_time=physical['reference_time'],
       physical_reference_full_state_error=physical['full_wavefunction_error'],
       physical_reference_error_account=physical['trajectory']['analytic_state_error_account'],
       small_complete_wavefunctions=6,unequal_rank_Fock_checks=5,cross_history_checks=3,
       emission_step_residual_checks=3,common_frame_S5_checks=3,large_matrixfree_checks=3,
       coherent_action_density_checks=4,Krylov_enclosure_orders=[1,2,3],
       action_enrichment_used_for_primary_formation=False,
       S5_unchanged=True,original_S6_preserved=True,arithmetic='ordinary double precision',
       outward_rounded=False,external_independent_validation=False,
       uniform_active_rank_proved=False,all_size_independent_disorder_formation_proved=False,
       six_event_target_completed=any(r['N']==6 and r['status']=='COMPLETE' for r in extension['attempts']))


if __name__=='__main__':print(json.dumps(audit(),indent=2))
