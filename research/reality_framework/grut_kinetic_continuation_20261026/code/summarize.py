"""Package actual control results; leave the unanswered physical outputs null."""
import json,hashlib,subprocess
import numpy as np
from bootstrap import ROOT,REPO,BASE,COMMIT
from moment_controls import save
import ks_common as K

def main():
    cases=json.loads((ROOT/'MOMENT_CONTROLS.json').read_text())['cases']
    follow=json.loads((ROOT/'FOLLOWUP_CONTROLS.json').read_text())['cases']
    rows=[]
    for label,r in cases.items():
        rows.append(dict(label=label,R_total=r['R'],formation=r['formation'],**r['final'],stop=r['stop']))
    obj=K.THM['objects'];qbind=obj['E_over_Q_equals_m_at']['Q']
    save('RESULTS.json',dict(status='Conditional controls; full kinetic result not obtained',
        final_physical_f_sym=None,physical_kinetic_threshold=None,formation_kernel_derived=False,
        primary_denominator_GeV_per_visible_entropy=K.DM_PER_S,
        allowance=K.F_ALLOW,rows=rows,
        diagnostic_gas_scan=[dict(R_total=r['R'],f_at_stop=r['final']['f_at_stop'],stop=r['stop'])
            for label,r in follow.items() if label.startswith('gas_scan')],
        object_scale=dict(locally_stable_branch_floor=obj['Q_min']['Q'],E_over_Qm_at_floor=obj['Q_min']['E_over_QmX'],
            absolute_binding_crossing=qbind,fastest_fragment_charge_approx=25000.,
            minimum_equal_fragment_count_by_charge_only=int(np.ceil(qbind/25000.)),
            statement='Three equal primary fragments would exceed the binding crossing by charge addition only. No merger probability, retained charge, or stable product is inferred.'),
        scope='R is total extra conformal energy including mediator and vector terms. The gas diagnostic R=40.5626/48.5 is total R, not the full source-plus-mediator row. All starts are forced chemical/kinetic equilibrium.'))
    paths=['HANDOFF_KINETIC.md','BRIEF_KINETIC.md',
        'tracks/kin_source/results/source_interface.json',
        'tracks/kin_rates/results/rates_interface.json',
        'tracks/kin_thermo/results/thermo_interface.json',
        'tracks/kin_thermo/code/thermo.py','tracks/kin_thermo/code/formation.py',
        'tracks/kin_thermo/code/kernels.py','tracks/kin_thermo/code/objects.py',
        'tracks/kin_solve/code/ks_common.py','tracks/kin_solve/code/thermo_state.py','tracks/kin_solve/code/stage2_moments.py']
    pins=[]
    for rel in paths:
        path=BASE+'/'+rel;raw=subprocess.check_output(['git','show',COMMIT+':'+path],cwd=REPO)
        pins.append(dict(commit=COMMIT,path=path,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw)))
    save('SOURCES.json',dict(internal=pins,issue_handoff='https://github.com/ryangrvr/TestingGRUT/issues/13#issuecomment-6428006911',
        attribution='Claude completed build tracks and partial solver; Codex ports and conditional continuation. Same-session AI work, not external replication or observation.'))

if __name__=='__main__':main()
