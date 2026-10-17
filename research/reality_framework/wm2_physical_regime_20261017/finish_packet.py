"""Close the two-worker record; never launches or changes a trajectory."""
from pathlib import Path
import json
import summarize
import result_tables
import archive

BASE=Path(__file__).resolve().parent


def run():
    rows=[json.loads((BASE/f'results/S9_N5_root{i}.json').read_text()) for i in range(2)]
    s=summarize.run();assert s['numerical_branch_closed']
    result_tables.run()
    report=(BASE/'REPORT.md').read_text()
    old='**Analytic result; primary numerical workers still pending at this checkpoint.**'
    new='**Completed S9: an analytic operational regime is derived; the bounded finite recovery fails.**'
    assert s['A']['status']=='FAIL','Review report wording for any different scientific outcome'
    report=report.replace(old,new)
    before=('Exactly the two original five-event root components run at T=100, with\n'
        'rank 96 and a predeclared 1200 seconds each. No fallback or N=6 run is\n'
        'authorized. The new worker limit differs from S8\'s 900 and S7\'s 600;\n'
        'no equal-budget or equal-error superiority claim is made. A completion\n'
        'or stop will be recorded separately and cannot pass the physical ledger.')
    after=('Exactly the two original five-event root components were run at T=100,\n'
        'with rank 96 and a predeclared 1200 seconds each. Full original-weight\n'
        'recovery fails. No fallback or N=6 trajectory is launched; this numerical\n'
        'branch is closed. The worker limit differs from S8\'s 900 and S7\'s 600;\n'
        'no equal-budget or equal-error superiority claim is made.\n\n'
        '| Remaining counts | Status | Completed / target model time | Wall s | Peak MiB |\n'
        '|---|---|---:|---:|---:|\n')
    after+='\n'.join(f"| {r['nonroot_type_counts']} | {r['status']} | {r['completed_physical_steps_time']:.6f} / 100 | {r['whole_worker_wall_seconds']:.3f} | {r['peak_resident_MiB']:.3f} |" for r in rows)
    after+=('\n\nAll stop reasons and distinct error components are in RESULT_TABLES.md and\n'
        'the raw JSON. The recorded stops are runtime limits, not exhausted\n'
        'projection/error allowances. Neither is a negative formation observation.\n'
        'The production certificate passes its controls but its evaluated cost\n'
        'does not complete this finite recovery protocol. Another solver revision\n'
        'is not opened. The analytic theorem does not apply at N=5 and cannot\n'
        'substitute for the missing target-time numerical result.')
    if before in report:report=report.replace(before,after)
    else:assert 'this numerical\nbranch is closed' in report
    (BASE/'REPORT.md').write_text(report)
    methods=['# Numerical-method ledger','',
        'Production monotone-sweep residual implemented; all prerecovery controls pass.',
        'Full physical reference error .01359811 is inside .19647972.','',
        'The bounded recovery **fails** and this numerical branch is closed.',
        'Exactly the two frozen workers ran; no fallback, new allocation, larger',
        'budget or six-event trajectory follows. No target-time interval is',
        'inferred from a stopped state. See RESULT_TABLES.md and summary.json.','']
    for r in rows:
        methods.append(f"* Counts {r['nonroot_type_counts']}: {r['status']}, model t={r['completed_physical_steps_time']:.6f}, wall {r['whole_worker_wall_seconds']:.3f}s, peak {r['peak_resident_MiB']:.3f} MiB; {r['stop_reason']}.")
    methods+=['','The new account is mathematically checked and small-state validated,',
        'but the direct residual implementation is too costly for this frozen',
        'finite benchmark. Source pins, costs and failures are preserved.',
        'The separate all-N analytic proof neither validates the numerical',
        'method nor fills the missing five-event target. S5 remains the completed',
        'finite-formation frontier; S9 adds a distinct asymptotic derivation.','']
    (BASE/'METHODS_LEDGER.md').write_text('\n'.join(methods).rstrip()+'\n')
    archive.run();return s


if __name__=='__main__':print(json.dumps(run(),indent=2))
