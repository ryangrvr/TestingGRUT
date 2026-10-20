"""Verify saved diagnostics and summarize the three bounded follow-up studies."""
from pathlib import Path
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text())


def window_summary(data, start):
    rows = [r for r in data['trace'] if r['time'] >= start - 1e-7]
    result = {}
    for side in ['left', 'right', 'central']:
        selected = []
        for row in rows:
            balls = [b for b in row['balls'] if b['radius'] == 12 and not b['overlaps_other_ball']]
            if side == 'left':
                balls = [b for b in balls if b['center_z'] < -1]
            elif side == 'right':
                balls = [b for b in balls if b['center_z'] > 1]
            else:
                balls = [b for b in balls if abs(b['center_z']) <= 1]
            if balls:
                selected.append(max(balls, key=lambda b: abs(b['Qbar'])))
        if selected:
            result[side] = dict(samples=len(selected), Qbar_min=min(x['Qbar'] for x in selected),
                                Qbar_max=max(x['Qbar'] for x in selected),
                                E_per_Q_min=min(x['energy_per_absolute_charge'] for x in selected),
                                E_per_Q_max=max(x['energy_per_absolute_charge'] for x in selected),
                                first_center=selected[0]['center_z'], last_center=selected[-1]['center_z'])
    return result


def compare(a, b, cutoff=None, start=0):
    bmap = {round(r['time'], 7): r for r in b['trace']}
    dz = dmiddle = dq = deq = 0.0
    mode_differences = [0.0, 0.0]
    for ra in a['trace']:
        if ra['time'] < start - 1e-7:
            continue
        if cutoff is not None and ra['time'] > cutoff + 1e-7:
            continue
        rb = bmap.get(round(ra['time'], 7))
        if rb is None:
            continue
        dz = max(dz, abs(ra['z_rms'] - rb['z_rms']))
        dmiddle = max(dmiddle, abs(ra['middle_norm_fraction'] - rb['middle_norm_fraction']))
        for ba in [x for x in ra['balls'] if x['radius'] == 12]:
            bs = [x for x in rb['balls'] if x['radius'] == 12 and abs(x['center_z'] - ba['center_z']) < 3]
            if bs:
                bb = min(bs, key=lambda x: abs(x['center_z'] - ba['center_z']))
                dq = max(dq, abs(ba['Qbar'] - bb['Qbar']))
                deq = max(deq, abs(ba['energy_per_absolute_charge'] - bb['energy_per_absolute_charge']))
        for m, (ma, mb) in enumerate(zip(ra.get('modes', []), rb.get('modes', []))):
            mode_differences[m] = max(mode_differences[m], abs(ma['amplification'] - mb['amplification']))
    return dict(max_z_rms_difference=dz, max_middle_fraction_difference=dmiddle,
                max_matched_ball_charge_difference=dq, max_matched_ball_EperQ_difference=deq,
                max_mode_amplification_differences=mode_differences)


def main():
    checks = []
    def check(name, passed, detail=None):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))

    for dirname in ['', 'formation', 'record']:
        root = ROOT / dirname
        manifest = load(root / 'FREEZE_MANIFEST.json')
        for relative, digest in manifest.items():
            actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
            check('Frozen source ' + relative, actual == digest)
    origin = load(ROOT / 'kernel_origin.json')
    check('Archived physical kernel unchanged', hashlib.sha256((ROOT.parent / 'kernel.cjs').read_bytes()).hexdigest() == origin['parent_sha256'])
    check('Tangent implementation controls', load(ROOT / 'controls_preflight.json')['pass'])
    check('Gaussian setup and ablation controls', load(ROOT / 'formation/controls.json')['pass'])
    studies = {}
    all_data = {}
    for dirname, label in [('', 'persistence'), ('formation', 'formation'), ('record', 'record')]:
        folder = ROOT / dirname
        protocol = load(folder / 'protocol.json')
        datasets = {}
        for spec in protocol['runs']:
            file = folder / 'results' / (spec['id'] + '.json')
            check(label + ': ' + spec['id'] + ' exists', file.exists())
            if not file.exists():
                continue
            data = load(file)
            datasets[data['id']] = data
            energy = data['maximum_relative_energy_drift']
            charge = data['maximum_charge_error']
            finite = all(all(math.isfinite(v) for v in row.values() if isinstance(v, (int, float))) for row in data['trace'])
            mode_finite = all(all(math.isfinite(m['norm_squared']) for m in row.get('modes', [])) for row in data['trace'])
            check(label + ': ' + spec['id'] + ' numerical consistency', finite and mode_finite and energy < .001 and charge < 1e-9,
                  dict(energy_drift=energy, charge_error=charge))
            if label == 'record':
                check('Record/control decoder ' + spec['id'], data['readout']['all_samples_match'], data['readout'])
        studies[label] = dict(planned=len(protocol['runs']), completed=len(datasets),
                              max_energy_drift=max(d['maximum_relative_energy_drift'] for d in datasets.values()),
                              max_charge_error=max(d['maximum_charge_error'] for d in datasets.values()))
        all_data[label] = datasets

    persistence = {}
    for name, data in all_data['persistence'].items():
        mode_end = data['trace'][-1]['modes']
        persistence[name] = dict(final_time=data['trace'][-1]['time'],
                                 window=window_summary(data, 160 if data['configuration']['T'] == 240 else 120),
                                 final_modes=mode_end,
                                 max_mode_amplifications=[max(row['modes'][i]['amplification'] for row in data['trace']) for i in range(2)],
                                 linear_warning=data['first_linear_warning_time'])
    comparisons = {}
    p = all_data['persistence']
    for name in ['single_static', 'pair_phase_0', 'pair_phase_halfpi']:
        comparisons['coarse_medium_' + name] = compare(p['coarse_' + name], p['medium_' + name], 240)
        comparisons['medium_fine_' + name] = compare(p['medium_' + name], p['fine_' + name], 160)
    comparisons['half_timestep'] = compare(p['medium_pair_phase_halfpi'], p['medium_half_dt_pair_phase_halfpi'], 240)
    comparisons['larger_box_phase0'] = compare(p['medium_pair_phase_0'], p['medium_larger_box_pair_phase_0'], 240)
    comparisons['larger_box_halfpi'] = compare(p['medium_pair_phase_halfpi'], p['medium_larger_box_pair_phase_halfpi'], 240)
    comparisons['late_half_timestep'] = compare(p['medium_pair_phase_halfpi'], p['medium_half_dt_pair_phase_halfpi'], 240, 160)
    comparisons['late_medium_fine_halfpi'] = compare(p['medium_pair_phase_halfpi'], p['fine_pair_phase_halfpi'], 160, 120)

    formation = {}
    for name, data in all_data['formation'].items():
        final = data['trace'][-1]
        balls = [b for b in final['balls'] if b['radius'] == 12]
        core = max(balls, key=lambda x: abs(x['Qbar'])) if balls else None
        formation[name] = dict(final_time=final['time'], maxF=final['max_F'], minS=final['min_S'],
                               core12=core, initial_EperQ=data['trace'][0]['Ebar'] / data['trace'][0]['Qbar'],
                               core_fraction_of_total_charge=core['Qbar'] / final['Qbar'] if core else None,
                               last_window=window_summary(data, 160 if final['time'] == 240 else 120))
    formation_comparison = compare(all_data['formation']['medium_sigma6'], all_data['formation']['medium_larger_box_sigma6'], 240)
    records = {name: data['readout'] for name, data in all_data['record'].items()}
    record_protocol = load(ROOT / 'record/protocol.json')
    for name, data in all_data['record'].items():
        window = [row for row in data['trace'] if row['time'] >= data['readout']['start'] - 1e-7]
        strict = []
        for row in window:
            two_sided = any(x['z'] < -1 for x in row['centers']) and any(x['z'] > 1 for x in row['centers'])
            balls = [b for b in row['balls'] if b['radius'] == record_protocol['receiver_radius'] and not b['overlaps_other_ball'] and b['center_z'] > 1]
            q = max(balls, key=lambda b: b['center_z'])['Qbar'] if two_sided and balls else None
            threshold = record_protocol['reference_charge']
            margin = threshold * record_protocol['undecided_fraction']
            decoded = None if q is None else (1 if q > threshold + margin else (-1 if q < threshold - margin else 0))
            strict.append(decoded == data['readout']['target'])
        records[name]['strict_separated_receiver_pass'] = all(strict)
        records[name]['receiver_review'] = 'Post-run clarification: require distinct left and right regions, and nonoverlapping receiver ball. The frozen raw code can label a central component before the readout window. All registered window rows satisfy the stricter condition; no evolution or registered-window outcome changed.'
        check('Strict separated receiver ' + name, all(strict))
    result = dict(studies=studies, all_numerical_and_source_checks_pass=all(x['passed'] for x in checks), checks=checks,
                  persistence=persistence, persistence_comparisons=comparisons,
                  formation=formation, formation_domain_comparison=formation_comparison, records=records,
                  interpretation='Conditional scalar-model mechanisms and finite diagnostics. Selected m=1 response is resolution-sensitive; no full nonlinear non-axisymmetric stability, generic cosmology, implemented detector or consciousness result.')
    (ROOT / 'SUMMARY.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(studies=studies, checks=len(checks), all_pass=result['all_numerical_and_source_checks_pass']), indent=2))
    if not result['all_numerical_and_source_checks_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
