"""Adapted summarize.py for one case (reproduction): rebuild the RESULTS.json row of <tag> from the regenerated
RUN/WINDOW files with summarize.py's own field mapping, and compare with the published row.
Usage: python3 -B summarize_one.py <packet_dir (regenerated)> <published RESULTS.json> <tag>"""
import sys, json
pk, pubres, tag = sys.argv[1:4]
r = json.load(open(f'{pk}/RUN_{tag}.json')); w0 = json.load(open(f'{pk}/WINDOW_{tag}_14_24.json')); w1 = json.load(open(f'{pk}/WINDOW_{tag}_24_34.json'))
row = dict(case=r['case'], initial_q_RMS_over_mean=r['initial'].get('q_RMS_over_mean'), windows=[])
for w in [w0, w1]:
    row['windows'].append(dict(s=w['window_s'], tau=w['window_tau'], coherent_negative_action=w['coherent_negative_frequency_action'],
        wave_negative_action=w['wave_negative_frequency_action'], total_negative_action=w['total_negative_frequency_action'],
        spectral_charge_ratio=w['spectral_charge_ratio'], reference_mass_means={m: w['equal_time'][f'ref_{m}']['hann2_mean'] for m in ['0', '7', '12', '20', '26', 'H']},
        omitted_Hartree_reference_above8=w['omitted_k_gt8_Hartree_negative_reference'],
        fundamental_shell_wave_fraction=w['fundamental_shell_wave_fraction'], low_wave_fraction_kcenters_le1p6=w['low_wave_fraction_kcenters_le1p6']))
pub = next(x for x in json.load(open(pubres))['rows'] if x['case']['tag'] == tag)
worst = 0.
def walk(a, b, p=''):
    global worst
    if isinstance(a, dict):
        for k in a: walk(a[k], b[k], p + '/' + k)
    elif isinstance(a, list):
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, p + f'[{i}]')
    elif isinstance(a, float):
        worst = max(worst, abs(a - b) / max(abs(a), 1e-300))
walk(pub, row)
line = ' | '.join(f"{w['coherent_negative_action']:.4f} + {w['wave_negative_action']:.4f} = {w['total_negative_action']:.4f}" for w in row['windows'])
print(json.dumps(dict(tag=tag, max_relative_difference_vs_published_row=worst, table_line=line, row_windows=row['windows']), indent=1))
