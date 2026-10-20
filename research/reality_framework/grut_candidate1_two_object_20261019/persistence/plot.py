"""Generate the scientific figure from saved diagnostic outputs."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'svg.hashsalt': 'grut-persistence-20261019'})


def read(relative):
    return json.loads((ROOT / relative).read_text())


def ball(row):
    choices = [b for b in row['balls'] if b['radius'] == 12]
    return min(choices, key=lambda b: abs(b['center_z'])) if choices else None


fig, axes = plt.subplots(2, 2, figsize=(13, 8.8), constrained_layout=True)
fig.suptitle('Candidate 1: assembly, persistence and a conditional charge record', fontsize=16)
colors = ['#2563eb', '#d97706', '#16a34a', '#6b7280']
for name, label, color in zip(['medium_sigma4', 'medium_sigma6', 'medium_sigma8', 'medium_rigid_sigma6'],
                              ['width 4', 'width 6', 'width 8', 'width 6, rigid field'], colors):
    data = read('formation/results/' + name + '.json')['trace']
    t = [r['time'] for r in data]
    style = '--' if 'rigid' in name else '-'
    axes[0, 0].plot(t, [r['max_F'] for r in data], style, lw=1.5, color=color, label=label)
    axes[0, 1].plot(t, [ball(r)['energy_per_absolute_charge'] if ball(r) else np.nan for r in data],
                    style, lw=1.5, color=color, label=label)
axes[0, 0].set(title='A. Assembly from prepared Gaussian clouds', xlabel=r'$\tau=gvt$', ylabel=r'Maximum $F=\sqrt{u^2+w^2}$')
axes[0, 0].legend(frameon=False, ncol=2, fontsize=9)
axes[0, 1].axhline(1, color='#111827', ls=':', lw=1, label='free-carrier value')
axes[0, 1].set(title='B. Energy per charge in the measured core', xlabel=r'$\tau$', ylabel=r'$\bar E_{R=12}/|\bar Q_{R=12}|$', ylim=(.85, 1.13))
axes[0, 1].text(.98, .08, 'Finite moving regions; not isolated subsystems', transform=axes[0, 1].transAxes, ha='right', fontsize=8, color='#4b5563')

record = read('record/protocol.json')
ax = axes[1, 0]
ax.axhspan(.9, 1.1, color='#e5e7eb', label='undecided range')
ax.axvspan(160, 240, color='#dcfce7', alpha=.5)
for name, label, color in [('medium_plus', r'$+\pi/2$, encounter', '#2563eb'), ('medium_minus', r'$-\pi/2$, encounter', '#d97706'), ('medium_miss', 'missed encounter', '#6b7280')]:
    data = read('record/results/' + name + '.json')['trace']
    values = []
    for row in data:
        separated = any(x['z'] < -1 for x in row['centers']) and any(x['z'] > 1 for x in row['centers'])
        choices = [b for b in row['balls'] if b['radius'] == 12 and not b['overlaps_other_ball'] and b['center_z'] > 1]
        values.append(max(choices, key=lambda b: b['center_z'])['Qbar'] / record['reference_charge'] if separated and choices else np.nan)
    ax.plot([r['time'] for r in data], values, color=color, lw=1.8, label=label)
ax.set(title='C. Phase input leaves distinguishable charge', xlabel=r'$\tau$', ylabel=r'Right-core charge / initial object charge', xlim=(80, 240), ylim=(.5, 1.6))
ax.legend(frameon=False, fontsize=8, loc='upper center', ncol=2)
ax.text(.98, .04, 'Green: tested readout window; ideal decoder', transform=ax.transAxes, ha='right', fontsize=8, color='#4b5563')

ax = axes[1, 1]
for level, h, style, width in [('coarse', '.5', ':', 1.3), ('medium', '1/3', '-', 1.8), ('fine', '.25', '--', 1.8)]:
    rows = read('results/' + level + '_pair_phase_halfpi.json')['trace']
    ax.plot([r['time'] for r in rows], [r['modes'][0]['amplification'] for r in rows], style,
            color='#2563eb', lw=width, label=r'$m=1$, $h=' + h + '$')
for level, style in [('medium', '-'), ('fine', '--')]:
    rows = read('results/' + level + '_pair_phase_halfpi.json')['trace']
    ax.plot([r['time'] for r in rows], [r['modes'][1]['amplification'] for r in rows], style,
            color='#d97706', lw=1.5, label=r'$m=2$, ' + level)
ax.set(title='D. Linear azimuthal response needs convergence', xlabel=r'$\tau$', ylabel='Canonical-norm amplification', yscale='log')
ax.legend(frameon=False, fontsize=8, ncol=2)
for ax in axes.flat:
    ax.grid(alpha=.18)
    ax.set_xlim(left=0 if ax is not axes[1, 0] else 80, right=240)
svg = ROOT / 'RESULTS.svg'
fig.savefig(svg, metadata={'Date': None})
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--qa', type=Path, default=ROOT / 'RESULTS.png')
qa = parser.parse_args().qa
qa.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(qa, dpi=145)
print(qa)
