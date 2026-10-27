"""Track pt_relax (Claude, own code): phi frequency action (number) from the saved shell spectra of the long runs.
N_phi/nQ = sum_shell E_phi(shell) / omega_phi(shell) / nQ, omega_phi^2 = k_shell^2 + m_phi^2 (Hartree, same sample).
X total action change is 2 x the X-bar change (charge fixed). Diagnostic [C] in the classical lattice model only."""
import os, json, glob, numpy as np
RUN = os.environ['RUNDIR']; OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
res = {}
for tag in ['IR12_s1_w1', 'IR12_fine_s1_w1']:
    rows, specs, nQ = [], [], None
    for f in sorted(glob.glob(os.path.join(RUN, tag, 'part*.npz'))):
        d = np.load(f); rows += json.loads(str(d['rows'])); specs += json.loads(str(d['specs']))
    nQ = rows[0]['Q']; mphi = {round(r['s'], 6): r['mphi2'] for r in rows}
    kk = np.arange(len(specs[0]['Ephi'])) * 2 * np.pi / 12.
    out = []
    for sp in specs:
        if abs(sp['s'] - round(sp['s'] / 50) * 50) > 1e-6: continue
        m2 = mphi[round(sp['s'], 6)]; E = np.array(sp['Ephi'])
        Nphi = float(np.sum(E[1:] / np.sqrt(kk[1:] ** 2 + m2)) / nQ)
        out.append(dict(s=sp['s'], N_phi_per_net=Nphi, mphi2=m2))
    res[tag] = out
    print(tag, [(o['s'], round(o['N_phi_per_net'], 4)) for o in out])
json.dump(res, open(os.path.join(OUT, 'phi_number.json'), 'w'), indent=1)
