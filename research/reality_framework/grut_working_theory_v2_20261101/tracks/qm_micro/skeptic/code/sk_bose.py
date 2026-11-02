"""qm_micro SKEPTIC (Claude, 2 Nov 2026): Bose enhancement of the linear exchange coefficient.  Own script; uses
qm_micro's flux() (read-only import) and its coexistence chemical potentials.
qm_micro's coarsen.py writes the off-coexistence exchange as Phi(mu_c) [e^{d/T} - e^{g}] (exponential mu scaling,
exact for Boltzmann).  Its linear coefficient is Phi/T, whereas detailed balance with Bose weights gives
dPhi/dmu = int w s f(1+f)/T (the L of qm_micro's own K_W).  B = T (dPhi/dmu)/Phi is the missing factor, per species."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
import qm_vapour as qv, bessel_vapour_copy as bv, thermo_copy as th
rows = []
for T in np.geomspace(0.002, 12.0, 40):
    T = float(T)
    mc = bv.mu_coex(T) if T < 0.5 else th.mu_coex(T)
    FX, dFX = qv.flux(mc, T, 'X', True, True); FB, dFB = qv.flux(-mc, T, 'Xbar', True, True)
    rows.append(dict(T=T, mu_c=mc, B_X=T * dFX / FX, B_Xbar=(T * dFB / FB) if FB > 0 else 1.0,
                     coded_over_true_L=(FX + FB) / T / (dFX + dFB)))
    print('T=%.4g B_X=%.4f B_Xbar=%.4f coded/true L=%.4f' % (T, rows[-1]['B_X'], rows[-1]['B_Xbar'], rows[-1]['coded_over_true_L']))
json.dump(dict(notes=__doc__, rows=rows), open(sys.argv[1], 'w'), indent=1)
