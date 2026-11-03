"""e2e_card: independent spot check of e2e_sm's non-thermal leptogenesis (eta_B per unit delta).

Own minimal model (not e2e_sm's code):
  - source: every unit of condensate energy lost at rate Gamma_cond becomes N4 pairs (chi -> N4 N4, E_pair = m_chi);
    N4 decays at once (Gamma_N4 = 0.12 H_i >> H late), each decay gives eps4 units of L.
    Only rows with f_sig = 0 (sigma condensate gone) are used for the chi -> N4 attribution [A].
  - washout: inverse decays of all four N, Gamma_W = (1/4) sum_j Gamma_j z_j^2 K1(z_j), z_j = M_j / T_SM
    (BDP hep-ph/0401240 W_ID = K z^3 K1 / 4 rewritten as a rate), plus Delta L = 2,
    Gamma = 0.186 mbar^2 (1.66 sqrt(g*)) T^3 / (1e10 GeV M_Pl) (BDP's W_DeltaL=2 rewritten as a rate).
  - eta_B = 7.04 * (28/79) * |Y_L| with s = (2 pi^2/45) g* T_SM^3, g* = 106.75.
Inputs: any fluid_*.json with the e2e_kinetic column format (t, a, H, f_sig, f_chi, T_SM, Gamma_cond).
Everything is a model choice; [C] on [E] rates.
"""
import json, sys, math
import numpy as np
from scipy.special import k1

MP = 2.435e18  # reduced Planck mass, GeV
GSTAR = 106.75
M = np.array([250.0, 400.0, 200.0, 120.0])          # H_i units (mf_hybrid)
GAM = np.array([0.5333, 2.4906, 0.3413, 0.12287])   # rest widths, H_i units (e2e_sm S1; N2 at mt = 0.0912 eV)
EPS = np.array([2.069e-2, 3.167e-3, 9.613e-3, 6.674e-3])  # per unit delta (e2e_sm S1), N4 = index 3
MCHI = 249.995
MBAR2 = 2.588e-3  # eV^2, sum m_i^2 for NO, m1 = 0


def run(path, Hi_GeV=None, wash_scale=1.0, eps4=EPS[3]):
    d = json.load(open(path))
    h = d['history']
    Hi = Hi_GeV or d.get('H_i_GeV', 1.3e11)
    t = np.array(h['t'], float); a = np.array(h['a'], float); H = np.array(h['H'], float)
    fs = np.array([x or 0.0 for x in h['f_sig']]); fc = np.array([x or 0.0 for x in h['f_chi']])
    T = np.array([x or 0.0 for x in h['T_SM']], float)
    # comoving condensate energy in units M_P^2 H_i^2 * (comoving volume)
    Ec = (fs + fc) * 3.0 * H**2 * a**3
    if 'Gamma_cond' in h:
        gc = np.array([x if x is not None else 0.0 for x in h['Gamma_cond']], float)
        src_E = gc * Ec  # comoving energy loss rate of the condensates
    else:  # cycle-mean condensates redshift as matter: loss rate = -d(Ec)/dt
        src_E = np.maximum(-np.gradient(Ec, t), 0.0)
    late = fs <= 1e-12
    srcN4 = np.where(late, 2.0 * src_E / MCHI, 0.0)  # comoving N4 number production rate
    # washout rate (H_i units)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        GW = np.zeros_like(T)
        for j in range(4):
            z = np.where(T > 0, M[j] / np.maximum(T, 1e-30), np.inf)
            GW += 0.25 * GAM[j] * np.where(z < 700, z**2 * k1(np.minimum(z, 700)), 0.0)
        cdl2 = 0.186 * MBAR2 * 1.66 * math.sqrt(GSTAR) * Hi**2 / (1e10 * 1.22e19)  # Gamma/H_i per That^3
        GW += cdl2 * T**3
    GW *= wash_scale
    # exponential integrator for dL/dt = eps4 * S - GW L   (L in units M_P^2 H_i^2 comoving / H_i)
    L = 0.0; Lnw = 0.0
    contrib = []
    for i in range(len(t) - 1):
        dt = t[i + 1] - t[i]
        g = 0.5 * (GW[i] + GW[i + 1]); S = 0.5 * (srcN4[i] + srcN4[i + 1]) * eps4
        if g * dt > 1e-12:
            L = L * math.exp(-g * dt) + S / g * (1 - math.exp(-g * dt))
        else:
            L = L + S * dt
        Lnw += S * dt
    # convert to yield: L * (MP/Hi)^2 [H_i^4 units comoving] / (s a^3) with s in H_i^3 units
    sfin = (2 * math.pi**2 / 45) * GSTAR * T[-1]**3 * a[-1]**3
    # L is in M_P^2 H_i units (comoving); n/H_i^3 needs a factor (M_P/H_i)^2
    YL = L * (MP / Hi)**2 / sfin
    YLnw = Lnw * (MP / Hi)**2 / sfin
    etaB = 7.04 * (28 / 79) * YL
    # source-weighted temperature and survival
    w = srcN4 * np.gradient(t)
    Tw = float(np.sum(w * T) / np.sum(w)) if np.sum(w) > 0 else None
    return dict(file=path.split('/')[-1], H_i_GeV=Hi, wash_scale=wash_scale, eps4=eps4,
                YN4_nonthermal=float(np.sum(w) * (MP / Hi)**2 / sfin),
                etaB=float(etaB), etaB_over_eps4=float(etaB / eps4), survival=float(YL / YLnw) if YLnw else None,
                etaB_per_delta=float(etaB / eps4 * EPS[3]), T_src_weighted_GeV=Tw * Hi if Tw else None,
                K=float(a[-1]**2 * H[-1]))


if __name__ == '__main__':
    out = []
    for p in sys.argv[1:]:
        Hi = None
        if 'H84' in p:
            Hi = 8.4e10
        for ws in (1.0, 0.5, 2.0):
            r = run(p, Hi_GeV=Hi, wash_scale=ws)
            out.append(r)
            print(json.dumps(r))
    json.dump(out, open('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/e2e_card/lepto_check.json', 'w'), indent=1)
