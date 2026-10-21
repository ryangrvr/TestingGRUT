"""Kibble-Zurek consistency, quench time, thermalisation threshold, zero-T check, percolation spot check.
Independent of the module.  Writes results/kz_therm_v.json."""
import json
import os
import numpy as np
from scipy import special, ndimage
from scipy.optimize import brentq, minimize_scalar
import vt_core as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
PI = np.pi
MP = 2.435323e18
GS = 106.75
VGEV = 1000.0
res = {}

# ---------------------------------------------------------------- zero-T check
P = V.BENCH


def V0T(phi):
    out = V.V0_phi(np.array([phi]), P)
    for n, m2, _ in V.masses_phi(np.array([phi]), P):
        out = out + n * V.cw(m2, V.MU0)
    return float(out[0])


r = minimize_scalar(V0T, bracket=(0.9, 1.0, 1.1), tol=1e-14)
h = 1e-4
res['zeroT'] = dict(phi_min=r.x, curvature=(V0T(r.x + h) - 2 * V0T(r.x) + V0T(r.x - h)) / h**2)
print(res['zeroT'])

# ---------------------------------------------------------------- DR kappa and T_c
tr = json.load(open(os.path.join(OUT, 'transition_v.json')))
Tc_DR = tr['running']['mu0']['DR']['Tc']
Tc_AE = tr['running']['mu0']['AE']['Tc']


def curvDR(T):
    return 2 * float(V.dVds(np.array([0.0]), T, P, V.MU0, 'DR')[0])


dT = 1e-6
dcdT = (curvDR(Tc_DR + dT) - curvDR(Tc_DR - dT)) / (2 * dT)
kappa = dcdT / (2 * Tc_DR)            # m^2 = 2 kappa Tc^2 eps
res['DR'] = dict(Tc=Tc_DR, kappa=kappa, slope_over_Pi=kappa / V.debye(P)['phi'], xi0_Tc=1 / np.sqrt(2 * kappa))
lp = P['lp']
tG = (3 * lp / (4 * PI))**2 / (2 * kappa)
res['DR']['t_G'] = tG
print(res['DR'])

# ---------------------------------------------------------------- Hubble, tau_Q, KZ
TcG = Tc_DR * VGEV


def H_of(Tsm, g=GS):
    return np.sqrt(PI**2 * g / 90) * Tsm**2 / MP


nu_I = 1 / (3 - 1.412625)
res['nu'] = nu_I


def kz(xi, z, mode, tau0_fac=1.0, Tc=TcG, kap=kappa):
    tauQ = 1 / H_of(Tc / xi)
    xi0 = 1 / (Tc * np.sqrt(2 * kap))
    tau0 = tau0_fac * xi0
    if mode == 'naive':
        a = nu_I / (1 + nu_I * z)
        eps = (tau0 / tauQ)**(1 / (1 + nu_I * z))
        xihat = xi0 * eps**(-nu_I)
    else:
        # matched: xi = xi0 eps^-1/2 (eps>=tG), xi0 tG^-1/2 (eps/tG)^-nu (eps<tG); tau = tau0 (xi/xi0)^z
        def corr(e):
            return xi0 * e**-0.5 if e >= tG else xi0 * tG**-0.5 * (e / tG)**(-nu_I)

        def f(le):
            e = np.exp(le)
            return np.log(e * tauQ) - np.log(tau0 * (corr(e) / xi0)**z)
        le = brentq(f, np.log(1e-40), np.log(1.0))
        eps = np.exp(le)
        xihat = corr(eps)
    NH = (tauQ / xihat)**3
    return dict(tauQ_GeVinv=tauQ, H=1 / tauQ, eps_hat=eps, xi_hat_GeVinv=xihat, log10_NH=np.log10(NH))


kzrows = []
for xi in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    row = dict(xi=xi)
    for z in [1, 2]:
        for mode in ['matched', 'naive']:
            row[f'z{z}_{mode}'] = kz(xi, z, mode)
    kzrows.append(row)
    print(xi, row['z1_matched']['H'], row['z1_matched']['log10_NH'], row['z1_naive']['log10_NH'],
          row['z2_matched']['log10_NH'], row['z2_naive']['log10_NH'], row['z1_matched']['xi_hat_GeVinv'],
          row['z1_matched']['eps_hat'])
res['kz'] = kzrows
# closed form check of naive formula and slopes
for z in [1, 2]:
    a = nu_I / (1 + nu_I * z)
    tauQ = 1 / H_of(TcG)
    xi0 = 1 / (TcG * np.sqrt(2 * kappa))
    closed = 3 * np.log10(tauQ**(1 - a) * xi0**a / xi0)
    res[f'naive_closed_form_xi1_z{z}'] = closed
    res[f'slope_z{z}'] = 6 * (1 - a)
    print('closed', z, closed, 6 * (1 - a))
# tau0 sensitivity at xi = 1
res['tau0_band_xi1'] = {f'z{z}_{m}': [kz(1.0, z, m, f)['log10_NH'] for f in (1e-2, 1, 1e2)] for z in (1, 2) for m in ('matched', 'naive')}
# collision-limited tau0 ~ 1/Gamma_phi (alternative, not in module's band)
zeta3 = special.zeta(3)
n_over_T3 = zeta3 / PI**2
Esum = 18 * lp**2 + 12 * P['c']**2 + 6 * P['d']**2
s_avg = 2 * (PI**4 / (30 * zeta3))**2
Gam_mod = n_over_T3 * Esum / (16 * PI * s_avg)       # module-style estimate, Gamma/T
tau0_coll_fac = (1 / (Gam_mod * TcG)) / (1 / (TcG * np.sqrt(2 * kappa)))
res['tau0_collision_factor'] = tau0_coll_fac
res['kz_xi1_tau0_collision'] = {f'z{z}_{m}': kz(1.0, z, m, tau0_coll_fac)['log10_NH'] for z in (1, 2) for m in ('matched', 'naive')}
print('tau0 coll fac', tau0_coll_fac, res['kz_xi1_tau0_collision'])
res['hidden_energy_H_factor_xi1'] = np.sqrt(1 + 4 / GS)

# ---------------------------------------------------------------- thermalisation threshold (phi kinetic)
# (a) module-style kinematics (reproduce): sigma v ~ |M|^2 S/(16 pi <s>), <s> = 2<E>^2
# (b) my kinematics: massless contact, sigma v_Mol = |M|^2 S/(32 pi E1 E2) exactly (angle independent);
#     thermal average <1/E> for Bose (pi^2/6)/(2 zeta3)/T, Boltzmann 1/(2T)
inv_E_bose = (PI**2 / 6) / (2 * zeta3)
inv_E_boltz = 0.5
Gam_b_bose = n_over_T3 * Esum / (32 * PI) * inv_E_bose**2
Gam_b_typ = n_over_T3 * Esum / (32 * PI) * inv_E_bose / (PI**4 / (30 * zeta3))  # particle at <E>
Gam_b_boltz = (1 / PI**2) * Esum / (32 * PI) * inv_E_boltz**2                   # Boltzmann n = T^3/pi^2


def ximin(Gam_over_T, Tc=TcG):
    return np.sqrt(H_of(Tc) / (Gam_over_T * Tc))


res['therm'] = dict(coupling_sum=Esum, Gamma_over_T_module_style=Gam_mod, ximin_module_style=ximin(Gam_mod),
                    Gamma_over_T_bose_avg=Gam_b_bose, ximin_bose_avg=ximin(Gam_b_bose),
                    Gamma_over_T_typical=Gam_b_typ, ximin_typical=ximin(Gam_b_typ),
                    Gamma_over_T_boltz=Gam_b_boltz, ximin_boltz=ximin(Gam_b_boltz),
                    H_at_Tc_GeV=H_of(TcG),
                    Gamma_D_exch_over_T=n_over_T3 * 6 * P['d']**2 / (16 * PI * s_avg),
                    ximin_D_exch=ximin(n_over_T3 * 6 * P['d']**2 / (16 * PI * s_avg)),
                    Gamma_X_over_T=n_over_T3 * (24 * P['lx']**2 + 6 * P['c']**2) / (16 * PI * s_avg),
                    ximin_X=ximin(n_over_T3 * (24 * P['lx']**2 + 6 * P['c']**2) / (16 * PI * s_avg)))
print(res['therm'])
# fugacity factors
for zf in [0.9, 0.5, 0.1]:
    res[f'fugacity_Tc_factor_z{zf}'] = np.sqrt(special.zeta(2) / float(np.real(special.spence(1 - zf))))
print({k: v for k, v in res.items() if k.startswith('fugacity')})

# ---------------------------------------------------------------- percolation spot check
rng = np.random.default_rng(12345)
perc = []
for L in [32, 64]:
    for seed in range(3):
        rng = np.random.default_rng(1000 + seed)
        s = rng.random((L, L, L)) < 0.5
        lab_p, n_p = ndimage.label(s)
        lab_m, n_m = ndimage.label(~s)
        sz_p = np.bincount(lab_p.ravel())[1:]
        sz_m = np.bincount(lab_m.ravel())[1:]
        bp = np.argmax(sz_p) + 1
        bm = np.argmax(sz_m) + 1

        def spans(lab, b):
            m = lab == b
            return all(m.take(0, axis=a).any() and m.take(L - 1, axis=a).any() for a in range(3))
        # walls: unlike nearest-neighbour pairs (open boundaries)
        tot = 0
        sep = 0
        for a in range(3):
            A = [slice(None)] * 3
            B = [slice(None)] * 3
            A[a] = slice(0, L - 1)
            B[a] = slice(1, L)
            sa, sb = s[tuple(A)], s[tuple(B)]
            unlike = sa != sb
            tot += unlike.sum()
            pa = np.where(sa, lab_p[tuple(A)], -lab_m[tuple(A)])
            pb = np.where(sb, lab_p[tuple(B)], -lab_m[tuple(B)])
            sep += (unlike & (((pa == bp) & (pb == -bm)) | ((pa == -bm) & (pb == bp)))).sum()
        perc.append(dict(L=L, seed=seed, frac_plus_in_largest=sz_p.max() / s.sum(),
                         frac_minus_in_largest=sz_m.max() / (~s).sum(), span_plus=spans(lab_p, bp),
                         span_minus=spans(lab_m, bm), wall_frac_between_spanning=sep / tot))
        print(perc[-1])
res['percolation'] = perc

with open(os.path.join(OUT, 'kz_therm_v.json'), 'w') as f:
    json.dump(res, f, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
