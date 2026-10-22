"""v_scenario / v6: independent SN and BAO chi^2 for QV1 (Einstein-de Sitter at late times) and QV1-cal.

Data downloaded by this verifier (scratchpad v_scenario_dl; byte-identical to the scenario's copies):
  DESI DR2 BAO Gaussian likelihood, CobayaSampler/bao_data desi_bao_dr2 (arXiv:2503.14738): 13 numbers + cov
  Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov (arXiv:2202.04077)
Own choices: SN sample zHD > 0.01 (all rows; calibrator flags not used), m_b_corr with STAT+SYS covariance, absolute
magnitude marginalised analytically (profile); d_L = (1 + zHEL) D_M(zHD). BAO "shape": amplitude c/(H0 r_d)
profiled analytically; "absolute": r_d = 147.09 Mpc (Planck 2018 VI table 2 [R]; pre-recombination physics is
identical in QV1 and LCDM with the same omega_b, omega_c).
Output: ../results/v6_data.json
"""
import json
import math
import os

import numpy as np
from scipy.optimize import minimize_scalar

DL = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_scenario_dl'
C = 299792.458
OMR = 4.1842525717429235e-05     # omega_r from T0 and g_rho = 3.383 (v3)
OM = 0.1424


def Ez(z, Om, Or=0.0):
    return np.sqrt(Or * (1 + z) ** 4 + Om * (1 + z) ** 3 + (1 - Om - Or))


def DM_grid(Om, Or, zmax):
    zz = np.linspace(0, zmax, 400001)
    f = 1 / Ez(zz, Om, Or)
    cum = np.concatenate([[0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(zz))])
    return zz, cum


# ---------------------------------------------------------------- SN
names = open(os.path.join(DL, 'pplus.dat')).readline().split()
tab = np.genfromtxt(os.path.join(DL, 'pplus.dat'), names=True, dtype=None, encoding=None)
cov = np.loadtxt(os.path.join(DL, 'pplus_cov.cov'))
n = int(cov[0])
cov = cov[1:].reshape(n, n)
sel = tab['zHD'] > 0.01
zHD, zHEL, mb = tab['zHD'][sel], tab['zHEL'][sel], tab['m_b_corr'][sel]
Csel = cov[np.ix_(sel, sel)]
Ci = np.linalg.inv(Csel)
one = np.ones(len(mb))


def chi2_sn(Om, Or=0.0):
    zz, cum = DM_grid(Om, Or, zHD.max() + 0.01)
    dm = np.interp(zHD, zz, cum)                 # units c/H0
    mu = 5 * np.log10((1 + zHEL) * dm)           # + const absorbed in M
    r = mb - mu
    a = r @ Ci @ r
    b = one @ Ci @ r
    c = one @ Ci @ one
    return float(a - b * b / c)


# ---------------------------------------------------------------- BAO
rows = [l.split() for l in open(os.path.join(DL, 'desi_mean.txt')) if l.strip() and not l.startswith('#')]
zb = np.array([float(r[0]) for r in rows])
vb = np.array([float(r[1]) for r in rows])
kind = [r[2] for r in rows]
Cb = np.loadtxt(os.path.join(DL, 'desi_cov.txt'))
Cbi = np.linalg.inv(Cb)


def bao_shape_vec(Om, Or=0.0):
    zz, cum = DM_grid(Om, Or, 2.5)
    out = []
    for z, k in zip(zb, kind):
        dm = np.interp(z, zz, cum)
        dh = 1 / Ez(z, Om, Or)
        if k == 'DM_over_rs':
            out.append(dm)
        elif k == 'DH_over_rs':
            out.append(dh)
        elif k == 'DV_over_rs':
            out.append((z * dm * dm * dh) ** (1 / 3))
        else:
            raise ValueError(k)
    return np.array(out)       # in units of c/H0 (to be divided by r_d)


def chi2_bao_shape(Om, Or=0.0):
    f = bao_shape_vec(Om, Or)
    A = (f @ Cbi @ vb) / (f @ Cbi @ f)
    r = vb - A * f
    return float(r @ Cbi @ r), float(A)


def chi2_bao_abs(h, Om, Or, rd=147.09):
    f = bao_shape_vec(Om, Or) * C / (100 * h) / rd
    r = vb - f
    return float(r @ Cbi @ r)


out = dict(n_SN=int(sel.sum()), n_BAO=len(vb))
# QV1: EdS with radiation
h_qv1 = math.sqrt(OM + OMR)
Om_qv1, Or_qv1 = OM / h_qv1 ** 2, OMR / h_qv1 ** 2
h_cal = 0.6736418228945811
Om_cal, Or_cal = OM / h_cal ** 2, OMR / h_cal ** 2
sn_fit = minimize_scalar(lambda Om: chi2_sn(Om), bounds=(0.1, 0.6), method='bounded', options=dict(xatol=1e-6))
bao_fit = minimize_scalar(lambda Om: chi2_bao_shape(Om)[0], bounds=(0.1, 0.6), method='bounded',
                          options=dict(xatol=1e-6))
out['SN'] = dict(QV1=chi2_sn(Om_qv1, Or_qv1), QV1cal=chi2_sn(Om_cal, Or_cal), flatLCDM_best=float(sn_fit.fun),
                 flatLCDM_Om=float(sn_fit.x), EdS_no_radiation=chi2_sn(1.0))
out['SN']['dchi2_QV1'] = out['SN']['QV1'] - out['SN']['flatLCDM_best']
out['SN']['dchi2_QV1cal'] = out['SN']['QV1cal'] - out['SN']['flatLCDM_best']
out['BAO'] = dict(QV1_shape=chi2_bao_shape(Om_qv1, Or_qv1)[0], QV1cal_shape=chi2_bao_shape(Om_cal, Or_cal)[0],
                  flatLCDM_best=float(bao_fit.fun), flatLCDM_Om=float(bao_fit.x),
                  QV1_abs=chi2_bao_abs(h_qv1, Om_qv1, Or_qv1), QV1cal_abs=chi2_bao_abs(h_cal, Om_cal, Or_cal))
out['BAO']['dchi2_QV1_shape'] = out['BAO']['QV1_shape'] - out['BAO']['flatLCDM_best']
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v6_data.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
