#!/usr/bin/env python3
"""e2e_sm: UV-completion radial-mode test, Delta N_eff of N_s dark radiation, sphaleron/ BBN timing,
N_s soft-mass variant estimate  -> results/checks.json.  All [C] on stated inputs or [E]."""
import json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KIN = os.path.join(HERE, '..', '..', 'e2e_kinetic', 'results')
RES = os.path.join(HERE, '..', 'results')
H_GEV = 1.3e11
MP = 2.435e18 / H_GEV
GST = 106.75
out = {}

d = json.load(open(os.path.join(KIN, 'fluid_HFa.json')))
h = d['history']
t = np.array(h['t']); a = np.array(h['a']); H = np.array(h['H'])
A1 = np.array([x or 0.0 for x in h['A_N1']]); A2 = np.array([x or 0.0 for x in h['A_N2']])
fcond = np.array(h['f_sig']) + np.array(h['f_chi'])

# ---- 1. UV completion: Phi_a = (rho_a/sqrt2) e^{i theta_a}, sigma = f1 theta_1, chi = f2 theta_2 ----
# theta-dot = 2 mu5 = 2 A m_parent:  sigma: 2*600/1200 A1 = A1 ; chi: 2*250/2000 A2 = A2/4
u_sig = 2 * 599.99 / 1200.0 * A1
u_chi = 2 * 249.995 / 2000.0 * A2
uv = {}
for name, u in [('sigma_f1200', u_sig), ('chi_f2000', u_chi)]:
    umax = float(np.max(u))
    rec = dict(theta_dot_over_f_max=umax,
               lambda_needed_m_rho_ge_2thetadot=2 * umax**2,
               lambda_needed_drho_over_f_le_0p1=umax**2 / 0.21)
    for lam in [1.0, 4 * math.pi, (4 * math.pi)**2 / 2]:
        mrho = math.sqrt(2 * lam)
        ok = np.where(u <= mrho / 2)[0]
        # first time after the maximum when theta-dot <= m_rho/2 stays true
        imax = int(np.argmax(u))
        idx = [i for i in ok if i > imax]
        t_ok = float(t[idx[0]]) if idx else None
        x = 0.5 * (1 + math.sqrt(1 + 4 * umax**2 / lam))       # rho^2/f^2 at the maximum, fixed kinetic energy
        rec['lam_%.3g' % lam] = dict(m_rho_over_f=mrho, thetadot_over_mrho_max=umax / mrho, t_ok=t_ok,
                                     condensate_fraction_left_at_t_ok=(float(np.interp(t_ok, t, fcond)) if t_ok else None),
                                     rho2_over_f2_at_max=x, mu5_reduction_at_max=math.sqrt(x))
    uv[name] = rec
v_over_Hi = 0.019595917942265423 * MP
uv['vacuum_phase_chi'] = dict(v_over_Hi=v_over_Hi, v_over_f=v_over_Hi / 2000.0,
                              v_over_f_mod_2pi=(v_over_Hi / 2000.0) % (2 * math.pi),
                              note='N2-4 mass phase in the +v vacuum (opposite sign in the -v vacuum); N1: sigma -> 0, phase 0')
uv['thermal_restoration'] = dict(T_max_GeV=d['summary']['T_max_GeV'], f1_GeV=1200 * H_GEV, f2_GeV=2000 * H_GEV,
                                 Tmax_over_f1=d['summary']['T_max_GeV'] / (1200 * H_GEV),
                                 note='lambda-dominated restoration needs T >~ sqrt(3) f')
out['uv'] = uv

# ---- 2. N_s dark radiation: Delta N_eff on the e2e_kinetic history (radiation extension after t_end) ----
fs_turnon = 1.35e-5            # (s_*/M_P = 3e-3), mf_chain G5 verified [Q]
Gs = 1.3e-6                    # H_i, mf_hybrid y_s [Q]; scales ~ H_i
rho = 3 * H**2 * MP**2
fSM = np.array(h['f_SM'])
rs_a3 = fs_turnon * rho[np.searchsorted(t, 1.0)] * a[np.searchsorted(t, 1.0)]**3     # rho_s a^3 (matter-like)
t_end, a_end = t[-1], a[-1]
rhoSM_a4_end = fSM[-1] * rho[-1] * a_end**4
dn = {}
for sname, sfac in [('s_star_3e-3MP', 1.0), ('s_star_1.4e16GeV', (1.4e16 / 2.435e18 / 3e-3)**2)]:
    for Gname, G in [('Gs_1.3e-6', Gs), ('Gs_8.4e-7(H=8.4e10)', Gs * 8.4e10 / 1.3e11)]:
        # rho_DR a^4 = int G rho_s a^3 e^{-G t} a dt in radiation (a = a_end sqrt(t/t_end)), s decays after t_end
        tt = np.logspace(math.log10(t_end), math.log10(t_end) + 9, 20000)
        aa = a_end * np.sqrt(tt / t_end)
        integrand = G * rs_a3 * sfac * np.exp(-G * tt) * aa
        rDR = np.trapezoid(integrand, tt) / rhoSM_a4_end
        dn['%s|%s' % (sname, Gname)] = dict(rhoDR_over_rhoSM=float(rDR), DeltaNeff=float(2.873 * rDR),
                                            t_decay=1 / G, a_decay=float(a_end * math.sqrt(1 / G / t_end)))
out['Ns_dark_radiation'] = dn
out['Ns_dark_radiation_note'] = ('DeltaN_eff = 2.873 rho_DR/rho_SM (g* = 106.75 -> today); f_s from mf_chain at turn-on '
                                 'scaled as s_*^2; the s energy is diluted by nothing before reheating ends (matter era), '
                                 'then grows ~a. Bounds [Q]: Planck 2018 N_eff = 2.99 +- 0.17 (1807.06209), ACT DR6 '
                                 '2.86 +- 0.13 (2503.14454).')

# ---- 3. timing: sphalerons, N lifetimes, BBN ----
alphaW = 1 / 30.0
T_lep = 20.0 * H_GEV
Gsph = 18 * alphaW**5 * T_lep        # [R] Gamma_sph/T ~ 18 alpha_W^5 (Moore) -- order of magnitude
Hlep = 1.66 * math.sqrt(GST) * T_lep**2 / 1.22e19
out['timing'] = dict(T_lepto_GeV=T_lep, sphaleron_rate_over_H_at_T_lepto=Gsph / Hlep,
                     T_sph_eq_GeV_est=18 * alphaW**5 * 1.22e19 / (1.66 * math.sqrt(GST)),
                     note='B-L is conserved by sphalerons; conversion Y_B = (28/79) Y_{B-L} applies once they are '
                          'in equilibrium (T <~ 1e12 GeV), long before BBN',
                     N_lifetimes_s=[1 / (g * H_GEV) * 6.582e-25 for g in [0.533, 2.49, 0.341, 0.123]],
                     s_lifetime_s=1 / (Gs * H_GEV) * 6.582e-25)

# ---- 4. N_s soft-mass variant (not the baseline): late out-of-equilibrium leptogenesis estimate [E] ----
Msoft, mts = 1.0, 1e-5          # H_i, eV (choices inside the mf_higgs-skeptic window 3e-7 .. 0.018 eV)
y2 = mts * 1e-9 * Msoft * H_GEV / 174.1**2
GNs = y2 * Msoft / (8 * math.pi)      # H_i
eps_s_max = 3 / (16 * math.pi) * Msoft * H_GEV * 0.0501e-9 / 174.1**2
r0 = dn['s_star_1.4e16GeV|Gs_1.3e-6']['rhoDR_over_rhoSM'] / 0.886
Tdec = 5.24 * 1868.9 / (a_end * math.sqrt(1 / Gs / t_end))  # bath T at s decay (H_i), radiation extension
YNs = 2 * (3 / 4) * r0 * Tdec / 3.0
out['Ns_soft_variant'] = dict(M_soft_Hi=Msoft, mtilde_s_eV=mts, Gamma_Ns_over_Hi=GNs, eps_s_max=eps_s_max,
                              Y_Ns=YNs, etaB_max=7.04 * 28 / 79 * eps_s_max * YNs,
                              note='needs its own CP phase <~ etaB_obs/etaB_max; breaks N_s -> i g5 N_s softly')
json.dump(out, open(os.path.join(RES, 'checks.json'), 'w'), indent=1)
print(json.dumps(out, indent=1)[:6000])
