# m1_transport: skeptic's pass

Adversarial check of `tracks/m1_transport/NOTES.md` (Claude's review of Codex's GRUT-M1 gas and transfer
packets, issue #13, BRIEF_M1.md), 1 November 2026. Codex packets read through git only (f8186fa, 94f660d,
d4747cb, bf8ab3b, c481530, 4fe963f, 6dc9556, ed73791). All code here is my own unless marked "repro";
it lives in `skeptic/code/`, results in `skeptic/results/`. One thread, foreground, no commits, no 3D runs.
Nothing here validates or excludes GRUT-M1; data comparisons remain retrospective.

Method: every key number of the track was recomputed by a different route (real-component Feynman rules
instead of the charge basis; a separately written background and quadrature; a coupled three-channel radial
scattering solver instead of a hard-sphere formula; an independently written two-fluid integrator and two
further ODE solvers for the toy). `refuted = true` only where a recomputation contradicts the reviewer.

## Verdict table (reviewer | skeptic)

| # | Claim | Reviewer | Skeptic (own) | Refuted? |
|---|---|---|---|---|
| T1 | sigma(XX->XX) threshold | 2.16627561e-7/m_X^2; M = -4 lambda_X + 2 g3^2/m_phi^2 = +0.009333 | 2.16627561e-7/m_X^2 from the (phi1, phi2) basis rotated to charge states: contact -0.004, t+u exchange +0.013333; finite-v correction -0.23 % at v_rel = 0.068 | No. But "correction 1" (not a contact cross section) corrects the task text, not Codex: f8186fa DERIVATIONS.md and c481530 DERIVATION.md write the amplitude as 4 lambda_X - g_h^2/3, contact plus exchange |
| T2 | Gamma(h->X Xbar), m_phi | 1.875659e-4 m_X; sqrt6 m_X | 1.875659e-4 m_X = 1.7333e-5 GeV (sum of two identical-real-daughter widths); sqrt6 m_X, tilt second derivative zero at v | No |
| T3 | sigma(X Xbar) threshold | 2.555e-6/m_X^2 | 2.5553e-6/m_X^2 (s +0.020, t +0.00667, contact -0.004) | No |
| T4 | Gamma_coll/H at 2 MeV | 16.00 (Maxwell) / 15.81 (Codex shape) | 15.99 / 15.81; H(2 MeV) = 1.78593e-24 GeV; n_X = 1.6013e-14 m_X^3 | No |
| T5 | collisions over x = 1-10 | 7.90 | 7.898 (frozen Maxwellian) | No |
| T6 | viscosity-proxy crossing | 0.6227 MeV, x = 3.256 | 0.6228 MeV, x = 3.255; total rate at x = 3.98 | No |
| T7 | free streaming to equality | Maxwell 36.58/39.70 pc, today 40.81; Codex shape 35.50/39.86; z_eq 3417 | 36.57/39.69, 40.79; 35.50/39.86 (repro); z_eq 3417.0; I_eq = 580.5 pc per unit v; Hubble radius 42.7 pc | No |
| T8 | BGK hierarchy values | NOT CHECKED; fluid toy 0.88-0.91 (1e5), 0.78-0.80 (1e6) | own two-fluid code: 0.8838/0.9116 at 1e5, 0.7985/0.7954 at 1e6 (erased/clustered), 0.8857/0.9905 at 3e4, 0.8691/0.8881 at 3e5; plateau (1-f)^2 = 0.8859 | No (toy numbers reproduced to 4 digits) |
| T9 | "testable, not distinctive" -> CORRECTED | M_hm(1e5) = 5.2e-3 Msun; probes stop at k ~ 25-50/Mpc; unobservable | M_hm(1e5) = 5.15e-3 Msun; thermal-equivalent 1.90 MeV; k_hm(T = 1/2) = 27636/Mpc; literature values (Irsic 2017 5.3 keV, Nadler 2021 6.5 keV, Gilman 2020 10^7.8 Msun, Banik 2021 3.6 keV, Boyarsky 2009 F < 0.35/0.6, Viel 2005 fit) check out from memory | No, but the verdict is OVERSTATED: Codex wrote "can now be challenged" and disclaimed observational status; "testable" is the reviewer's paraphrase, so this is NEW, not CORRECTED |
| T10 | Lyman-alpha toy bound | 1.4e-9 (10 h/Mpc), 1.3e-8 (30 h/Mpc) | 1.40e-9, 1.26e-8 | No |
| T11 | source of the gas; scalings | 0.845 % x 7.37; f ~ A^2; **lambda_fs ~ T_f: 200 MeV -> 366 pc, 2 MeV -> 3.7 pc** | 0.845 % x 7.37 and f ~ A^2.0-2.1 confirmed; **clock scaling inverted**: p(2 MeV) = p_f a_f/a(2 MeV) ~ p_f (2 MeV/T_f), so lambda_fs ~ 1/T_f (times a log): 394 pc (T_f = 2 MeV, v(2 MeV) = 0.56), 46 pc (20 MeV, including 6 pc before 2 MeV), 4.6 pc (200 MeV), 0.33 pc (2 GeV) | **Yes** (the T_f scaling and its numbers); the rest stands. Also 0.66 is p/m at 20 MeV, v = 0.57 |
| T12 | antiparticle content | 4.8e-5 of the gas | 0.559/11629 = 4.8e-5 | No |
| T13 | gas-core elastic scattering (NEW) | interior "classically forbidden"; hard sphere R = 3.6-4.65/m_X, sigma = 160-265/m_X^2; Gamma/H = 0.7-1.1e6; x_dec 830-1060 (T ~ 3 keV); two-stage toy 0.32-0.73 at 3e5, 0.006-0.60 at 1e6 | coupled three-channel solve at k = 0.0684: delta_0 = -0.356, delta_1 = -0.027, delta_2 = +0.003; sigma_el = 332/m_X^2 (a_s = 5.43/m_X); Gamma/H = 1.43e6; x_dec = 1194 (T ~ 2.3 keV); k-scan 0.02-0.1 gives 332-342, a resonance 1333 at k = 0.15; uncoupled diagonal potential would give 1561/m_X^2 (a_s = 17). Two-stage toy numbers reproduced exactly by independent code and three integrators | No on the magnitude (confirmed, 25-100 % larger); the MECHANISM is wrong: the interior is propagating, see below |

## What is refuted

**T11, formation-clock scaling.** The reviewer writes "v_rms at 2 MeV scales as T_f/20 MeV, so lambda_fs and
M_hm scale as T_f and T_f^3" and tabulates 366 pc / 40 Msun for T_f = 200 MeV and 3.7 pc / 4e-5 Msun for
T_f = 2 MeV (`results/source_scalings.json` `clock_scaling`, NOTES section 4). Momentum redshifts as
p ∝ 1/a: a gas released with p_f ≈ 0.68 m_X (set by the pulse energy, reviewer's own step 4) at a higher
temperature redshifts longer before 2 MeV, so p(2 MeV) ≈ p_f (2 MeV/T_f) and lambda_fs ∝ 1/T_f up to a log.
Own integration (exact relativistic v(a), g_* table above 2 MeV, `sk_rates.py`): lambda_rms to equality =
394 pc (2 MeV; the gas is then semi-relativistic, v = 0.56), 46 pc (20 MeV), 4.6 pc (200 MeV), 0.33 pc
(2 GeV). The qualitative conclusion ("the scale is an input, not an M1 output") survives; the direction and
the numbers do not. An earlier clock makes the dip smaller-scale and even less observable.

## What is corrected without changing the number

**T13, mechanism.** The reviewer's premise is that at k = 0.068 m_X the core interior is classically
forbidden because the charged diagonal S^2 + 2aF^2 - 1 = +0.40 at the centre, hence a hard sphere. In the
canonical three-channel operator (bf8ab3b DERIVATION.md, used by Codex for the antiparticle) the anomalous
coupling aF^2 = 0.69 at the centre is larger than that diagonal excess. The charged 2x2 block
[[0.3955, 0.6922], [0.6922, 0.8212]] has a negative eigenvalue (-0.12) at r = 0, and the uniform-condensate
dispersion (k^2 + aF^2 - nu^2)^2 = 4 nu^2 Omega^2 + a^2 F^4 is gapless: a slow X propagates inside the core
as a Bogoliubov phonon with k_in = 0.22 m_X (sound speed^2 0.31). The interior is a refracting medium, not a
wall. Solving the coupled radial problem (own solver, `sk_xcore.py`, validated on the vacuum, an analytic
s-wave square barrier and the hard-sphere limit; closed-channel leakage 2e-9; step and r_max converged)
gives nevertheless a geometric answer: sigma_el = 332/m_X^2 at k = 0.068 (scattering length 5.4/m_X,
effectively repulsive), 332-342 over k = 0.02-0.1, with a resonance (1333/m_X^2) near k = 0.15. So the
reviewer's Gamma/H ~ 1e6 and keV-scale momentum decoupling are confirmed and slightly strengthened
(Gamma/H = 1.4e6, x_dec ≈ 1200, T ≈ 2.3 keV), but by a different physics than stated. Taking the reviewer's
diagonal potential literally, including its attractive S^2 - 1 < 0 shell at r = 3.3-12/m_X, would have given
1561/m_X^2 (a near-threshold state, a_s = 17/m_X); the hard-sphere shortcut lands within a factor 1.3-2 of
the coupled result by accident. Codex's own antiparticle numbers already implied this size: at k = 1 the
charge-exchange cross section 2.76/m_X^2 is 1.30 % of events, i.e. sigma_tot(Xbar) ≈ 210/m_X^2.

**T13, two-stage toy.** The large suppressions (0.32-0.73 at k = 3e5/Mpc, 0.006-0.60 at 1e6) are not a
solver artefact: an independently written integrator and Radau/DOP853 at rtol 1e-10 reproduce them to
five digits, and run_two_stage with a_dec -> a_i reproduces the no-drag "clustered" run (ratio 1.0000).
Tracing the trajectory (`sk_trace_two_stage.py`) shows the mechanism: while gas and cores are tightly
coupled the mixture is a fluid with c_eff^2 = f c_s^2; in radiation domination the sub-horizon equation has
zero Hubble friction and only logarithmic (Meszaros) growth, so even (k c_eff/aH)^2 = 0.015 (k = 3e5,
0.12 k_J,eff) bends delta from 3.0 to 2.4 over the ~7.6 coupled e-folds, and for k = 1e6 a short coupled
phase reverses the Meszaros velocity so delta_c coasts negative. This is a dark-acoustic-like effect of the
gas pressure acting on the whole mixture, and it is phase-sensitive as the reviewer says. It strengthens the
"not robust for k >= 3e5/Mpc" conclusion within the toy, with the toy's limits: isothermal c_s, no
radiation perturbations, instantaneous decoupling, fluid gas.

## Overstated or misdirected verdicts

* **T9 CORRECTED -> should be CONFIRMED + NEW.** Codex's sentence is "a concrete conditional forecast that can
  now be challenged. It is not yet an earned, distinctive GRUT prediction", preceded by "not nonlinear
  structure predictions or an observational M1 fit". "Testable" is the reviewer's word. The observability
  analysis (5e-3 Msun, 3.5 decades below current probes, minihalo probes at O(1) sensitivity) is correct as
  far as I can check and is useful, but it corrects nobody. Digest correction 3 should be reworded.
* **T1 correction 1 is addressed to the task text.** Codex's packets state the amplitude as contact plus
  exchange and never use the lambda|X|^4/4 convention. The reviewer's NOTES say so; the digest does not.
* **Convention label.** `results/transport_check.json` describes the Viel/Schneider half-mode as "T^2 = 1/2"
  while the code uses 2^(nu/5), i.e. T = 1/2 (the standard Schneider et al. 2012 half-mode; T^2 = 1/2 would
  give 19926/Mpc instead of 27628/Mpc). The number is right; the label is wrong. The comparison with the BGK
  half-suppression point is nonetheless apples-to-apples: at P_mix/P_cold = (1 + (1-f)^2)/2 the warm transfer
  function is 0.51, so the 1.8x ratio compares two T_w = 1/2 points.
* **T11 wording.** "v_rms(20 MeV) = 0.66 c" is p_rms/m_X; the speed is 0.57 c. Harmless (only p redshifts).

## Checks that passed without remark

T1-T7, T10, T12 recomputed from scratch agree with the reviewer to all quoted digits (table above). The
reviewer's literature values used for T9/T10 (Irsic et al. 2017; Nadler et al. 2021; Gilman et al. 2020;
Banik et al. 2021; Boyarsky et al. 2009; Viel et al. 2005; Munoz et al. 2020; Dror/Ramani/Lee et al. PTA
papers; Dai & Miralda-Escude 2020) agree with what I recall; I could not fetch arXiv:2504.06367 and mark it
[R]. Units: sigma in m_X^-2, Gamma = n sigma v with n in m_X^3 and the trailing m_X in GeV, Y_Q as charge per
entropy, f as charge fraction for n_X and 5.877 % mass fraction of all matter in the toy, consistent with
Codex's FOLLOWUP_NOTES distinction. Signs: the h exchange is attractive and the contact repulsive, opposite
signs as stated; the net threshold amplitude is attractive, Born-valid (alpha = 8e-4, alpha/v = 0.012, no
Yukawa bound state: alpha m_X/m_phi = 3e-4 << 0.84).

## New issues for the record

1. T11 clock scaling inverted (above); `source_scalings.json` `clock_scaling` and NOTES section 4 need the
   1/T_f direction.
2. T13 mechanism: Bogoliubov phonon interior, geometric cross section anyway; resonance at k = 0.15.
3. The X-core elastic size was implicit in Codex's own antiparticle numbers (sigma_tot(Xbar, k = 1) ~
   210/m_X^2); the omission in f8186fa could have been quantified from bf8ab3b without a new solver.
4. With Gamma_Xcore/H ~ 1e6 the gas is diffusive at 2 MeV (mean free path 3e15/m_X versus a 62.8 pc comoving
   wavelength of 8e22/m_X physical), so Codex's statement that the gas "free streams" between 2 MeV and
   0.62 MeV and the whole BGK matching at 2 MeV presuppose a channel the packets did not evaluate; the
   reviewer's point, sharpened.
5. Nothing in this pass touches Codex's hierarchy code, CLASS runs or the 64/96-bin source dependence; those
   remain uncertified, as the reviewer says.

## Files

* `skeptic/code/sk_xsec.py` -> `skeptic/results/sk_xsec.json` (T1-T3, real-component basis)
* `skeptic/code/sk_rates.py` -> `skeptic/results/sk_rates.json` (T4-T7, T9-T11; own background; clock scaling)
* `skeptic/code/sk_xcore.py` -> `skeptic/results/sk_xcore.json` (T13 coupled three-channel X-core scattering; ~2 min)
* `skeptic/code/sk_mixture_check.py` -> `skeptic/results/sk_mixture_check.json` (T8/T13 toy: repro of the reviewer's
  mixture.py plus own two-fluid integrator; ~8 min)
* `skeptic/code/sk_trace_two_stage.py` (trajectory trace behind the two-stage interpretation)
* `skeptic/results/sk_claims.json` (this table in machine-readable form)

Reproduce from `tracks/m1_transport/skeptic/` with `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SP=<scratchpad>` after `git show 6dc9556:.../late_kin_h0025_96.npz` and `git show d4747cb:.../PROFILE_M1.csv`
into `$SP/review/m1_transport_skeptic/codex/`.
