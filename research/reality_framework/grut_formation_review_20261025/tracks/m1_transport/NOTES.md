# m1_transport: the free-X gas and the parsec-scale transfer forecast

Track of Claude's independent review of Codex's GRUT-M1 chain (issue #13, BRIEF_M1.md), 1 November 2026.
Packets read through git only: 94f660d (grut_m1_synthesis_20261029: REPORT.md, code/collision_rate.py,
code/free_stream.py, free_stream.json, collision_rate.json), f8186fa (grut_m1_stabilized_transport_20261029:
REPORT.md, DERIVATIONS.md, FOLLOWUP_NOTES.md, STABILIZED_CANDIDATE.json, TRANSFER_COMPARISON.json,
code/elastic.py, elastic_s96_n512_a16.json), and the source of the gas: ed73791 (dynamics, 10 % pulse),
bf8ab3b (scattering), c481530 (kinetics), 4fe963f (expansion clock), 6dc9556 (core matching,
SCENARIO_AUDIT.json, late_kin_h0025_96.npz), d4747cb (CANDIDATE_M1.json, PROFILE_M1.csv).

Everything below is a check of one supplied model. Nothing here validates or excludes GRUT-M1; data
comparisons are retrospective.

## 0. Verdict table

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled. "Own" = independent
code in `code/`; "repro" = Codex's checkpoint distribution read from the npz and re-integrated with our own
integrals (labelled reproduction, not independent).

| # | Claim (Codex) | Label | Codex | Ours | Verdict |
|---|---|---|---|---|---|
| 1 | Threshold elastic XX -> XX cross section at lambda_X = 1e-3: sigma = 2.16627561e-7 / m_X^2 | D/C | 2.16627561e-7 | 2.16627561e-7 (own: M = -4 lambda_X + 2 g3^2/m_phi^2 = -0.004 + 0.013333 = +0.009333, sigma = M^2/(128 pi m_X^2)) | CONFIRMED; but it is not a contact cross section: phi exchange is 3.33x the contact amplitude and of opposite sign; contact-only would be 3.98e-8 (M1's lambda_X abs(X)^4) or 2.49e-9 (the lambda abs(X)^4/4 convention of the task text) |
| 2 | Gamma(phi -> X Xbar) = 1.875659e-4 m_X | D | 1.875659e-4 | 1.875659e-4 (own: (2 g m_X)^2 beta/(16 pi m_phi)) | CONFIRMED |
| 3 | m_phi in M1 | D | sqrt(6) m_X (used as m_phi^2 = 6 in the kernels) | m_phi^2 = U''(v) = 2 lambda_phi v^2 = 6 m_X^2 = (0.2264 GeV)^2; the tilt -eps v (phi-v)^3 has zero second derivative at v | CONFIRMED |
| 4 | X Xbar -> X Xbar threshold cross section (not in Codex's transport) | D | - | 2.555e-6 / m_X^2 (s-channel +0.020, t +0.00667, contact -0.004), 11.8x the XX value; pole at sqrt(s) = m_phi needs v = 0.58 c per particle, the 2 MeV gas has v = 0.07 c | NEW; irrelevant for the transport because the checkpoint Xbar fraction is 5e-5 |
| 5 | Gamma_coll/H = 15.72 at 2 MeV (96-bin relativistic kernel; elastic.py NR value 15.877) | C | 15.72 | 16.00 (own: Maxwellian with v_rms = 0.0684, n_X = f Y_Q s(2 MeV), g_*s = 10.75) ; 15.81 (repro: Codex's checkpoint shape, own NR integral) | CONFIRMED |
| 6 | 7.91173 collisions per particle over a factor 10 in a | C | 7.91173 | 7.901 (own, Maxwellian shape frozen) ; 7.811 (repro shape frozen) | CONFIRMED (0.1 %; the 1.3 % repro difference is the frozen non-Maxwellian shape) |
| 7 | Viscosity-proxy crossing Gamma_visc = H at 0.62156 MeV (x = 3.2619) | C | 0.62156 MeV | 0.6227 MeV, x = 3.256 (own) ; 0.6262 MeV (repro) | CONFIRMED; the proxy is 2/3 of the total rate by construction; the total rate itself crosses H at 0.51 MeV (x = 3.98) |
| 8 | Collisionless comoving free-streaming distance to equality: mean 35.51 pc, rms 39.86, median 34.06, p10-p90 12.12-61.36; mean through today 39.61 | C | 35.51 / 39.86 / 34.06 / 12.1-61.4 / 39.61 | own Maxwellian (same v_rms): 36.58 / 39.70 / 35.26 / 17.5-57.3 / 40.81 ; repro (Codex distribution, own background): 35.50 / 39.86 / 34.06 / 12.12-61.35 / 39.61 | CONFIRMED (the Maxwellian differs by 3 % in the mean and more in the percentiles because the checkpoint distribution has kurtosis 1.92 and 2.8 % in its floor bin) |
| 9 | Starting after decoupling: delayed start x = 2, 4, 10 gives 33.69, 31.89, 29.56 pc | C | 33.69 (x=2), 31.89 (x=4) | own Maxwellian from x = 3.26 (0.62 MeV): 33.41 pc; from x = 10: 30.45 pc | CONFIRMED |
| 10 | Background: z_eq = 3420.5 | C | 3420.5 | 3417.0 (Omega_m h^2 = 0.143, T_CMB = 2.7255 K, N_eff = 3.044; Codex uses Omega_m = 0.315, h = 0.674) | CONFIRMED |
| 11 | BGK hierarchy: P_mix/P_cold = 0.997879, 0.979835, 0.893210, 0.849806, 0.764169 at k = 1e4, 3e4, 1e5, 3e5, 1e6 /Mpc with 5.877 % gas | C | as listed | own two-fluid toy [E] (cold + gas fluid, c_s from v_rms, standard background): 0.98-0.99 (3e4, clustered start), 0.88-0.91 (1e5), 0.86-0.89 (3e5), 0.78-0.80 (1e6); asymptotic plateau (1-f)^2 = 0.886 before the growth deficit a^(-6f/5) | NOT CHECKED as a hierarchy (evolve.c / model.py not rerun); the magnitudes are those of a generic cold + warm mixture and the toy agrees to 1-3 % at 1e5 and 1e6 |
| 12 | "A concrete conditional forecast that can now be challenged; not yet an earned, distinctive GRUT prediction" | A | - | k = 1e5/Mpc is M_hm = 5.2e-3 Msun (0.041 Msun in a full-wavelength sphere); the half-suppression k ~ 5.0e4/Mpc is 0.04 Msun; current probes stop at k ~ 25-50/Mpc (M_hm ~ 4e7-3e8 Msun); planned minihalo probes (PTA, caustic microlensing, 21-cm) target 1e-13-1 Msun at O(1) sensitivity to abundance, not 10 % | CONFIRMED (not distinctive); CORRECTED on "testable": a 10.7 % dip at 63 pc is not observable with any current or planned probe; the shape is the standard mixed cold + warm signature (Boyarsky et al. 2009, arXiv:0812.0010) |
| 13 | Lyman-alpha: "no 8-12 % step on the tested range" | C | ratios differ from 1 by ppm at 1-10 h/Mpc | own: absolute(dP/P) <= f k^2 <lambda^2>/3 = 1.4e-9 (10 h/Mpc), 1.3e-8 (30 h/Mpc); fraction-weighted bound: f = 0.059-0.070 is below every Lyman-alpha mixed-DM threshold (F_WDM < 0.35 Bayesian / 0.6 frequentist, arXiv:0812.0010; f < 0.16 even at m = 1 keV, arXiv:2504.06367) and the gas's thermal-equivalent mass is 1.9 MeV versus the 5.3 keV pure-WDM bound (arXiv:1702.01764) | CONFIRMED |
| 14 | Source of the gas: "6.595 % of the original charge in free X", v_rms = 0.0684 c | C/A | 6.595 %, 0.06838 c | the 6.6 % is 0.845 % direct Xbar emission from a supplied 10 % radial pulse (bf8ab3b) multiplied by a 7.37x catalytic chain on the cores (c481530: 1695 Xbar -> 13189 free X) run at the supplied 20 MeV density; v_rms(2 MeV) = v_rms(20 MeV)/10 with v_rms(20 MeV) = 0.66 c set by the pulse energy (2.3 % of the core mass) shared among the chain products; Xbar emission scales as A^2.0-2.1 across the 5-40 % pulses | CORRECTED in interpretation: the gas fraction and its temperature are set by supplied inputs (pulse amplitude, formation clock, reference Q), not by the model |
| 15 | Gas composition includes antiparticles | A/C | checkpoint N_minus = -1.4e-13 (solver residual); late window emits 1.5 Xbar + 3.2 h per core per 1e4/m_X | Xbar fraction of the gas 4.8e-5 at tau = 1000; the gas is a processed antiparticle plume; the late source is not in the transport | CONFIRMED as stated by Codex; negligible for the transport, material for the ledger |
| 16 | Gas-core elastic scattering (omitted by Codex's transport, which has only XX self-scattering) | E | not computed ("gas-core reactions remain open") | own: for X at k = 0.068 m_X the core interior (S^2 + 2 a F^2 > 1 + k^2) is classically forbidden; hard-sphere sigma = 4 pi/k^2 sum (2l+1) sin^2 delta_l with R = 3.6-4.65/m_X (between the barrier edge 3.15/m_X and the half-F radius 4.65/m_X of the Q = 2e5 profile) gives sigma ~ 160-265/m_X^2, Gamma/H ~ 0.7-1.1e6 at 2 MeV (vs 16 for XX), momentum decoupling at x ~ 830-1060, T ~ 3 keV instead of 0.62 MeV | NEW: the decoupling temperature and the k >= 3e5/Mpc part of the curve are not robust until this channel is solved (two-stage fluid toy: 0.89 at 1e5 unchanged; 0.3-0.7 at 3e5 and 0.4-0.6 at 1e6, phase dependent) |

## 1. Cross sections [item 1]

Model (CANDIDATE_M1.json, d4747cb DERIVATION.md, bf8ab3b): with X = v F e^{-i omega t}/sqrt(2), phi = v S, the
dimensionless potential U = (b/4)(S^2-1)^2 - c (S-1)^3 + S^2 F^2/2 + (a/4) F^4 with a = lambda_X/g^2 = 0.1,
b = lambda_phi/g^2 = 3, c = eps/g^2 = 0.2 corresponds to

    V = lambda_phi/4 (phi^2 - v^2)^2 - eps v (phi - v)^3 + g^2 phi^2 |X|^2 + lambda_X |X|^4,

g = 0.1, lambda_X = 1e-3, lambda_phi = 0.03, eps = 0.002, v = 0.924090 GeV, m_X = g v = 0.0924090 GeV.
So M1 (and Candidate 1) normalise the quartic as lambda_X |X|^4, not lambda_X |X|^4/4 as in the task text.
Consequences (own derivation, `code/xsec.py`):

* m_phi^2 = V''(v) = 2 lambda_phi v^2 = 6 m_X^2, m_phi = 2.4495 m_X = 0.22635 GeV. The tilt's second derivative
  vanishes at v, so m_phi is the Candidate 1 value. The phi cubic changes from 6 lambda_phi v to
  6(lambda_phi - eps) v (-6.7 %); it does not enter tree-level XX -> XX.
* phi X X* vertex: 2 g^2 v = 2 g m_X (Codex's `GH = 2 G = 0.2` in m_X units). Gamma(phi -> X Xbar) =
  (2 g m_X)^2 sqrt(1 - 4/m_phi^2)/(16 pi m_phi) = 1.875659e-4 m_X = 1.7333e-5 GeV. Confirms bf8ab3b.
* XX -> XX: contact iM = -i lambda_X x 4 (two X, two X* contractions) -> M_c = -4 lambda_X = -0.004; phi exchange
  in t and u: M = -(2 g m_X)^2/(t - m_phi^2) each, = +4 g^2/6 = +0.0066667 each at threshold (t = u = 0).
  Total M = +0.0093333 (Codex's `AMP = 4 lambda - g_h^2/3 = -0.0093333` is the same modulus). With
  dsigma/dOmega = |M|^2/(64 pi^2 s), s = 4 m_X^2 and the identical-particle 1/2:
  sigma = |M|^2/(128 pi m_X^2) = 2.16627561e-7/m_X^2. CONFIRMED to 10 digits.
* Is the pure contact dominant? No. |M_phi|/|M_c| = 3.33 with opposite sign; sigma_total/sigma_contact = 5.44;
  contact-only sigma would be 3.98e-8/m_X^2 (M1 convention) or 2.49e-9/m_X^2 (lambda|X|^4/4 convention).
  The threshold (constant) amplitude is accurate at the gas velocities: |t|_max/m_phi^2 = 7.8e-4 at
  v_rel = 0.068, changing sigma by 0.1 % at the backward angle.
* X Xbar -> X Xbar at threshold: contact -0.004, t-channel +0.0066667, s-channel -(2 g m_X)^2/(s - m_phi^2) =
  +0.020 (s = 4 below the pole at 6): M = +0.022667, sigma = |M|^2/(16 pi s) = 2.555e-6/m_X^2 (distinct
  particles), 11.8x the XX value. The pole needs sqrt(s) = m_phi, i.e. v = 0.577 c per particle: reached only at
  Codex's 20 MeV clock (gas v_rms ~ 0.66 c, where c481530 includes the resonant redistribution), not at 2 MeV.
  Since the checkpoint gas is 99.995 % X^+, this channel does not affect the 2 MeV transport.
* For orientation, sigma/m = 6.0e-8 cm^2/g: no halo-scale self-interaction.

## 2. Collisions and decoupling [item 2] (own code `code/transport_check.py`)

Number density: n_X = f Y_Q s(T) with f = 13189.13/200000 = 0.065946, Y_Q = 5.0794142e-9, s(2 MeV) with
g_*s = 10.75: n_X = 1.2636e-17 GeV^3 = 1.6013e-14 m_X^3. This equals Codex's inherited
n_core(20 MeV)/1000 x N_+ = 1.6013e-14 m_X^3, so the two clocks are the same clock (constant g = 10.75 from 20 to
2 MeV). Background: own entropy-conserving photon + e+- + decoupled neutrino model, H(2 MeV) = 1.7859e-24 GeV
(Codex 1.7856e-24), H/m_X = 1.93e-23.

Rate: Gamma = n_X sigma <|v_1 - v_2|>, Maxwellian <|v_rel|> = sqrt(2) sqrt(8/3pi) v_rms = 0.08912:
Gamma/H = 16.00 (Codex 15.72 relativistic, 15.88 NR threshold). With Codex's checkpoint shape and our own
isotropic NR relative-speed integral: 15.81. Integrated collisions, Gamma/H ~ x^-4/H(x): 7.90 (frozen Maxwellian)
over x = 1-10; Codex evolves the shape and gets 7.91. Viscosity proxy (2/3 of the total rate by Codex's
definition, not a solved eigenvalue) crosses H at x = 3.256, T = 0.6227 MeV (Codex 3.262, 0.62156). The
total rate crosses H at 0.51 MeV. All CONFIRMED.

Caveat on the decoupling claim (NEW, see item 16 and section 4): the XX self-scattering rate is not the gas's
dominant momentum-transfer channel. Scattering of the gas X on the cold cores has a cross section of order the
core's geometric size, 1e2-1e4/m_X^2 versus 2e-7/m_X^2, and even with n_core = n_X/13189 that is Gamma/H ~
1e6 at 2 MeV.

## 3. Free streaming, scales and observability [item 3]

Own background (Omega_m h^2 = 0.143, h = 0.674, Omega_gamma h^2 = 2.474e-5, N_eff = 3.044, flat Lambda):
a(2 MeV) = 8.3945e-11, a_eq = 2.9257e-4 (z_eq = 3417). The comoving distance per unit velocity at 2 MeV,
I = c a_i int dln a/(a^2 H), is 580.5 pc to equality and 647.5 pc to today; the comoving Hubble radius at 2 MeV
is 42.7 pc. A Maxwellian with v_rms = 0.0684 (mean v = 0.0630) gives mean 36.6 pc, rms 39.7 pc, p10-p90
17.5-57.3 pc to equality, 40.8 pc to today. Codex's distribution integrated with exact relativistic v(a) in our
background reproduces its numbers to 4 digits (35.50 / 39.86 / 34.06 / 12.12-61.35 / 39.61 pc). Starting after
the 0.62 MeV crossing instead of 2 MeV only lowers the mean to 33.4 pc: the RD contribution is logarithmic in a.

Suppression scales. 1/lambda_mean = 2.7e4/Mpc; 2 pi/lambda_mean = 1.7e5/Mpc; the "k_fs = aH/v_rms at
2 MeV" convention gives 3.4e5/Mpc; a thermal relic with the same v_rms today (m = 1.9 MeV) has a Viel-fit
half-mode k_1/2 = 2.76e4/Mpc. Codex's BGK curve reaches half of its (1-f)^2 plateau at k = 5.0e4/Mpc, i.e.
1.8x the thermal-equivalent k_1/2. These are mutually consistent conventions for the same ~30-60 pc scale.

Halo masses (M = (4 pi/3) rho_m (pi/k)^3, half-mode convention of Schneider et al. 2012, arXiv:1112.0330;
rho_m = 3.97e10 Msun/Mpc^3; the sphere of a full wavelength is 8x larger):

| k [/Mpc] | 2 pi/k | M_hm [Msun] | M(2pi/k sphere) | Codex P_mix/P_cold |
|---|---|---|---|---|
| 1e4 | 628 pc | 5.2 | 41 | 0.9979 |
| 3e4 | 209 pc | 0.19 | 1.5 | 0.9798 |
| 5.0e4 (half suppression) | 126 pc | 0.041 | 0.33 | 0.943 |
| 1e5 | 62.8 pc | 5.2e-3 | 0.041 | 0.8932 |
| 3e5 | 20.9 pc | 1.9e-4 | 1.5e-3 | 0.8498 |
| 1e6 | 6.3 pc | 5.2e-6 | 4.1e-5 | 0.7642 |

Mass inside the mean free-streaming sphere: 8e-3 Msun.

What is probed today (thermal-WDM bounds converted with the Viel et al. 2005 fit, astro-ph/0501562):
* Lyman-alpha forest: k <~ 10-30 h/Mpc; m_WDM > 5.3 keV (3.5 keV with a free thermal history), Irsic et al. 2017,
  PRD 96, 023522, arXiv:1702.01764: k_1/2 = 40 (25)/Mpc, M_hm = 8e7 (3e8) Msun.
* Milky Way satellites: m_WDM > 6.5 keV, Nadler et al. 2021, PRL 126, 091101, arXiv:2008.00022: k_1/2 = 51/Mpc,
  M_hm = 4e7 Msun.
* Strong-lensing flux ratios: M_hm < 10^7.8 Msun, Gilman et al. 2020, MNRAS 491, 6077, arXiv:1908.06983.
* Stellar streams: subhalos of 1e6-1e9 Msun, m_WDM > 3.6 keV, Banik et al. 2021, JCAP 10 (2021) 043,
  arXiv:1911.02663.
Proposed probes of minihalos (1e-3-1 Msun is k = 1.7e5-1.7e4/Mpc in the half-mode convention):
* Pulsar timing: Dror, Ramani, Trickle, Zurek 2019, arXiv:1901.04490; Ramani, Trickle, Zurek 2020, JCAP 12
  (2020) 033, arXiv:2005.03030 (SKA-like PTA sensitive to subhalos down to ~1e-13 Msun at CDM concentration,
  as an O(1) constraint on their abundance); Lee, Mitridate, Trickle, Zurek 2021, arXiv:2012.09857.
* Caustic-crossing microlensing of highly magnified stars: Dai & Miralda-Escude 2020, AJ 159, 49,
  arXiv:1908.01773 (minihalo detection through flux fluctuations; abundance-level sensitivity).
* 21-cm: Munoz, Dvorkin, Cyr-Racine 2020, PRD 101, 063526, arXiv:1911.11144: k = 3-80/Mpc at 1-30 %
  precision (HERA), m_WDM up to 8-14 keV; three decades in k short of 1e5/Mpc.
* WIMP-CDM reference cutoffs, M_cut ~ 1e-11-1e-3 Msun: Green, Hofmann & Schwarz 2005, JCAP 08 (2005) 003,
  astro-ph/0503387; Loeb & Zaldarriaga 2005, PRD 71, 103520, astro-ph/0504112; Bringmann 2009, New J. Phys.
  11, 105027, arXiv:0903.0189 [R, numbers from memory].

Verdict on "distinctive, testable": the forecast's dip (0.2 % at 5 Msun, 2 % at 0.2 Msun, 10.7 % at 5e-3 Msun)
sits 3.5 decades in k and 10 decades in mass below the smallest scale any current probe reaches. The proposed
minihalo probes cover the right mass window but at order-unity sensitivity to the minihalo abundance, with the
abundance itself controlled by the primordial spectrum extrapolated five decades beyond the CMB/Lyman-alpha
range, by the survival of 1e-3 Msun minihalos to z = 0, and, in this model, by the fact that the cold cores (17 TeV
objects with no SM portal) have no WIMP-like kinetic-decoupling cutoff of their own. A gradual 10 % linear
suppression at 63 pc is not observable with any current or planned probe. It is also not distinctive: for a fixed
gas fraction f and velocity the curve is the standard cold + warm mixture (Boyarsky, Lesgourgues, Ruchayskiy &
Viel 2009, JCAP 05 (2009) 012, arXiv:0812.0010; also Anderhalden et al. 2013, arXiv:1212.2967 and Diamanti
et al. 2017, arXiv:1701.03128 [R]): plateau (1-f)^2 = 0.886, deepened by the growth deficit a^(-6f/5) of the cold
part. Our two-fluid toy (`code/mixture.py`, [E]) reproduces Codex's magnitudes with no hierarchy at all: 0.88-0.91
at 1e5, 0.78-0.80 at 1e6.

Lyman-alpha with a 6 % warm fraction at v = 0.07 c (2 MeV): the no-gravity bound |dP/P| <= f k^2 <lambda^2>/3
is 1.4e-9 at 10 h/Mpc and 1.3e-8 at 30 h/Mpc. The fraction-weighted statement is stronger than needed: at
f = 0.059 (all matter) or 0.070 (dark matter) the mixed-model Lyman-alpha analyses allow any warm mass
(arXiv:0812.0010: F_WDM < 0.35 at 95 % Bayesian, < 0.6 at 99.7 % frequentist; Garcia-Gallego et al. 2025,
arXiv:2504.06367: f_WDM < 0.16 at 2 sigma even for m = 1 keV), and the gas is in any case colder than a 1.9 MeV
thermal relic, 360x the 5.3 keV pure-WDM bound. Codex's "no 8-12 % step at Lyman-alpha scales" is CONFIRMED; no
present constraint applies.

## 4. The source of the gas and what would change the forecast [item 4] (`code/source_scalings.py`)

The 6.6 % free X is not a formation-ensemble result. Its chain, from Codex's own ledgers:
1. ed73791: a prepared Q = 200000 reference core is given a supplied 10 % radial Gaussian pulse; the open
   spherical receiver records net negative charge 0.00799 Q (5 %: 0.0020, 20 %: 0.0304, 40 %: 0.0949).
2. bf8ab3b: the spectral split gives direct Xbar/Q_0 = 0.002088, 0.008448, 0.03447, 0.1459 for the four pulses,
   i.e. Xbar emission ~ A^2.0-2.1; nothing in M1 fixes A.
3. c481530: at the supplied 20 MeV density (tau = 1000, H dt = 0.0016) the 1695 Xbar + 95 X per core become
   11629 free X + 0.56 Xbar by charge exchange on cores plus a neutral-production chain (X + core -> h + core,
   h -> X Xbar, Xbar + core -> X + core(Q-2)); through the extended source window N_+ = 13189 (6dc9556). The
   chain multiplies the free charge 7.37x relative to one exchange per Xbar; eta_equivalent = 0.934, p = 0,
   frozen-mass test 0.0668(1-eta) = 0.0044 < 0.01.
4. v_rms(2 MeV) = 0.0684 c is v_rms(20 MeV)/10, and v_rms(20 MeV) = 0.66 c (2526 m_X of gas kinetic energy for
   11629 particles) is the pulse's 2.3 % energy excess shared among the chain products.
5. The checkpoint gas is 99.995 % X^+ (Xbar fraction 4.8e-5); Codex's statement that the composition includes
   antiparticles refers to the emitted plume (mostly Xbar) and to the continuing late source (1.5 Xbar + 3.2 h per
   core per 1e4/m_X in the 10000-20000 window, 103 m_X of core excess still unreleased), which the transport omits.

What changes the forecast:
* Gas fraction f: linear in the Xbar emission for a fixed chain factor; 5 % pulse -> f ~ 1.6 % and a 3 % dip;
  20 % pulse -> f >~ 27 % (the fixed-chain extrapolation exceeds the mass budget at 40 %, where Codex's own
  completion fails). Alternatively f = 1 - eta up to the 15 % allowed by the frozen mass test at p = 0 gives up to
  a ~28 % dip. The suppression amplitude is therefore an input (pulse amplitude, or eta), not an M1 output.
* Formation clock T_f: v_rms at 2 MeV scales as T_f/20 MeV, so lambda_fs and M_hm scale as T_f and T_f^3:
  T_f = 200 MeV -> 366 pc, k_1/2 ~ 2.7e3/Mpc, M_hm ~ 40 Msun; T_f = 2 MeV -> 3.7 pc, 4e-5 Msun. Codex states the
  20 MeV clock is an explicit input (H/m_X = 1.9e-21 there, not the handoff's 3e-14).
* Velocity distribution: a fragmentation ensemble (Kusenko & Shaposhnikov; Enqvist & McDonald; Kasuya &
  Kawasaki) emits a spectrum of Q-balls and free quanta set by the condensate instability, not by one pulse on one
  reference core; the chain factor 7.4 depends on the density clock through tau.
* Gas-core elastic scattering (NEW, [E], `code/core_barrier.py`, `code/mixture.py`): at k = 0.068 m_X the
  interior of the Q = 2e5 core (central F = 2.63, S = 0.126; X "barrier" S^2 + 2 a F^2 - 1 = +0.38 m_X^2 at the
  centre, barrier edge at r = 3.15/m_X, half-F radius 4.65/m_X, F below 1 % of F(0) beyond 11.5/m_X) is
  classically forbidden to the X, so the core scatters nearly as a hard sphere: sigma ~ 4 pi/k^2 sum (2l+1)
  sin^2 delta_l ~ 160-265/m_X^2 for R = 3.6-4.65/m_X (1100-4000/m_X^2 if the 10-20/m_X field tail counts).
  With n_core(2 MeV) = 1.21e-18 m_X^3 the rate is Gamma/H ~ 0.7-1.1e6 (XX: 16), falling as a^-2, so the gas
  stays momentum-locked to the cold cores until x ~ 830-1060, T ~ 3 keV, not 0.62 MeV. Codex's bf8ab3b packets show 99 % elastic reflection for Xbar at
  k = 0.4-1 with ~1 % inelastic channels; the X case at k ~ 0.07 is not computed anywhere. Consequences in the
  two-stage fluid toy (tight coupling as one fluid with c_eff^2 = f c_s^2 until Gamma = H, then two fluids):
  k = 1e5 stays at 0.887-0.904 for sigma = 0.01-1000/m_X^2 (the plateau is set by f), but k = 3e5 moves to
  0.32-0.73 and k = 1e6 to 0.006-0.60 depending on the acoustic phase at decoupling. After decoupling the gas
  still free-streams 16-26 pc (logarithmic RD integral). The sub-30 pc part of Codex's curve and its decoupling
  temperature are therefore not robust within the model; the 60-600 pc part is set by f alone.

Could the parsec-scale curve ever become distinctive relative to generic mixed cold + warm dark matter? Not in
its shape: any model with a cold fraction 1-f and a warm fraction f with the same v_rms gives the same linear
curve, and a mixed-DM fit (f, k_1/2) absorbs both of M1's free inputs. The only M1-specific ingredients are
the non-Maxwellian shape (kurtosis 1.92 -> 1.67, invisible in P(k)) and the gas-core coupling, which would
replace the free-streaming cutoff by a dark-acoustic-like feature at k >~ 3e5/Mpc (M <~ 2e-4 Msun), again at
unobservable scales. The forecast becomes a prediction only if f and T_f are derived from a formation ensemble;
even then it would be a prediction of a 1e-3-1 Msun minihalo deficit at the 10 % level, below every projected
probe's sensitivity.

## 5. What this track did not check

* Codex's BGK hierarchy (evolve.c, model.py, linear_*.json), its collision-operator closure and its
  initial-condition matching were not rerun; item 11 is a magnitude check with a different (fluid) method.
* The CLASS runs, the holder and the charge-power numbers are outside this track (m1_holder, m1_vacuum).
* The 64- versus 96-bin source dependence and the C-code conservation checks were only read, not reproduced.
* The X-core elastic cross section is an estimate; the three-channel radial problem at k ~ 0.07 m_X with the
  anomalous and neutral couplings (Codex's own bf8ab3b machinery) would settle it.
* No statement here validates or excludes M1; all data comparisons are retrospective and the forecast is
  conditional on supplied inputs.

## 6. Files and reproduction

All code is own code unless marked. Run from this directory with one thread:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    python3 -I code/xsec.py                # item 1: results/xsec.json (own)
    python3 -I code/transport_check.py     # items 2-3: results/transport_check.json (own; the
                                           #   'checkpoint_reproduction' blocks read Codex's late_kin_h0025_96.npz,
                                           #   extracted with git show 6dc9556:... into $SP/review/m1_transport/codex/)
    python3 -I code/mixture.py             # item 3-4 toy: results/mixture_estimate.json (own, [E]; ~10 min)
    python3 -I code/source_scalings.py     # item 4 ledger arithmetic: results/source_scalings.json
    python3 -I code/core_barrier.py        # item 4 X-core estimate: results/core_barrier.json (reads d4747cb PROFILE_M1.csv)

`results/claims.json` holds the verdict table in machine-readable form. Scratch and Codex data copies live in
/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/m1_transport/.
