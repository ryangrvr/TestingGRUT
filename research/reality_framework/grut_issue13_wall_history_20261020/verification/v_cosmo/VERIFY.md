# VERIFY: module cosmo_observable (verifier key v_cosmo)

21 October 2026. Adversarial verification. I wrote my own code from scratch, using numpy/scipy only, with no
random numbers. I did not import or copy the module's code. I read the module's results JSON only after I had my
own numbers, and only to find where we disagree.

Files in this directory:
- `vc_independent.py` writes `vc_results.json`. Parts (1)-(7): T_dom, R_H, stress tests, LMS units, v_max,
  BBN, and the helper numbers for the conservative argument. Runtime is about 15 s.
- `vc_eos_check.py` writes `vc_eos_results.json`. It computes the period-averaged w of homogeneous phi, X and D
  motions in the unchanged potential.
- `vc_run.log` is the console output.

How to reproduce: `OMP_NUM_THREADS=2 python3 vc_independent.py; python3 vc_eos_check.py`

## Independent method (short)

**SM thermodynamics.**
- Photons plus e+- use exact Fermi-Dirac integrals with mu = 0. Neutrinos decouple instantaneously, with
  (T_nu/T)^3 = g_s,ge(T)/(11/2).
- Two variants: N_eff = 3 (g_rho0 = 3.3626) and N_eff = 3.044, where the neutrino energy density is rescaled
  (g_rho0 = 3.3826).
- A variant with decoupling at 2 MeV changes T_dom by 1e-3.

**Background.**
- a/a0 follows from entropy conservation, and t(ln T) is integrated with DOP853 at rtol 1e-11.
- The walls scale as rho_w = A sigma/t, either inside H (backreaction) or outside it.
- A variant uses rho_w = 2 A sigma H, with H solved self-consistently.
- The integration starts at 5 MeV with t = 1/(2H). The error from this start decays as (a_start/a)^2, which is
  below 1e-6 by the keV range.
- The switch happens where rho_w = f rho_r. After it, rho_w ∝ a^-n.

**Units.** CODATA constants. M_P = 2.435323e18 GeV (reduced), used as H^2 = rho/(3 M_P^2).

## Claims

| # | Claim (module) | My method | My number | Verdict |
|---|---|---|---|---|
| 1a | T_dom = 6.367 keV (A = 0.8, sigma/t, backreaction) | ODE above, g_rho0 = 3.3826 | 6.3670 keV (6.3765 with N_eff = 3); t_dom = 2.7146e4 s; z_dom = 2.711e7 | CONFIRMED |
| 1b | T_dom without backreaction = 5.822 keV; closed form [40 A^2 sigma^2/(pi^2 g M_P^2)]^(1/4) | ODE and closed form | 5.8222 (ODE), 5.8224 (closed form); t_dom = 3.8823e4 s vs Saikawa eq. 2.19, 3.8817e4 s | CONFIRMED |
| 1c | 2 A sigma H convention: 6.924 keV = 2^(1/4) × no-BR; full range 5.45-7.34 keV | ODE | 6.9242; range 5.446 (A = 0.7, no BR) to 7.344 (A = 0.9, 2 A sigma H) | CONFIRMED |
| 1d | g_* at keV is photons plus neutrinos after annihilation | analytic and integral | g_rho = 3.3626 (N_eff = 3) or 3.3826 (3.044); the module uses 3.383. A wrong g = 10.75 would give 4.36 keV | CONFIRMED |
| 2a | R_H (n = 1, frozen) = 1.304e9 vs 70; 1.362e9 / 1.234e9 vs 67 / 74 | ODE to T_CMB; analytic rho_w0 = rho_r(T_dom) T0/T_dom | 1.3041e9; 1.3625e9 / 1.2336e9; H0 = 9.128e10 km/s/Mpc; age 21.42 yr | CONFIRMED |
| 2b | R_H (n = 2 bracket) = 2.505e5 | same | 2.5046e5 (2.617e5 / 2.369e5 vs 67 / 74); age 5.58e4 yr | CONFIRMED |
| 2c | 1 - Omega_w at T_rec = 2.68e-13 (n = 1), 6.67e-9 (n = 2) | same | 2.685e-13, 6.672e-9 | CONFIRMED |
| 2d | R_H range "varying A, switch point and convention": 1.14e9-2.06e9 (n = 1), 2.32e5-2.76e5 (n = 2); "at least 5.4 orders"; "2.3e5-2.1e9 for every variant" | joint grid: A ∈ {0.7, 0.8, 0.9} × {sigma/t, 2 A sigma H} × f_switch ∈ {0.1, 0.3, 1, 3} | **n = 1: 1.036e9-2.296e9; n = 2: 2.167e5-3.436e5**. The minimum gap is 10^5.34. The module varied one factor at a time (its runs list) | CORRECTED (cosmetic; conclusion unchanged) |
| 3a | All wall energy into radiation at T_dom gives Delta N_eff ≈ 7.45 at recombination | rho_dr = rho_r at conversion; ratio constant after e+- annihilation; Delta N_eff = (1 + a N_eff)/a with a = (7/8)(4/11)^(4/3) | 7.447 (N_eff = 3.044); 7.403 (N_eff = 3) | CONFIRMED |
| 3b | Into stable matter gives Omega h^2 ~ 1.1e3 (1134) | rho_r(T_dom)(T0/T_dom)^3/rho_c100, cross-checked by rho/s | 1133.9 (rho/s = 4.13e-6 GeV vs 4.37e-10 GeV for 0.12) | CONFIRMED |
| 4a | 1 GeV^3 = 4.114693e21 J/m^2 = 4.578214e4 kg/m^2 | GeV/(hbar c)^2, then divide by c^2 | 4.114693e21 J/m^2; 4.578214e4 kg/m^2 | CONFIRMED |
| 4b | LMS print 4.22e-9 kg/m^2 (R+M) and 3.85e-9 kg/m^2 (R+M+Lambda), which literally give sigma^(1/3) = 45 keV, inconsistent with their 0.96 / 0.93 MeV by about 1e4 in sigma | fetched arXiv:1505.03673v1 and the published PLB 747 (2015) 426 (SCOAP3 copy); text extraction | Printed values (both versions): "σ < 4.22×10−9 kg/m2 ... 0.96 MeV" and "σ < 3.85 × 10−9 kg/m2 and 0.93 MeV"; Table 1 row "σ < (95%)": 0.96, 0.93 (no unit). Literal reading: 9.218e-14 / 8.409e-14 GeV^3, cube roots **45.17 / 43.81 keV**, which is 21.25 / 21.23 below the labels. sigma is off by 9598 / 9565 | CONFIRMED |
| 4c | Reading kg/cm^2 reconciles them: 0.973 / 0.944 MeV, within 1.5% | conversion | 9.218e-10 / 8.409e-10 GeV^3, which gives **0.9732 / 0.9439 MeV** (1.4% / 1.5% above the labels) | CONFIRMED (as an interpretation, not an erratum) |
| 4d | The kg/cm^2 reading gives G sigma t0 = 4e-6, matching the delta T/T ~ 1e-5 agreement; the literal reading gives 4e-10 | G_SI sigma t0/c, t0 = 13.797 Gyr | 4.09e-6 / 3.73e-6 (literal 4.09e-10 / 3.73e-10). (1 MeV)^3 gives G sigma t0 = 4.4e-6 | CONFIRMED numerically; mild overclaim in wording (see below) |
| 4e | "Residual is a common conversion-constant offset" | none | The common residual is a factor 0.986 in scale (4% in sigma). I found no constant that produces it | NOT CHECKED (speculative statement) |
| 5a | v_max (no domination before today) 47.3-68.1 MeV; 63.32 at A = 0.8 with Planck-like Lambda | root of rho_w(T0) = rho_r0 + rho_m0 + rho_Lambda, scaling with backreaction | 63.317 (A = 0.8); 68.100 (0.643); 66.199 (0.7); 60.879 (0.9); 47.298 (rho_Lambda = 0). Analytic, no BR: sigma_max = 3.04e-5 GeV^3, v = 68.6 MeV | CONFIRMED. The "frozen after acceleration" 61.85 sub-variant is NOT CHECKED |
| 5b | v_max (subdominant CMB, reconciled LMS R+M+Lambda) = 2.07 MeV; literal 96 keV | v = (sigma/0.0942809)^(1/3) | 2.0739 MeV (R+M: 2.138); literal 96.26 keV | CONFIRMED (conditional on the kg/cm^2 reading) |
| 6 | BBN weakly affected: Delta N_eff-equivalent 1.2e-4 at 1 MeV, 0.021 at 0.1 MeV (Omega_w 1.92e-5, 2.93e-3) | rho_w/rho_1nu with T_nu/T from entropy | 1 MeV: Omega_w = 1.909e-5, Delta N_eq = 1.188e-4. 0.1 MeV: 2.928e-3 and 0.02139. A range at 0.1 MeV: 0.0187-0.0241 | CONFIRMED numbers (0.4% difference at 1 MeV). Wording caveat below |
| 7 | Conservative argument: after parity every hidden component has w ≤ 1/3, so rho_h/rho_rad cannot fall below about 1 by recombination | analytic virial identity plus numerical homogeneous-field w | see the section below | SOUND within the action, with stated caveats |

### Overclaims and gaps

1. **The R_H ranges are one-at-a-time ranges** (§0.2 and §4 "every persistent-network variant"). Varying the
   factors jointly gives 1.04e9-2.30e9 (n = 1) and 2.17e5-3.44e5 (n = 2), so the minimum gap is 5.3 orders, not
   5.4. The conclusion is unaffected.
2. **BBN wording.** The quoted Delta N_eff-equivalents at 1 and 0.1 MeV are correct, but rho_w/rho_r grows as
   T^-2. At the deuterium-burning epoch the equivalents are 0.049 (70 keV), 0.100 (50 keV) and 0.28 (30 keV),
   with S - 1 = rho_w/(2 rho_r) = 0.33%, 0.67% and 1.9%.
   - Using Steigman's fit y_D ∝ [6/(eta10 - 6(S-1))]^1.6 (quoted from memory, not re-read), D/H would shift by
     roughly +1-2% [E]. That is comparable to the about 1% observational error on D/H.
   - "Not appreciably" (CHC-4) is plausible, but no BBN network was run. "≤ 0.024 at T ≥ 0.1 MeV" omits the
     epoch that matters for D/H.
3. **LMS agreement statement.** LMS phrase their agreement as an energy scale: "less then 1 Mev (the original
   Zel'dovich bound)". They do not state it in terms of G sigma t0. The module's G sigma t0 argument is its own
   physical inference. The inference is sound: the literal value would put the wall delta T/T near 1e-9, which
   is undetectable.
4. **Published version not cited.** The module cites only v1. The published PLB 747, 426 (2015) prints the same
   SI numbers, so there is no erratum in the journal version (checked here).
5. "Energy scale = sigma^(1/3) in natural units" is the only simple reading that fits to 1.5%. Two alternatives
   fail:
   - the label read as a vev with lambda = 0.01 gives sigma = 8.3e-11 GeV^3, 11× below the reconciled value;
   - reproducing the literal SI value with the label as a vev needs lambda = 1.2e-8, which is implausible.

   Arithmetically, an hbar c unit slip in cm versus m when converting the label gives the same factor 1e4 as the
   kg/cm^2 reading. The paper does not say which of its two numbers is the faithful one. The physical-amplitude
   argument favours the MeV labels.

### Literature checked (read level)

- **Lazanu, Martins, Shellard**, arXiv:1505.03673v1 and the published version (PLB 747 (2015) 426,
  doi:10.1016/j.physletb.2015.06.034, SCOAP3 copy). Read level: section (secs. 2, 4 eqs. 27-28, 6, Table 1, 7;
  text extraction).
  - eq. 27: rho = (1 - lambda)(A eta/V) sigma/(c t).
  - eq. 28: coefficient (6 pi G t0 sigma/c)^2.
  - The published version adds "σ is the surface tension of the domain walls".
  - Matter era: A eta/V = 1.93 with lambda = 2/3, which gives rho t/sigma = 0.643. This matches the module.
  - Subdominance assumption (sec. 2: the walls "always have a small contribution"): confirmed.
- **Hiramatsu, Kawasaki, Saikawa**, arXiv:1309.5001. Read level: section (eqs. 2-4, 29, 30).
  - V = (lambda/4)(phi^2 - eta^2)^2 for a real field.
  - sigma_wall = 2 sqrt(2 lambda) eta^3/3, which matches the BRIEF tension. There is no factor-2 complex-field
    issue.
  - A ≡ rho_wall t/sigma_wall (eq. 4); A ≃ 0.8 ± 0.1 (eq. 29); eps_gw ≃ 0.7 ± 0.4 (eq. 30); H^-1 = 2t (radiation).
  - The module's use is correct.
- Saikawa 1703.02576 and Saikawa-Shirai 1803.01038: **not re-read**. Their numbers were reproduced
  independently: t_dom eq. 2.19, and g_rho0 = 3.383 is consistent with N_eff = 3.044.

## Conservative argument (claim 7): assessment

**Statement.** After parity (rho_h ≈ rho_rad), every hidden component has w ≤ 1/3, so rho_h/rho_rad cannot fall
below about 1 by recombination, because g_* does not change below about 20 keV.

**Logic.** With no non-gravitational hidden→SM energy flow, d ln(rho_h a^4)/d ln a = 1 - 3 w_h ≥ 0 when
p_h ≤ rho_h/3, and conversions inside the hidden sector conserve rho_h. The SM side behaves as follows:
- rho_e±/rho_r = 6.9e-9 at 20 keV, so SM radiation scales as a^-4 to about 1e-8 below 20 keV;
- if parity is reached above 20 keV, e+- annihilation can lower the ratio by at most 1.21, or by 2.60 from
  g = 106.75;
- the argument bounds rho_h/rho_rad, not rho_h/rho_tot: baryons give rho_b/rho_r ≈ 0.5 at z_*.

**w ≤ 1/3 for each component of the unchanged action.**
- **Particle gases** (phi, X, D quanta, any distribution): p ≤ rho/3 by kinetic theory.
- **Thin walls:** w = v^2 - 2/3 ≤ 1/3.
- **Gravitational waves:** w = 1/3.
- **Q-balls and non-relativistic condensates:** w ≈ 0.
- **General field configurations.** For scalar fields, rho - 3p = 4<U> - <phi_i ∂_i U> (virial identity). Take
  the origin at the true vacuum (s = phi - v, X, D). Then the sum of (k - 4) U_k equals
  -lambda v s^2 (2v + s) - 2 g_X^2 v (v + s)|X|^2 - D^2 (m_D0^2 + h_D^2 v (v + s)).
  This is ≤ 0 pointwise wherever phi > 0. So any motion within one vacuum has w ≤ 1/3.
- **Numerical check** (vc_eos_results.json):
  - phi oscillations inside the +v vacuum: w = -0.08 to 0;
  - rotating X at radius 0.1-10 v: w = -0.11 to 0.3331, approaching 1/3 from below;
  - D oscillations: w = 0.006 to 0.331.

**The only w > 1/3 in the action.**
- Homogeneous phi oscillations that cross both vacua with <phi^2> > v^2 satisfy
  3p - rho = lambda v^2 (<phi^2> - v^2) > 0.
- Numerically, w = 0.442 at amplitude 2v and 0.357 at 5v.
- This needs rho ≳ 0.0225 v^4 ≈ 2e10 GeV^4. At parity, rho_tot(T_dom) = 3.7e-21 GeV^4, about 30 orders of
  magnitude less. So it is impossible after parity.
- Running couplings that make U steeper than quartic give δw ~ beta/(18 lambda), negligible, and also only at
  large amplitude.

**No kination (w = 1).**
- There is no flat direction or massless field in the vacuum: m_X^2 = g_X^2 v^2 > 0, and X has no wine-bottle
  potential.
- U(1)_X is global and unbroken in the vacuum, so there is no fixed-radius rotation.

**Loophole (gravitational hidden→SM transfer).**
- Hawking evaporation before recombination needs M ≲ 5.2e9 kg, i.e. closed walls with R ≲ 1 cm. Black-hole
  formation needs R > 4.2 mm (module). Such objects do not exist in a horizon-scale network after parity.
- The graviton-mediated decay of h has an estimated lifetime of about 4e38 s.
- Both are negligible.

**Verdict:** the argument is sound for the unchanged action.
