# Second-agent review: Blind second implementation of the floor and six-maser envelope

Separate Claude agent spawned by the main session for issue #13; it returned structured findings, rendered here verbatim. Scripts and outputs are in this directory; paths inside them refer to the session scratchpad. Reviewer: blind_envelope (Claude Code second-agent retrospective review; I did not read Codex's audit.py, report.py, DERIVATION.md, COMPARISON.json, DISTANCE_ENVELOPE.csv, maser_table1.json or any Claude script. I froze my outputs at 2026-10-20T20:44:41Z (out/BLIND_FREEZE_sha256.txt) and only then read Codex's REPORT.md at 2ee7f82 and issue #13 to compare.)

## Findings

### 1. [PASS] Charge energy floor m* = v min(gX,(lambda_phi lambda_X)^(1/4)) with rho >= m*|n_Q| (canonical, -+++)

**Deciding evidence.** Derived independently. Minimising A/x^2 + lX x^2 over x = |X| gives 2 sqrt(A lX). This leaves f(s) = gX^2 s + sqrt(lp lX)|v^2 - s|, which is piecewise linear in s = phi^2, so min f = v^2 min(gX^2, sqrt(lp lX)). With n_Q = -2 Im(X* Xdot), the identity |Xdot|^2 + m^2|X|^2 - m|n_Q| = |Xdot + i sgn(n_Q) m X|^2 holds numerically to a maximum relative error of 6.7e-15. Benchmark: m* = 56.23413252 GeV, mu = 0.562341325. A direct 2D/3D Nelder-Mead minimisation of U/|X|^2 matches m* to a relative error of 1e-16, and a brute-force grid matches it to 8.6e-6. Pointwise sampling found 0 violations in 2e6 generic samples (minimum rho/(m*|nQ|) = 19.2) and 0 in 1e6 adversarial samples near saturation (minimum 1.00000000000004). A minimiser of rho/|nQ| returns 1.0000000000000. The bound is therefore sharp: it saturates at phi = 0, |X| = (lp/(4 lX))^(1/4) v = 1.2574 v, Xdot = -i sgn m* X. Six random coupling sets in both classes also agree to <= 2e-16.

**Scope or correction.** The proof additionally needs U_D >= 0, i.e. mD0^2 >= 0 and lambda_D >= 0. The benchmark's mD0 = gX v satisfies this. A tachyonic mD0^2 < 0 with m_D^2 > 0 (example: mD0^2 = -25 GeV^2) makes U_D negative at phi = 0 (-1.56e5 GeV^4) and voids U >= m*^2|X|^2. The floor is classical; one-loop terms are not covered.

### 2. [PASS] omega_h,min = (m* r/m_b) omega_b and Hfloor = 100 sqrt(omega_b + omega_h,min): 116.66353672 nominal; 113.25181752 with omega_b - 3 sigma and 1.01 m_b

**Deciding evidence.** I confirmed omega_b = 0.02233 +- 0.00036 myself from the Mossa et al. accepted manuscript R2 (Padua repository, dated 2 Sep 2020): text lines 118-119 and Table 1 row 'D+3nu (with new LUNA data) 0.02233 +- 0.00036'. My extracted file mossa/mossa_R2.txt has these at lines 277-281 and 927. IDEAS/Nature abstract: '1.6 percent level' (0.00036/0.02233 = 1.61%). Results: ratio m* r/m_b = 59.9511, omega_h,min = 1.33870808, Hfloor = 116.663537, conservative Hfloor = 113.251818. These match REPORT.md to every printed digit. Using m_b = 0.9315 (amu) raises omega_h,min by 0.698% and Hfloor by 0.343% (117.0636), which strengthens the exclusion. The BBN eta conversion (eta10/Omega_b h^2 = 273.3-274.0, with n_gamma = 410.727 cm^-3) implies 0.9363-0.9387 GeV per baryon, so 0.938 is consistent and 1.01 m_b = 0.947 covers that range. n_B = 2.50837e-7 cm^-3 matches Codex's 2.5083693e-7.

**Scope or correction.** None needed.

### 3. [CONDITIONAL] Angular-distance envelope D_A(z) <= 2c/(Hfloor(1+z)) [1 - (1+z)^(-1/2)] follows from H(z) >= Hfloor (1+z)^(3/2)

**Deciding evidence.** Derived: chi = int c dz/H <= (c/Hfloor) int (1+z)^(-3/2) dz = (2c/Hfloor)[1 - (1+z)^(-1/2)]. The closed form matches numerical Einstein-de Sitter quadrature to <= 3.6e-14 for z in {0.002...1}. 3000 random flat models with extra Lambda, radiation, w in [-1,1] and matter (all rho >= 0) stay below the envelope (maximum ratio - 1 = -0.151). Open curvature (Omega_k > 0, using sinh) also stays below (maximum ratio - 1 = 9.8e-15, equality at K = 0). The bound fails only through components with negative energy density.

**Scope or correction.** Holds for flat or open FRW with every non-floor component nonnegative. I quantified the two escapes for the six-maser check at the conservative floor. (1) Total rho_Lambda <= -0.4759 rho_c/h^2 = -3.85e-47 GeV^4 = -(2.49 meV)^4, i.e. Omega_Lambda <= -0.590. (2) Closed curvature with Omega_k <= -0.563. If the 3 sigma boxes on both sides are required, the negative-Lambda window is Omega_Lambda in [-3.30, -0.590], H0 in [54.6, 89.8] and q0 in [1.39, 5.46].

### 4. [PASS] Pesce et al. 2020 Table 1 transcription (six galaxies, CMB-frame optical velocities)

**Deciding evidence.** I transcribed the table from pdftotext -layout of the supplied PDF, which is arXiv:2001.09213v2 (sha256 c9f2d2b6...; layout lines 208-226). I cross-checked it against the arXiv v2 LaTeX source main.tex lines 124-131 (identical). Rows: UGC 3789 51.5 +4.5/-4.0, 3319.9 +- 0.8; NGC 6264 132.1 +21/-17, 10192.6 +- 0.8; NGC 6323 109.4 +34/-23, 7801.5 +- 1.5; NGC 5765b 112.2 +5.4/-5.1, 8525.7 +- 0.7; CGCG 074-064 87.6 +7.9/-7.2, 7172.2 +- 1.9; NGC 4258 7.58 +- 0.11, 679.3 +- 0.4. Table note (main.tex line 131): 'CMB reference frame using the optical convention'.

**Scope or correction.** None.

### 5. [PASS] Six-row comparison: primary corner (v + 3 sigma_v, v_pec = -750, multiplicative relativistic Doppler) mismatches UGC 3789 and NGC 5765b; the -1500 stress case leaves only NGC 5765b

**Deciding evidence.** Blind values (D_low vs D_A,max in Mpc, -750 / -1500). UGC 3789: 39.50 vs 35.197 / 41.578 (X / o). NGC 6264: 81.10 vs 90.967 / 96.995 (o / o). NGC 6323: 40.40 vs 72.080 / 78.228 (o / o). NGC 5765b: 96.90 vs 77.841 / 83.953 (X / X). CGCG 074-064: 66.00 vs 67.028 / 73.208 (o / o). NGC 4258: 7.25 vs 12.549 / 19.073 (o / o). These are identical to all 18 REPORT.md table entries at 0.01 Mpc. Flipping NGC 5765b to consistent needs v_pec <= -3105 km/s (Pesce Table 4's most negative flow-model value for it is -460), or a distance excursion of 6.74 lower widths (-750) or 5.54 (-1500). UGC 3789 flips at v_pec = -1255 km/s. CGCG 074-064 passes by only 127 km/s (margin -1.03 Mpc).

**Scope or correction.** The corners are sensitivity corners, not confidence limits, as Codex states. Using 3x the 16-84 width assumes the tails are not much heavier than Gaussian. I could not test that without posterior samples.

### 6. [PASS] Necessary product bound mX|YQ0/YB0| <~ 63.9344 GeV (primary) or 74.6429 GeV (stress); benchmark 100 GeV excluded

**Deciding evidence.** The most constraining row is NGC 5765b. At -750: H_req = 90.98 km/s/Mpc, m*|r| <= 35.953 GeV, mX|r| <= 63.934 GeV. At -1500: H_req = 98.12, m*|r| <= 41.975, mX|r| <= 74.643. Both match Codex. Other floor choices: nominal floor 60.16; amu m_b 59.74; omega_b - 5 sigma with 1.02 m_b 66.89. All are below 100.

**Scope or correction.** The coupling-free statement is on m*|r| (<= 35.95 / 41.97 GeV). The mX form divides by mu = 0.5623 and so holds only at fixed benchmark couplings, as REPORT says.

### 7. [PASS] No defensible alternative velocity-frame or redshift-composition convention overturns the exclusion

**Deciding evidence.** I computed shifts from the Planck dipole (369.82 km/s toward l,b = 264.021, 48.253; RA/Dec 167.942/-6.944) and SIMBAD ICRS coordinates. v_CMB - v_hel = +60.5 / -0.4 / -32.7 / +206.9 / +259.0 / +205.5 km/s (UGC 3789 ... NGC 4258). Local Group frame (KM96, quoted from memory) shifts v_LG - v_CMB = +51.5 / +170.5 / +254.1 / -232.7 / -310.5 / -146.0. Effect on z_cos*c for NGC 5765b at -750: Pesce Eq. 9 composition +1 km/s; additive composition -22; relativistic reading of v +125; radio reading +250; heliocentric -207; LG -233. The flip needs +2431. The most favourable combination (radio reading + best of CMB/helio/LG frame + additive composition) still mismatches UGC 3789 (3.63 Mpc) and NGC 5765b (17.26 Mpc), with bounds of 67.01 GeV (-750) and 77.57 GeV (-1500).

**Scope or correction.** The only row that flips under any convention is CGCG 074-064. The heliocentric and LG frames turn it into a mismatch (+1.07 and +1.48 Mpc at -750), which strengthens the exclusion. The CMB frame is the physically appropriate one for FRW z_cos.

### 8. [CONDITIONAL] The 100 GeV, r = 1 preparation is excluded (conditional on leading classical order, flat FRW, standard BBN and nonnegative other mean densities)

**Deciding evidence.** The exclusion survives a total rho_Lambda down to -3.85e-47 GeV^4 (Omega_Lambda >= -0.590), so the 'nonnegative Lambda' premise can be relaxed to this quantified threshold. For comparison, the s08 one-loop new-scalar vacuum term quoted in issue #13 (-1.22652632365e-6 v^4) is -1.2265e6 GeV^4 at v = 1000 GeV. That is 3.18e52 times the escape threshold. Any residual total vacuum energy with Omega_Lambda in [-3.30, -0.590] lets all six maser boxes intersect the model, but it requires q0 >= 1.39.

**Scope or correction.** To make the classical exclusion order-consistent with the one-loop ledger, the total vacuum counterterm must be stated as a matching choice that cancels the new-scalar vacuum energy to below 3.85e-47 GeV^4 (about 1 part in 3e52), or it must be calibrated independently. Data that independently require q0 < 0 (SN Ia deceleration) can calibrate 'Omega_Lambda > -0.59' without using the masers. I did not run that check.

## Strongest agreement

My blind reimplementation reproduces every number Codex publishes for this step, from the action and primary sources, with no shared code. That covers m* = 56.23413252 GeV (mu = 0.562341325), the two Hfloor values, all 18 envelope entries in the REPORT.md table and both product bounds (out/bcd_envelope.json; values in the numbers list). I wrote the m* bound's derivation myself, and it is sharp: rho/|nQ| reaches m* exactly at phi = 0, |X| = 1.2574 v, with 0 violations in 3e6 samples. The decisive row, NGC 5765b, is robust. Flipping it needs v_pec <= -3105 km/s (the flow models give at most -460) or a 6.7-width distance excursion. No frame or redshift convention I tested (heliocentric, Local Group, radio, relativistic, additive or Pesce Eq. 9) moves it by more than +250 km/s.

## Strongest defect or narrowing

I found no numerical defect. The strongest narrowing concerns the premise that all other mean densities are nonnegative, which is load-bearing, and I quantified it. The six-maser exclusion of 100 GeV, r = 1 disappears if the total vacuum density is <= -3.85e-47 GeV^4 = -(2.49 meV)^4 (Omega_Lambda <= -0.590, H0 = 89.8, q0 >= 1.39), or if space is closed with Omega_k <= -0.563. Open curvature and any positive component do not weaken it (checked to 1e-14). This answers Codex's question in issue #13: the classical exclusion is order-consistent with the s08 one-loop term (-1.2265e6 GeV^4 at v = 1 TeV) only if the vacuum counterterm is stated explicitly as a matching choice that cancels it to roughly 1 part in 3e52, or if Omega_Lambda > -0.59 is calibrated independently, for example from SN Ia evidence that q0 < 0. A second, minor premise: the floor proof also needs mD0^2 >= 0 and lambda_D >= 0, which the benchmark satisfies.

## Key numbers

- m*: 56.23413252 GeV (mu = m*/mX = 0.5623413252) (a_charge_floor.py analytic; optimizer relative difference -1.1e-16; grid 8.6e-6)
- saturating configuration: phi = 0, |X|/v = 1.2574334, Xdot = -i sgn(nQ) m* X; min rho/(m*|nQ|) = 1.0000000000000 (out/a_charge_floor.json)
- pointwise violations: 0 / 2,000,000 generic; 0 / 1,000,000 adversarial (out/a_charge_floor.json)
- omega_b (Mossa+2020, D+3nu with LUNA): 0.02233 +- 0.00036 (mossa/mossa_R2.txt lines 277-281 (ms lines 118-119) and line 927 (Table 1))
- omega_h,min nominal: 1.3387080801 (out/bcd_envelope.json floor)
- Hfloor nominal / amu m_b / conservative: 116.66353672 / 117.06357717 (+0.343%) / 113.25181752 km/s/Mpc (out/bcd_envelope.json floor)
- mass per baryon implied by eta10/Omega_b h^2 in 273.28-274.0: 0.93874-0.93627 GeV (n_gamma = 410.727 cm^-3) (out/bcd_envelope.json)
- D_low vs D_A,max at -750 (Mpc): UGC3789 39.50 > 35.197 X; NGC6264 81.10 < 90.967; NGC6323 40.40 < 72.080; NGC5765b 96.90 > 77.841 X; CGCG074-064 66.00 < 67.028; NGC4258 7.25 < 12.549 (out/bcd_envelope.json nominal_vpec_m750)
- D_A,max at -1500 (Mpc): 41.578, 96.995, 78.228, 83.953 (NGC5765b still mismatched by 12.95 Mpc), 73.208, 19.073 (out/bcd_envelope.json stress_vpec_m1500)
- bound from NGC 5765b at -750 / -1500: m*|r| <= 35.953 / 41.975 GeV; mX|r| <= 63.934 / 74.643 GeV (benchmark 100) (out/bcd_envelope.json, out/e_addenda.json)
- v_pec needed to flip NGC 5765b / UGC 3789 / CGCG 074-064: -3105 / -1255 / -626 km/s (Pesce Table 4 most negative: -460 / -142 / -372) (out/e_addenda.json)
- NGC 5765b distance-width multiple to flip: 6.74 (at -750), 5.54 (at -1500) (out/e_addenda.json)
- v_CMB - v_helio: +60.5, -0.4, -32.7, +206.9, +259.0, +205.5 km/s (table order) (out/bcd_envelope.json frame_shifts_kms)
- most favourable convention bound: mX|r| <= 67.01 GeV (-750), 77.57 GeV (-1500) (out/bcd_envelope.json variant_most_favourable_*)
- negative-Lambda escape threshold: H_L^2 = -4758.8 (km/s/Mpc)^2, Omega_Lambda = -0.590, rho_Lambda = -3.853e-47 GeV^4 = -(2.491 meV)^4, H0 = 89.82 (out/bcd_envelope.json premise_escape_negative_Lambda)
- negative-Lambda window (both 3 sigma boxes): Omega_Lambda in [-3.30, -0.590], H0 in [54.6, 89.8], q0 in [1.39, 5.46] (out/bcd_envelope.json negative_Lambda_window)
- closed-curvature escape: Omega_k <= -0.563 (H0 = 90.59) (out/bcd_envelope.json premise_escape_closed_curvature)
- ratio |one-loop rho_vac,new(v = 1 TeV)| / escape threshold: 3.18e52 (inline computation from issue #13 value -1.22652632365e-6 v^4)

## Scripts and outputs

- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/a_charge_floor.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/out/a_charge_floor.json (python3 -I a_charge_floor.py; derivation in docstring)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/bcd_envelope.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/out/bcd_envelope.json (floor, envelope checks, six rows at -750/-1500/central, 7 convention variants, negative-Lambda/closed-curvature escapes)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/e_addenda.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/out/e_addenda.json (flip thresholds, open-curvature check, floor variants)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/out/rows_blind.csv (all 114 evaluated rows)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/out/BLIND_FREEZE_sha256.txt (hashes of scripts and outputs, frozen 2026-10-20T20:44:41Z before reading REPORT.md)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/pdf/pesce_layout_mine.txt (my pdftotext -layout; Table 1 at lines 208-226)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/arxiv_src/x/main.tex (arXiv 2001.09213v2 source; Table 1 lines 124-131, Table 4 lines 352-363)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/mossa/mossa_R2.txt (Mossa R2 manuscript text; sha256 of PDF 3c1d6b36...)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/blind_envelope/simbad/coords.csv and simbad/cgcg.txt (ICRS coordinates)

## Not done

- I did not read the Nature published full text (paywalled). I used the accepted manuscript R2 from the Padua repository and the Nature/IDEAS abstract (1.6% precision); the published value was not compared digit by digit.
- I did not test the heavy-tail risk of the 3-width distance corners: Pesce posterior samples were not obtained, and the Gaussian-width extrapolation of the 16-84 percentiles is assumed.
- The CMB dipole (Planck 2018) and Local Group apex (Karachentsev and Makarov 1996) parameters come from memory and were not fetched. The conclusions do not depend on them, since flipping a decisive row needs a shift of at least 511 km/s for UGC 3789 and 2431 km/s for NGC 5765b.
- I did not run an independent SN Ia or BAO calculation to calibrate Omega_Lambda > -0.59 independently of the masers.
- I did not audit quantum or one-loop corrections to the floor, gravitational binding or backreaction, or the 3D dynamical durability claims (outside this key's scope).
- I did not commit, push or post anything to issue #13. I did not run Codex's verifier (audit.py) and did not read any of the forbidden files.
