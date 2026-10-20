# Second-agent review: Adversarial derivation and scope review

Separate Claude agent spawned by the main session for issue #13; it returned structured findings, rendered here verbatim. Scripts and outputs are in this directory; paths inside them refer to the session scratchpad. Reviewer: derivation_scope (Claude Code subagent). Adversarial audit of DERIVATION/ASSUMPTIONS at 2ee7f82a62287ae67d92b3b481646f8efed20fc0, checked against the in-tree Candidate 1 quantum and energy packets (identical at 2ee7f82 and fd0de7d; git diff --stat is empty) and the s08 ledger at fd0de7d.

## Findings

### 1. [PASS] (A) Classical normalization and pointwise bound, DERIVATION.md:9-84. The bound reads n_Q = q2 q1dot - q1 q2dot = omega f^2, Q = Qbar/gX^2, E = (v/gX) Ebar, min 2V/F^2 = min(1, sqrt(ab)), with square completion for both signs, giving rho_h >= m*|n_Q| and m* = v min(gX, (lambda_phi lambda_X)^(1/4)).

**Deciding evidence.** Checked by hand. j^mu = i(X dX* - X* dX) with X = (q1 + i q2)/sqrt2 gives j = -(q2 dq1 - q1 dq2), so n_Q = q2 q1dot - q1 q2dot, and X = f e^{-i omega t}/sqrt2 gives omega f^2. U = gX^2 v^4 Vbar and |X|^2 = v^2 F^2/2, so U/|X|^2 >= mX^2 min(1, sqrt(ab)). The infimum of a(s-1)^2/(2y) + b y/2 over y is sqrt(ab)|s-1|. The square-completion identity expands exactly. At the benchmark ab = 0.1 < 1, so the minimizer is S = 0, F^2 = sqrt(10). Independent code (cw_floor.py) gives tree m*^2/v^2 = 0.0031622776601684 at (phi = 0, R/v = 1.7782794100), so m* = 56.23413252 GeV, matching PREDICTION.json mstar_GeV = 56.23413251903491.

**Scope or correction.** None needed. This is a leading-order classical, minimally coupled statement, as DERIVATION.md:3-5 says.

### 2. [CONDITIONAL] (A) The step from the slice bound to H0 >= 100 sqrt(omega_b + omega_h,min) and to the D_A envelope (DERIVATION.md:86-184). This covers the FRW slice, gradients, averaging, conservation of charge and entropy, the omega_b convention and binding.

**Deciding evidence.** (1) The gradient terms h^ij d_i q d_j q are >= 0 on any spacelike slice, and the bound uses the same normal n for rho and n_Q, so it holds frame by frame. (2) Averaging gives <rho_h> >= m*<|n|> >= m*|<n>|. A^3 n_Q is constant under exact U(1) on a closed comoving domain. With present r, the entropy cancels exactly. (3) omega_b convention: the audit's mb = 0.938 corresponds to eta10/omega_b = 273.50 (n_gamma = 410.718 cm^-3). The common 273.9 and 273.279 factors imply mb = 0.93663 and 0.93876 GeV, which changes n_B by +0.146% and -0.081%. Both are inside the 1% upward-mb allowance, which lowers n_B (scope_checks_output.json). (4) If all hidden matter is clumped (Dyer-Roeder empty beam, D_A = c int dz/((1+z)^2 H)) with the same H floor, D_max for NGC 5765b (stress case) rises from 83.9531 to 83.9760 Mpc (+0.0273%). The ratio D_low/D_max stays 1.1539. (5) To remove the stress-case mismatch, 25.357% of the floor hidden energy would have to be missing from the z~0 Friedmann source. Weak-field binding cannot do that: the reference object has R_s/R90 = 9.6e-31 (EXCLUSION_AUDIT.json). It would take compactness near black-hole scale, or loss of charge into black holes that then evaporate.

**Scope or correction.** The derivation is correct under assumptions C1-C5. Binding/collapse and charge loss through black holes are the only physical routes that close a 25% gap, and they are outside the EFT, as DERIVATION.md:113-114 already states. Light-path representativeness (C6) is not load-bearing for focusing at z <= 0.034.

### 3. [PASS] (A) Nonminimal xi R|X|^2: the matching zeros hold only at mu0, so xi runs. Can the improvement terms push the averaged hidden energy below m*|n|?

**Deciding evidence.** With L containing -xi R|X|^2, T_munu gains 2 xi (G + g Box - nabla nabla)|X|^2. On flat FRW this gives rho_xi = 2 xi [3H^2 psi + 3H psi_dot - lap(psi)/a^2], with psi = |X|^2. The FRW average removes the Laplacian. Dilution then gives <rho_xi> = -12 xi H^2 <psi>, so <rho_xi>/(m* n) = 6 xi H^2/(omega m*). The one-loop result 16pi^2 beta_xi_ab = lambda_abcd (xi_cd - delta_cd/6) gives beta_xi_X = -2.955e-5 per e-fold, so |xi_X(m*)| = 1.70e-5. Even with xi = 1e-3 and H = H_floor = 2.42e-42 GeV, the averaged ratio is 7.8e-90 (Q-ball, omega = 0.8 mX) or 1.1e-89 (Q-matter). The Jordan-frame shift 2 xi |X|^2/MP^2 inside a core is 5.3e-34. Pointwise, the wall term -2 xi lap(psi) is about 6.6e-5 relative.

**Scope or correction.** The pointwise statement (DERIVATION.md:72-74, 88-89) is exact only at xi = 0. With running xi it can fail locally at about 1e-4 relative. The averaged floor that enters Friedmann is unaffected to better than 1e-88. DERIVATION.md:115-116 lists nonminimal terms as needing 'a new calculation'; this result closes that item for the averaged floor.

### 4. [PASS] (A) D-sector sign conditions at the benchmark.

**Deciding evidence.** GRUT_CANDIDATE_1.md:211 and PARAMETER_PROVENANCE.md:14 give mD0 = gX v = 100 GeV, so mD0^2 = +1e4 GeV^2, with hD^2 = 1e-4 and lambda_D = 1e-3, all >= 0. The generated portal |alpha| <= 1.758e-8 for mu/mu0 in [1/4, 4] (quantum REPORT.md:75-79) keeps the D quadratic positive for R/v < 754 and leaves quartic boundedness intact (alpha^2 = 3.1e-16 << lambda_X lambda_D = 1e-6). At one loop, a numerical minimization of the ratio with D left free returns D = -1.4e-8, i.e. zero (cw_crosscheck_output.json).

**Scope or correction.** ASSUMPTIONS.md A1 says 'positive couplings'. The bound also needs mD0^2 >= 0, which is a mass parameter rather than a coupling. Add it explicitly.

### 5. [CONDITIONAL] (B) Vacuum-energy matching. Which vacuum energy gravitates, and is 'nonnegative Lambda' (DERIVATION.md:119-131, ASSUMPTIONS.md C3) consistent with Candidate 1's own matching?

**Deciding evidence.** In Friedmann, the total renormalized vacuum density gravitates: rho_Lambda + C + V_eff,new(vac) + V_eff,SM(vac), and only that total is mu-independent (quantum REPORT.md:134: 'enter gravity through their sum'). The floor (U(v) = 0) and the fixed-charge energies are vacuum-subtracted (energy REPORT.md:51-55; quantum REPORT.md:138). If the benchmark matching rho_Lambda(mu0) = 0 is kept literally with fixed coefficients, the in-tree one-loop new-scalar term is -1.2265263237e-6 v^4 (reproduced: -1.2265263236516e-6), i.e. -1.2265e6 GeV^4 at v = 1 TeV. That alone is an AdS radius of 0.751 m (2.5 ns), with turnaround at radiation temperature T = 13.7 GeV, before BBN. SM one-loop terms in MSbar/Landau gauge at mu = 100 GeV, with tree V_SM(vac) = 0, add about +5.6e6 to +5.8e6 GeV^4 (top-dominated, mt = 163-172.5 GeV). The literal total is then about +4.4e6 GeV^4, i.e. de Sitter with H of order 5e-16 GeV. At mu = 56 GeV the SM sign flips (-8.4e6 to -1.2e7). Either way |rho_vac| is about 4.5e52 times the observed 2.7e-47 GeV^4. The audit packet never mentions this: grep for 'vacuum energy' or '1.22' in DERIVATION, ASSUMPTIONS and REPORT returns nothing.

**Scope or correction.** 'Nonnegative total vacuum energy' cannot come from the Candidate 1 benchmark (GRUT_CANDIDATE_1.md:55; A4). It is a separate cosmological matching choice, calibrated externally to the observed late acceleration and fine-tuned against the new-scalar and SM loops to about 1 part in 1e52. Re-label C3 as an E-class input: 'total renormalized vacuum energy rho_vac,tot calibrated to observed acceleration; the exclusion requires only rho_vac,tot >= -2.72e-47 GeV^4' (see the next finding). This narrows the exclusion's domain without overturning it.

### 6. [PASS] (B) Does the exclusion survive a small negative total Lambda? Critical value for the NGC 5765b stress case.

**Deciding evidence.** I reimplemented the envelope with H^2 >= Hf^2 (1+z)^3 - L, using Hf = 113.25181751823 (conservative) and the COMPARISON.json z_high corners and D_low values, and solved D_A,max(z_high; L) = D_low with brentq. At L = 0 it reproduces the analytic envelope to 1e-11. NGC 5765b, stress case: critical omega_Lambda = -0.33585, i.e. Omega_Lambda = -0.3547 relative to the model's H0 = 97.30 km/s/Mpc. That is rho_vac = -2.719e-47 GeV^4 = -(2.284 meV)^4, and needs q0 = +1.03 at the rescue point. NGC 5765b, primary case: omega_Lambda = -0.4759 (Omega_Lambda = -0.590). UGC 3789, primary case: omega_Lambda = -0.2696 (rho = -2.18e-47 GeV^4).

**Scope or correction.** The stress-case rejection holds for any rho_vac,tot > -2.72e-47 GeV^4, which includes every positive calibration (observed about +2.7e-47). A rescue needs a negative vacuum density of the same order as the observed dark energy, combined with q0 of about +1. Supernova data (q0 < 0) would exclude that, but that is an added external datum and is not used in the frozen protocol.

### 7. [PASS] (C) One-loop correction to the homogeneous floor m*_eff^2 (MSbar, mu = mX = 0.1 v, fields phi, x, y, D, with the benchmark D parameters).

**Deciding evidence.** Static Coleman-Weinberg, strict first order: delta m*^2 = +1.1787209183e-7 v^2 (+3.7274e-5 relative). This matches the in-tree quantum REPORT.md:109 value (1.17872091832e-7) to 10 digits, and a finite-difference Hessian cross-check gives 1.17872091e-7. Full numerical minimization of the truncated V0 + V1 gives +3.7282e-5. A better homogeneous estimate is the rotating Q-matter one-loop grand potential at chemical potential mu = m*. It needs no derivative expansion: it includes gyroscopic x-y mixing and the Goldstone mode, with a UV remainder falling as k^-5 (checked to 1e-9), and it gives exactly zero for a free complex scalar. Result: delta m*^2 = +1.0806e-7 v^2 (+3.417e-5), so delta m*/m* = +1.709e-5. Mass-squared signs: at the tree minimizer the Hessian eigenvalues (0.00316, 0.00949, 0.0100, 0.0216) v^2 are all positive. The rotating-frame matrix is diag(0.00632, 0) >= 0. Along the tree ratio-minimizing path R^2 = sqrt(10)(1 - phi^2), the phi-x block is negative for phi/v in [0.45, 0.895] (minimum -0.00335 v^2), but there the tree ratio is >= 0.004134 v^2 = 1.307 m*^2, well away from the minimum. Removing the stress-case mismatch needs m* lower by 25.357% (primary case: 36.07%; UGC 3789 primary: 20.95%). The loop shift is about 1.4e4 times smaller and has the opposite sign, raising the floor. Treating the generated m_X^2(4 mu0) = -1.757e-6 v^2 as an uncancelled shift gives -0.056% in m*^2.

**Scope or correction.** Both results are perturbative one-loop estimates (constant-field or homogeneous Q-matter). Neither is a quantum energy inequality: local renormalized energy densities can be negative. Two-loop, graviton and SM-loop contributions to the floor are not included. The quantum scope therefore narrows from 'uncomputed' to 'one-loop shift of about +2e-5 in m*, about 1e4 times below the margin'.

### 8. [DEFECT] (D) The audit packet calls loop and vacuum effects uncomputed, although in-tree packets at the same commit computed them.

**Deciding evidence.** UNCERTAINTY_BUDGET.md:27 says 'Loop/running/higher-operator changes and renormalized stress corrections are uncomputed'. PARAMETER_PROVENANCE.md:35-37 says 'Vacuum/loop effects at higher order are a separate, uncomputed question'. DERIVATION.md:115-117 says running/loop changes 'require a new calculation'. Yet research/reality_framework/grut_candidate_1_quantum_20261018/REPORT.md:97-98, 109 and 128-138 are present at 2ee7f82, unchanged (git diff --stat 2ee7f82 fd0de7d over the quantum, energy and candidate packets is empty). They give the one-loop floor shift (+0.0037%) and the -1.2265e-6 v^4 vacuum obligation with 'No adjustment has been applied here'.

**Scope or correction.** Cite the in-tree one-loop floor shift and quantify the EFT row: about +2e-5 in m* against the 0.2536 needed. Import the total-vacuum matching obligation into ASSUMPTIONS, as in the vacuum-matching finding above. This correction strengthens the exclusion's quantum robustness and narrows its cosmological premise.

### 9. [PASS] (D) Wording and over-claim checks in DERIVATION and ASSUMPTIONS: the envelope monotonicity, the product bound, peculiar-velocity composition and the hidden-pressure independence.

**Deciding evidence.** D_A,max is proportional to u^2(1 - u) with u = (1+z)^{-1/2}. Its derivative u(2 - 3u) changes sign at z = 1.25, so it increases for z < 1.25 (DERIVATION.md:173). The product bound at DERIVATION.md:177-179 follows algebraically from D <= 2c[...]/(Hf(1+z)). The Doppler composition at line 168 is the exact radial relativistic form. The z = 0 H0 floor needs no equation of state. H(z) >= Hf(1+z)^{3/2} needs only conserved n_Q and n_B plus Lambda_tot >= 0 and flatness, all stated. DERIVATION.md:128-131 already lists negative Lambda and positive curvature as outside scope.

**Scope or correction.** No over-claim found in the inequalities themselves. The one gap is that 'cosmological constant' at DERIVATION.md:120 and C3 should say 'total renormalized vacuum energy', with its calibration source named.

### 10. [CONDITIONAL] Reconciling the s08 ledger with the 2ee7f82 audit. Neither contradicts the other, but each omits the other's load-bearing item.

**Deciding evidence.** grep -i 'maser|H0|exclu|2ee7f82|479f111|Pesce|BBN' over ledger_s08_candidate1.md finds no empirical-exclusion record. The ledger cites the closure-fork branch at 96a7e7e/1489773 (s08.md:4). s08.md:110, 116 and 126 state that the gravitating total vacuum energy is 'not fixed by the model' and that its matching is free. That is exactly the premise the audit's C3 needs and does not import.

**Scope or correction.** The s08 addendum should record the 2ee7f82 conditional exclusion and the 3D survival results. The audit should record that its C3 is the external total-vacuum calibration that s08 identifies as unfixed. The critical bound rho_vac,tot >= -2.72e-47 GeV^4 (stress case) is the shared quantitative link.

## Strongest agreement

The classical chain from rho_h >= m*|n_Q| to H0 >= 100 sqrt(omega_b + omega_h,min) and to the D_A envelope is algebraically correct at the stated order. I checked current sign, normalization, potential minimization, square completion, FRW gradients, averaging, charge and entropy bookkeeping, the omega_b convention (n_B within -0.08%/+0.15%), the product bound and envelope monotonicity. The exclusion margin is large compared with every correction I could compute. The stress case needs m* (or the hidden floor energy) lower by 25.357%. The one-loop shift is +1.71e-5 (rotating Q-matter) or +1.86e-5 (static CW, reproducing the in-tree +1.17872091832e-7 v^2 to 10 digits); that is about 1.4e4 times smaller and in the direction that raises the floor. Running nonminimal xi changes the averaged floor by at most 1e-89. Clumped-matter (empty-beam) propagation changes D_max by +0.027%. The benchmark D sector is sign-safe (mD0^2 = +1e4 GeV^2, lambda_D = 1e-3).

## Strongest defect or narrowing

The exclusion holds only with an explicit total-vacuum matching that Candidate 1 does not supply, and the audit does not state this. The H0 floor and the envelope need total renormalized vacuum energy rho_vac,tot >= -2.72e-47 GeV^4 = -(2.28 meV)^4 (NGC 5765b stress case; omega_Lambda >= -0.336, Omega_Lambda >= -0.355 relative to the model H0 = 97.3). Candidate 1's own benchmark matching (rho_Lambda(mu0) = 0) with fixed one-loop coefficients gives -1.2265e6 GeV^4 from the new scalars at v = 1 TeV. That is an AdS radius of 0.75 m, and recollapse would come at about 14 GeV, before BBN. SM one-loop terms (about +5.6e6 GeV^4 at mu = 100 GeV in MSbar) flip the literal sign, so |total| is about 1e52-1e53 times the observed dark energy. Nonnegative total vacuum energy is therefore a separate cosmological matching, tuned and calibrated to the observed acceleration. It should be an E-class input in ASSUMPTIONS rather than an unlabeled 'nonnegative Lambda' under C3. Once that calibration is adopted, any rho_vac,tot > -2.72e-47 GeV^4 keeps the benchmark excluded. The packet also wrongly describes the in-tree one-loop floor and vacuum results as 'uncomputed' (UNCERTAINTY_BUDGET.md:27, PARAMETER_PROVENANCE.md:35-37).

## Key numbers

- m* (tree): 56.23413252 GeV; m*^2/v^2 = 0.0031622776601684 at phi = 0, R/v = 1.7782794100 (cw_floor_output.json; matches PREDICTION.json mstar_GeV)
- strict one-loop new-scalar vacuum energy: -1.2265263236516e-6 v^4 = -1.2265e6 GeV^4 at v = 1 TeV (in-tree value -1.22652632365e-6) (cw_floor_output.json; quantum REPORT.md:98,131)
- static CW delta m*^2 (strict): +1.1787209183e-7 v^2 = +3.7274e-5 relative (finite-difference check 1.17872091e-7; numerical truncated minimization +3.7282e-5) (cw_floor_output.json; cw_crosscheck_output.json; in-tree REPORT.md:109 gives 1.17872091832e-7)
- rotating Q-matter one-loop delta m*^2: +1.0806499e-7 v^2 = +3.4173e-5 relative; delta m*/m* = +1.7087e-5; finite remainder -1.5506e-8 v^4 (cw_floor_output.json rotating_Qmatter)
- required fractional reduction of m* (or hidden floor energy): stress NGC 5765b 0.25357; primary NGC 5765b 0.36066; primary UGC 3789 0.20949 (cw_floor_output.json from COMPARISON.json product limits 74.6429/63.9344/79.0507 GeV)
- margin / one-loop shift in m*: 13606 (static), 14840 (rotating) (cw_floor_output.json)
- critical negative vacuum, NGC 5765b stress: omega_Lambda = -0.33585; Omega_Lambda = -0.35473; model H0 = 97.301; q0 = 1.032; rho_vac = -2.7190e-47 GeV^4 = -(2.2835 meV)^4 (vacuum_lambda_output.json)
- critical negative vacuum, NGC 5765b primary: omega_Lambda = -0.47588; Omega_Lambda = -0.58991; rho_vac = -3.853e-47 GeV^4 (vacuum_lambda_output.json)
- critical negative vacuum, UGC 3789 primary: omega_Lambda = -0.26963; Omega_Lambda = -0.26617; rho_vac = -2.183e-47 GeV^4 (vacuum_lambda_output.json)
- literal-matching cosmology (new scalars only): L_AdS = 0.751 m (2.51e-9 s); radiation turnaround T = 13.7 GeV; |rho_new|/rho_obs = 4.48e52 (vacuum_lambda_output.json)
- SM one-loop vacuum (MSbar, Landau, tree V_SM(vac) = 0): +5.64e6 to +5.77e6 GeV^4 at mu = 100 GeV; -8.4e6 to -1.24e7 at mu = 56.2 GeV; +1.9e7 to +2.3e7 at mu = 172.5 GeV (vacuum_lambda_output.json)
- beta_xi at xi = 0 (per e-fold): X: -2.955e-5; phi: -1.058e-4; D: -6.54e-6; |xi_X(m*)| = 1.70e-5 (scope_checks_output.json)
- averaged xi-improvement relative to m* n (xi = 1e-3, H = H_floor): 7.8e-90 (Q-ball) / 1.1e-89 (Q-matter); Jordan shift 5.3e-34; pointwise wall about 6.6e-5 (scope_checks_output.json)
- empty-beam vs filled-beam D_max (NGC 5765b stress): 83.9760 vs 83.9531 Mpc (+0.0273%); D_low/D_max = 1.1539 (scope_checks_output.json)
- omega_b to n_B convention: eta10/omega_b = 273.50 for mb = 0.938; 273.9 corresponds to mb = 0.93663 (+0.146% n_B), 273.279 to mb = 0.93876 (-0.081% n_B) (scope_checks_output.json)
- tree path negative-curvature window: phi/v in [0.45, 0.895], minimum eigenvalue -0.003347 v^2; tree ratio there >= 0.004134 v^2 = 1.307 m*^2 (cw_floor_output.json)

## Scripts and outputs

- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/cw_floor.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/cw_floor_output.json (static CW floor, numerical minimization, negative-curvature map, rotating Q-matter one-loop grand potential, margin ratios)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/cw_crosscheck.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/cw_crosscheck_output.json (finite-difference Hessian cross-check, D left free)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/vacuum_lambda.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/vacuum_lambda_output.json (critical negative Lambda per row and case, literal-matching AdS/dS consequences, SM one-loop vacuum budget)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/scope_checks.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/scope_checks_output.json (beta_xi and averaged improvement term, D-sector signs, omega_b/eta convention, empty-beam distance, removal fraction)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/derivation_scope/s08.md (verbatim extract of fd0de7d ledger_s08_candidate1.md used for grep evidence)
- All run with: python3 -I <script> (numpy/scipy only; each under 1 minute). Source files read via git show 2ee7f82...:path and fd0de7d...:path; working tree untouched.

## Not done

- Did not rerun the author's audit.py/quality.py or recheck the Pesce Table 1 transcription or percentile-width corners (outside the derivation key); used COMPARISON.json z_high and D_low values as given.
- Did not compute the critical positive spatial curvature (Omega_k) that would evade the floor; only negative Lambda was computed.
- No two-loop, graviton-loop, curvature-squared or SM-loop contribution to the floor m*; no quantum energy inequality; the one-loop results are perturbative homogeneous estimates only.
- SM vacuum energy is a one-loop MSbar/Landau-gauge estimate with tree masses, for illustration only. It shows the sign and magnitude are scheme- and scale-dependent, not a precise SM vacuum value.
- No BBN recomputation with a hidden energy component, no primordial-black-hole formation or evaporation calculation, and no Buchert backreaction estimate (binding/collapse loopholes are bounded only by the required 25.36% energy removal).
- Did not propose the joint next calculation for issue #13 or post to GitHub; that is left to the orchestrating session.
