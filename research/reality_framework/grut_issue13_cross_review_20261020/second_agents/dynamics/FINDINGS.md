# Second-agent review: 3D dynamics evidence

Separate Claude agent spawned by the main session for issue #13; it returned structured findings, rendered here verbatim. Scripts and outputs are in this directory; paths inside them refer to the session scratchpad. Reviewer: Claude Code (Opus 5.5), key "dynamics": independent audit of research/reality_framework/grut_3d_empirical_20261020 at 2ee7f82a62287ae67d92b3b481646f8efed20fc0. I extracted the raw archives into my own scratch directory and verified all 190 raw files and all archives against RAW_MANIFEST.json (0 mismatches). I recomputed everything with my own scripts and did not import analyze.py. I compiled the frozen cartesian.cpp (sha256 3e8c04f2…) and ran one deciding box control at coarse resolution.

## Findings

### 1. [PASS] (1) Finite-time durability: 26/26 finite-object runs keep at least 90% of initial Q in R12 core balls over tau 120..160 (isolated about 99.85%, pairs about 97%), with resolution, timestep, box and rotation controls.

**Deciding evidence.** Recomputed from the raw CSVs (scripts/independent_dynamics.py, out/independent_dynamics.txt):
- 26/26 runs pass. Every late window is t=120..160 with 9 samples.
- Pairs: minimum joint R12 retention is 96.98103% (pair_negative_fine). Each core holds at least 32.5% of Q0 (criterion 15%). Minimum core separation is 27.05 (criterion >24).
- Perturbed single: 99.85098 / 99.84893 / 99.84842% at h=.75/.5/.375. Late E12/Q12 is at most 0.911367. At the endpoint there is exactly one 2%-of-peak component, holding 94.86-95.16% of Q0.
- Grid differences (max late |dQ12|/Q0): pair_reference 0.21427% then 0.06406% (inferred order p=2.44); pair_generic 0.33654% then 0.10203% (p=2.40). Richardson continuum retention is 97.517% (reference) and 97.103% (generic); the fine grid is within 0.035 points of it.
- Controls: single_box 0.02229%, single_dt 0.00015%, pair_box 0.06678%, pair_dt 0.00271%. The whole-event rotation difference falls from 0.06077% to 0.04318%.
- The result does not depend on the ball radius: minimum R8 retention is 96.228%.
- Every number in the REPORT section 2 tables is reproduced exactly.

**Scope or correction.** The pass applies to the D=0, flat, classical three-field sector, one finite tau window, and small perturbations.
- The single-object seeds inject only +0.186% (generic) and +0.184% (negative) of the energy, about 1.9% of the binding energy Q-E=35.7.
- Pair cores are still shedding charge late in the run, at -0.06 to -0.24 points per 40 tau. Survival beyond tau 240 (see the box-control finding) is extrapolation.

### 2. [PASS] Sub-claim of (1): the conservative fixed-vacuum boundary or box size does not fake retention (charge cannot leave the box, so reflected charge could re-enter the cores).

**Deciding evidence.** Why a control was needed: the frozen pair_box control covered only the reference pair, which moves along z. The transversely deflected negative-seed pair had no box control. Its core0 ends at (-10.31,-11.38,-21.33) moving at (-0.087,-0.095,-0.155), and its R12 ball is only 0.62 from the y face at tau 160.

My deciding control (ctrl/COMMANDS.txt, scripts/box_control_compare.py, out/box_control_compare.txt), coarse h=.75, eps=-0.05:
- Frozen box L24/Lz48: minimum late retention 96.86943%.
- Enlarged box L36/Lz60: 96.71650%, a difference of -0.153 points.
- Per-core late differences are 0.021% and 0.132% of Q0; maximum centre shift is 0.0056.
- Continued to tau 240 in the large box: 96.32335%, slope -0.183 points per 40 tau, wall gap at least 5.77.
- Consistent with this, charge outside the R12 balls at tau 160 is 2.544% in the frozen box and 2.637% in the larger box (pair_reference medium).

Build validation: my binary reproduces the published pair_generic_coarse run to a maximum absolute difference of 3.06e-11 across all columns.

**Scope or correction.** The walls inflate retention by at most about 0.15 points, against a margin of about 7 points above the 90% threshold. Codex's box control should be described as covering the axial pair only. Generic, negative and rotated pairs approach the walls (R12 gaps 5.25, 0.62 and 3.33).

### 3. [CONDITIONAL] Sub-claim of (1): the tested pair cores endure the collision ('outgoing collision cores survive').

**Deciding evidence.** In pair_reference_fine (scripts/extra_checks.py) the two incoming Q=399.4 objects merge before re-splitting:
- At tau 40-55, both half-space core searches converge on one lump: separation 0.00-0.17, Q12=798.79 (99.99% of Q0), maxF2 4.85 against 3.90 for a single object.
- By tau 70-80 it re-emerges as two cores holding 32.5% and 65.0% of Q0 (259.9 and 519.2).

The smaller product is weakly bound:
- Its rest-frame E12/Q12 is 0.9595. My stationary-branch solve gives 0.95991 at Q=258.7 (Omega=0.8483).
- Its binding margin is therefore about 4%, against 8.9% for the original object.
- The stationary branch reaches E/Q=1 near Q of about 180 (Omega about 0.90).

The larger product is still excited: rest-frame E/Q 0.8874 against 0.8828 on the stationary branch.

**Scope or correction.** Suggested wording: 'an equal-charge head-on pair (v=0.2, phase offset pi/2) merges, then re-splits 1:2 into two compact cores holding at least 96.9% of Q through tau 160'.
- The incoming objects do not survive as themselves.
- Only one collision geometry (velocity, phase offset, equal charges) was tested.
- The report should state the merger and the 1:2 re-split explicitly.

### 4. [CONDITIONAL] Sub-claim of (1): the retention criterion and windows were frozen before the runs.

**Deciding evidence.** What matches:
- The sha256 of protocol.json (2aa03cdb…) equals parent_protocol_sha256 in TRANSVERSE_PROTOCOL.json.
- The frozen hashes of cartesian.cpp, run.py, profile.dat, modes.dat and PREFLIGHT.json all match.
- The internal timestamps are consistent: the freeze is at 17:35:16Z and the interruption record at +764 s, with pair_reference_fine at tau 65. Running sequentially, the preceding 10 cases take 374.2 s and the interrupted case reaches tau 65 after about 204 s, so 578 s fits within 764 s.

Limits:
- protocol.json and all raw outputs are in a single commit (479f111, 19:27Z), so git does not attest the order. The freeze rests on the self-reported timestamp.
- The 'one dominant component' rule (>0.8 Q0, runner-up <0.05 Q0) exists only in analyze.py, which was not frozen. This is immaterial: every single run has exactly one component at 94.9-95.2%.
- SPECTRAL_PROTOCOL says it was 'proposed after coarse/medium runs', but its timestamp (+239 s) is before pair_reference_medium could finish running sequentially (+301 s). This is a minor inconsistency.

**Scope or correction.** The freeze is credible and internally consistent but self-attested. Report it as a local freeze only, as Codex already does.

### 5. [PASS] (2) The old transverse canonical amplification (about 9.064) is not a validated internal-instability exponent; the response is chiefly deflection of surviving cores (direct impact x/delta about -14.24 and 7.59, converging).

**Deciding evidence.** Cartesian impact runs (scripts/independent_dynamics.py section F):
- x/delta for coarse/medium/fine is -13.726052/7.333882, -14.114130/7.523045 and -14.241539/7.586161.
- The reference-run x is about 1e-15, so subtracting it introduces no bias.
- Halving the seed changes the response by 0.019% and 0.057%; the box control by 0.127% and 0.227%; the dt control by 0.031% and 0.018%.
- Richardson order is about 2.2, giving continuum values -14.388 and +7.661. The fine grid is within 1.0%.

The late response grows linearly in time, as for rigid deflection:
- d(x/delta)/dt is -0.1217 and +0.0643 per tau. A linear fit over tau 80-160 has maximum relative error 1.5e-3 and 4.3e-3.
- The best exponential-plus-constant fit has rate at most 6e-4 per tau. Reading 9.064 as exponential growth over 160 tau would require 0.0138 per tau.
- Momentum balances: the charge-weighted x_cm/delta is 0.30, against individual core responses of -14.24 and +7.59.

Cross-code check (scripts/old_m1_vs_cartesian.py, scripts/rigid_norm_estimate.py):
- The old cylindrical code's own core-displacement fits are 6.29/10.24/11.93 and -4.65/-6.20/-6.81. They extrapolate to 15.3 and -7.8, matching the Cartesian magnitudes; the sign convention is opposite.
- Rigid translation of the two outgoing cores gives amplifications of 5.371/7.887/8.947. The reported canonical amplifications are 5.395/7.971/9.064, so the difference is at most 1.3%.

**Scope or correction.** This strengthens Codex's argument: the old 9.064 is essentially the norm of the core displacements.

Correction: 9.064 is itself under-resolved. Its Richardson continuum is about 11.15 (p about 1.46). The old symmetry residuals (136.8%, 73.7%, 44.3%) converge only at order 1.5-1.8.

Separately, the frozen Cartesian R36 symmetry-residual criterion fails, which REPORT discloses as not counted. Late pair_reference translation/rotation residuals at h=.75/.5/.375 are 0.2288/0.124, 0.2264/0.116 and 0.2263/0.114, i.e. resolution-independent. In the larger box they are 0.187/0.089, so the residual is boundary-dominated.

### 6. [PASS] (3) The isolated branch is a local fixed-charge energy minimum (numerical Hessian, phase and translation zero modes, angular ordering).

**Deciding evidence.** I re-solved the profile from a generic tanh guess without reading profile.dat (scripts/independent_spectral.py):
- Q=399.438462, E=363.922999, E/Q=0.9110865, F(0)=1.979050, S(0)=0.049080.
- dQ/dOmega=-4093.2, against the author's -4094.03.

My own finite-difference Hessians:
- Fixed-Q l=0 lowest eigenvalue: 0.364348002 / 0.364348691 / 0.364348863 for dr=.1/.05/.025, and 0.36297 at R=60. These are continuum levels near the 0.36 threshold.
- Unconstrained l=0: -0.18886. The rank-one charge term removes this negative direction.
- Translation (l=1): -1.706e-4 / -4.263e-5 / -1.066e-5, shrinking by 4.00 per halving, i.e. O(dr^2).
- Phase: -6.49e-5 / -1.62e-5 / -4.06e-6.
- l=2: 0.285352 / 0.285457 / 0.285483. l=3: 0.37953.
- All match SPECTRAL_RESULTS.json to the printed digits.

Extension to the collision products (scripts/products_branch.py):
- Q about 260: Omega=0.8483, dQ/dOmega=-2029, fixed-Q l=0 lowest 0.2831 (threshold 0.2804), l=2 0.2896.
- Q about 519: Omega=0.7764, dQ/dOmega=-6075, fixed-Q l=0 lowest 0.4002, l=2 0.2494.

**Scope or correction.** This is a re-implementation of the same finite-difference method, not a method-independent or interval-rigorous proof.
- It covers stationary nodeless profiles only, with D=0 and no metric perturbations.
- It does not cover the time-dependent pair.
- The products sit on the stable (dQ/dOmega<0) part of the branch, but the smaller one has only about 4% binding margin.

### 7. [PASS] (4) Conservation quality: maximum charge drift 4.91e-13, maximum Hamiltonian drift 2.76e-5.

**Deciding evidence.** Recomputed over all 29 runs:
- Maximum |dQ/Q0| is 4.907e-13 (pair_negative_fine).
- Maximum |dE/E0| is 2.762e-5 (production_coarse, periodic). For the 26 object runs alone it is 5.469e-6 (pair_reference_fine).
- My enlarged-box control gives 2.24e-13 and 5.10e-6.
- Energy is sampled at synchronized kick-drift-kick times (a half-kick is applied before each diagnostic), and all fixed-vacuum face bonds are counted.

**Scope or correction.** Charge conservation to roundoff follows from the scheme's structure: kick-drift-kick with a zero Dirichlet halo for u and w, so summation by parts makes it exact. It confirms the code is correct but says nothing about retention. It also means charge cannot exit the box, which is why the box control above is the real leakage test.

## Strongest agreement

Claim (2) holds, and on stronger grounds than the report gives: the old transverse amplification is core deflection, not instability.
- Rigid translation of the two outgoing cores reproduces the old code's canonical relative-shift amplification to within 1.3% at every grid: 8.947 against 9.064 at h=.25, 7.887 against 7.971, and 5.371 against 5.395.
- The old code's own core-displacement fits extrapolate to 15.3 and -7.8. The Cartesian continuum values are -14.39 and +7.66 (opposite sign convention).
- The direct Cartesian impact response is linear in the seed (0.02% and 0.06% change when it is halved) and linear in time (best exponential rate at most 6e-4 per tau, against 0.0138 implied).
- It converges at second order and balances momentum.

The old number is a deflection norm, and it is itself under-resolved (Richardson continuum about 11.2). It is not an internal-instability exponent.

## Strongest defect or narrowing

No DEFECT survived in the four dynamical claims. Every reported number reproduces exactly from the raw data.

The strongest narrowing is what 'pair durability' means. The two equal Q=399.4 objects merge completely at tau 40-55 (one lump with 99.99% of Q) and then re-split 1:2 into cores holding 32.5% and 65.0% of Q. The surviving cores are fission products, not the incoming objects.
- The smaller product (Q about 260, Omega about 0.848) has a rest-frame E/Q of 0.9595, a binding margin of about 4%. The stationary branch reaches E/Q=1 near Q of about 180.
- Only one collision geometry was tested (v=0.2, phase offset pi/2, equal charges).
- The single-object seeds inject only 0.18-0.19% of the energy, and 0.27% for the negative pair seed.
- The pair cores are still shedding charge at up to 0.24 points per 40 tau.
- The retention freeze is self-attested: protocol and raw data share one commit, 479f111.

My deciding control shows the walls cannot be faking retention: enlarging the box changes the worst case by only -0.153 points (96.869% to 96.717% at coarse resolution), and the enlarged box still holds 96.32% at tau 240.

## Key numbers

- finite-object runs passing R12 retention: 26/26 (out/independent_dynamics.txt section B)
- minimum joint late R12 retention (pairs): 96.98103% (pair_negative_fine) (out/independent_dynamics.txt)
- perturbed single late R12 retention, h=.75/.5/.375: 99.85098 / 99.84893 / 99.84842 % (out/independent_dynamics.txt)
- Richardson continuum pair retention, reference / generic: 97.517% / 97.103% (p=2.22 / 2.25) (out/independent_dynamics.txt section C)
- pair_reference max late Q12 difference, coarse-medium / medium-fine: 0.21427% / 0.06406% (p=2.44) (out/independent_dynamics.txt)
- box-control change, worst-case negative-seed pair (coarse, L24/Lz48 to L36/Lz60): 96.86943% to 96.71650% (-0.153 points) (out/box_control_compare.txt)
- negative-seed pair retention at tau 240 in the large box: 96.32335%, slope -0.183 points per 40 tau (out/box_control_compare.txt)
- closest approach of an R12 ball to a fixed-vacuum face (pair_negative_fine): 0.62 (out/extra_checks.txt)
- merged-lump charge during merger (tau 45-50): Q12=798.79 (99.99% of Q0), core separation 0.00 (out/extra_checks.txt)
- product core charges, fine reference: 259.87 (32.5%) / 519.20 (65.0%) (out/independent_dynamics.txt)
- smaller product rest-frame E/Q against the stationary branch: 0.9595 vs 0.95991 at Q=258.7, Omega=0.8483 (out/products_lab_EQ_and_shape.txt; out/products_branch.txt)
- injected energy of the single generic seed: +0.1858% of E (about 1.9% of binding energy 35.7) (out/extra_checks.txt)
- impact x/delta, h=.75/.5/.375: -13.726052/7.333882; -14.114130/7.523045; -14.241539/7.586161 (out/independent_dynamics.txt section F)
- impact x/delta continuum (Richardson): -14.388 / +7.661 (p=2.18 / 2.13) (out/independent_dynamics.txt)
- best exponential growth rate fitted to impact response, tau 80-160: at most 6e-4 per tau (0.0138 per tau implied by 9.064 over 160 tau) (out/independent_dynamics.txt)
- old m=1 amplification against the rigid-translation estimate: 9.064 vs 8.947; 7.971 vs 7.887; 5.395 vs 5.371 (out/rigid_norm_estimate.txt)
- old m=1 amplification, Richardson continuum: about 11.15 (p=1.46) (out/old_m1_vs_cartesian.txt)
- old m=1 core-displacement fits, Richardson continuum: 15.28 / -7.82 (out/old_m1_vs_cartesian.txt)
- Cartesian late R36 symmetry residual, pair_reference t160, c/m/f (translation): 0.2288 / 0.2264 / 0.2263 (resolution-independent; 0.187 in the larger box) (out/misc_checks.txt)
- fixed-Q l=0 lowest eigenvalue (dr=.025): 0.364348863 (out/independent_spectral.txt)
- l=2 lowest eigenvalue (dr=.025): 0.285483333 (out/independent_spectral.txt)
- dQ/dOmega at Omega=0.8: -4093.2 (author -4094.03) (out/independent_spectral.txt)
- maximum charge drift / Hamiltonian drift, 29 runs: 4.907e-13 / 2.762e-5 (object runs only 5.469e-6) (out/independent_dynamics.txt section A; out/misc_checks.txt)
- build validation against published pair_generic_coarse: maximum absolute difference 3.06e-11 (Bash run; ctrl/COMMANDS.txt)

## Scripts and outputs

- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/independent_dynamics.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/independent_dynamics.txt (conservation, retention, convergence/Richardson, controls, impact x/delta)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/old_m1_vs_cartesian.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/old_m1_vs_cartesian.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/rigid_norm_estimate.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/rigid_norm_estimate.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/extra_checks.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/extra_checks.txt (wall gaps, injected energy, components, merger phase)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/independent_spectral.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/independent_spectral.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/products_branch.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/products_branch.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/products_lab_EQ_and_shape.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/products_lab_EQ_and_shape.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/misc_checks.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/misc_checks.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/scripts/box_control_compare.py -> /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/box_control_compare.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/ctrl/COMMANDS.txt (exact build and run commands for the deciding control)
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/ctrl/control_pair_negative_coarse_bigbox_T240.csv
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/ctrl/repro_pair_generic_coarse_neg.csv
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/ctrl/repro_pair_generic_coarse.csv
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/out/SHA256SUMS.txt
- /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/issue13/dynamics/raw/ex/raw/ (verified extracted raw data)

## Not done

- I did not re-run any medium or fine production campaign. The only new evolutions are three coarse runs: a build-validation reproduction, a small-box baseline and the enlarged-box control.
- I did not run an absorbing-boundary (sponge) control; I used box enlargement instead.
- I did not test other collision velocities, phase offsets, unequal charges, repeated collisions or larger perturbations, so durability of a population of objects is untested.
- I did not evolve finite D or metric perturbations; everything is D=0 flat classical.
- The Hessian check re-implements the same finite-difference method. It is not method-independent and not interval-rigorous.
- I did not audit the old cylindrical solver (kernel_modes.cjs) itself. I used only its archived large-domain traces.
- I did not inspect the float32 spatial volumes beyond the endpoint component lists.
- I did not audit the production/condensate runs, the H0 exclusion or the empirical files; they belong to other review keys.
- The grut_reality_audit_20261020 directory was skimmed only: it covers 1D plane-wave full-field tests, not 3D object dynamics.
- I made no commits, pushes or issue comments, per this task's scope.
