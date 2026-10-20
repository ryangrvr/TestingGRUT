# Independent re-implementation of the closed-whole demo (from specification only)

Date: 2026-10-20  
Python 3.13 / numpy 2.5.3 / scipy 1.18.1 / cvxpy 1.9.3 (CLARABEL)

## Code path

- Re-implementation (written before reading the reference script): `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/reimpl.py` (sha256 `a4740e24fb2c1efd300baf1b86b6f0b103611655b83c04059e99017bc89948ad`)
- Independent expm cross-check of the eigendecomposition evolution: `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/xcheck_expm.py`
- My outputs: `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/reimpl_results.json`, time series `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/series.npz`
- Comparison script and output: `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/compare.py`, `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl/comparison.json`
- Reference results read AFTER computing mine: `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/results/closed_whole_demo.json`, `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/results/marginal_compatibility.json`
- The reference script `closed_whole_demo.py` was NOT read at any point of this task (only its sha256 was taken for provenance).

Commands:
```
cd /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/reimpl
python3 reimpl.py          # ~5 s evolution loop + SDP
python3 xcheck_expm.py     # expm cross-check
python3 compare.py         # comparison + this report
```

## Method

- Qubit 0 is the leftmost tensor factor; operator on qubit k = I⊗…⊗P_k⊗…⊗I. H built term by term from the spec; asserted Hermitian and equal to 0.5 Z_0 + I⊗H_Abar(7q) + 0.5 X_0X_1.
- Bath state exp(−H_Abar)/Z from eigh of the 7-qubit H_Abar (asserted equal to scipy expm version to 1e-12).
- Exact evolution U(t)=V e^{−iwt} V† from eigh(H) (256×256); 9-qubit record protocol uses eigh of H⊗I_R (512×512); CNOT = P0_0⊗I + P1_0⊗X_8, asserted CNOT²=I.
- Partial trace by reshape+einsum; S_A = −Σλ ln λ; trace norm via SVD; grid t = 0.1·k, k=0..400.
- Werner compatibility: single SDP `maximize p s.t. rho_ABC ⪰ 0, Tr=1, Tr_C rho = W(p), Tr_A rho = W(p)`; solution verified (min eigenvalue −1.4e−9, marginal residuals 6e−15) and checked against the closed form p_max = 2/3 from the two-qubit symmetric-extension criterion Tr ρ_B² ≥ Tr ρ_AB² − 4√det ρ_AB (root found by brentq at 0.6666666666666667).

## My numbers
```
{
  "1_energy": {
    "E_tot_t0": -10.228788512441923,
    "max_abs_drift": 1.7763568394002505e-14
  },
  "2_entropy": {
    "S_A_t0": -6.661338147750941e-16,
    "S_A_t10": 0.39287118329981924,
    "S_A_t40": 0.41983791686476823,
    "S_A_max": 0.48933887859496117,
    "t_of_max": 37.5,
    "ln2": 0.6931471805599453
  },
  "3_BLP": {
    "D_t0": 1.0000000000000007,
    "sum_positive_increments": 1.2035373576160984,
    "D_min": 0.5418014611825231,
    "t_of_min": 37.6,
    "n_positive_increments": 199
  },
  "4_record": {
    "p_Z0_plus_t1_before_CNOT": 0.3549095082640668,
    "p_Z0_plus_t2_with_record": 0.3385311237554523,
    "p_Z0_plus_t2_no_intervention": 0.449147357617469,
    "p_ZR_plus_t2_with_record (record persists)": 0.35490950826406714,
    "Z0_t1_expectation": -0.2901809834718664
  },
  "5_initial_correlations": {
    "marginal_A_trace_norm_diff": 0.0,
    "marginal_Abar_trace_norm_diff": 0.0,
    "max_trace_distance_rhoA": 0.10248731291636753,
    "t_of_max": 0.7,
    "trace_distance_t0": 8.095563357601182e-19
  },
  "6_werner": {
    "p_max_sdp": 0.6666666737712065,
    "sdp_status": "optimal",
    "min_eig_solution": -1.4460216420750407e-09,
    "residual_AB_marginal": 6.453171330633722e-15,
    "residual_BC_marginal": 6.494804694057166e-15,
    "p_max_analytic_2_over_3": 0.6666666666666667
  }
}
```

expm cross-check (independent of eigh): S_A(10)=0.392871183300, S_A(37.5)=0.489338878595, S_A(40)=0.419837916865, D(10)=0.674532015003, D(37.5)=0.543930303492, D(40)=0.674858229496, TD(0.7)=0.10248731291636759, pZ0(3)=0.3549095082640673 — all equal to the eigh numbers to ≤1e-12.

## Comparison with the packet JSONs

| quantity | mine | packet | agree | note |
|---|---|---|---|---|
| 1. E_tot = Tr(H rho(0)) | -10.228788512441923 | -10.2287885124 | yes | rel diff 4.10e-12; packet value printed to 12 significant figures |
| 1. max |E_tot(t)-E_tot(0)| on grid | 1.7763568394002505e-14 | 7.1054273576e-15 | yes | abs diff 1.07e-14; both are double-precision round-off (<1e-13); a 1e-6 RELATIVE criterion is meaningless between two numbers that are themselves numerical zero, so compared as |a-b|<1e-12 |
| 1b. full 401-point E_tot(t) series: max |mine-packet| | 4.1941561335079314e-11 | 0.0 | yes | abs diff 4.19e-11; series compared pointwise (packet stores 12 sig figs) |
| 2b. full 401-point S_A(t) series: max |mine-packet| | 5.755951271169124e-13 | 0.0 | yes | abs diff 5.76e-13; series compared pointwise (packet stores 12 sig figs) |
| 2. S_A(t=10) | 0.39287118329981924 | 0.3928711833 | yes | rel diff 4.60e-13; also packet S5 S_A/10.0/a = 0.3928711833 |
| 2. S_A(t=40) | 0.41983791686476823 | 0.419837916865 | yes | rel diff 5.52e-13 |
| 2. max_t S_A | 0.48933887859496117 | 0.489338878595 | yes | rel diff 7.94e-14 |
| 2. argmax_t S_A | 37.5 | 37.5 | yes | abs diff 0.00e+00 |
| 3. BLP sum of positive increments of D | 1.2035373576160984 | 1.20353735762 | yes | rel diff 3.24e-12 |
| 3. min_t D(t) | 0.5418014611825231 | 0.541801461183 | yes | rel diff 8.80e-13 |
| 3. argmin_t D(t) | 37.6 | 37.6 | yes | abs diff 0.00e+00 |
| 3. D(t=0) | 1.0000000000000007 | 1.0 | yes | rel diff 6.66e-16 |
| 3. D(t=40) | 0.6748582294962724 | 0.674858229496 | yes | rel diff 4.04e-13 |
| 3. largest single increase of D | 0.01454521447272017 | 0.0145452144727 | yes | rel diff 1.39e-12 |
| 4. p(Z_0=+1) at t1=3 (before CNOT) | 0.3549095082640668 | 0.354909508264 | yes | rel diff 1.88e-13 |
| 4. p(Z_0=+1) at t2=6, record written | 0.3385311237554523 | 0.338531123755 | yes | rel diff 1.34e-12 |
| 4. p(Z_0=+1) at t2=6, no intervention | 0.449147357617469 | 0.449147357617 | yes | rel diff 1.04e-12 |
| 4. register population |0>_R after write (= p at t1) | 0.35490950826406714 | 0.354909508264 | yes | rel diff 1.89e-13 |
| 5. max marginal deviation rho1 vs rho2 (A and Abar) | 0.0 | 0.0 | yes | abs diff 0.00e+00 |
| 5. max_t trace distance rho_A(t): rho1 vs rho2 | 0.10248731291636753 | 0.102487312919 | yes | rel diff 2.57e-11 |
| 5. argmax_t of that trace distance | 0.7 | 0.7 | yes | abs diff 0.00e+00 |
| 6. max Werner p, AB & BC share B (SDP) | 0.6666666737712065 | 0.666666673771 | yes | rel diff 3.10e-13; both CLARABEL optima; exact value 2/3 (two-qubit symmetric-extension criterion), solver overshoot ~7e-9 |
| 6. max Werner p vs analytic 2/3 | 0.6666666737712065 | 0.6666666666666666 | yes | rel diff 1.07e-08; SDP optimum vs closed form; rel diff 1.07e-8 is solver tolerance |

**All 23 comparisons agree: True.** The packet stores values to 12 significant figures, so the attainable agreement on scalars is ~1e-12 relative, well inside the 1e-6 criterion. Two quantities are numerical zeros (energy drift 1.8e-14 vs 7.1e-15, marginal deviation 0 vs 0) and are compared absolutely; relative agreement is undefined for them, which is the only "explain" item.

## Provenance check of the packet

- sha256(code/closed_whole_demo.py) = `ee19642f84a4b943dffdda56b3ba1f1cbcd2302b7ce5f5428b4a44bcc180bdb1`; JSON claims `ee19642f84a4b943dffdda56b3ba1f1cbcd2302b7ce5f5428b4a44bcc180bdb1`; match = True
- sha256(code/marginal_compatibility.py) = `9a6333833d7b70ed10f729e253a1ce6d3a8b0c9a678850997a7876dffa600695`; JSON claims `9a6333833d7b70ed10f729e253a1ce6d3a8b0c9a678850997a7876dffa600695`; match = True

## What the numbers demonstrate (standard QM only, no new physics)

- Closure forces: E_tot constant to 1e-14 on a mixed global state (unitarity). S_tot/purity constancy and S_A=S_Abar for pure states are the same theorem.
- Closure does not force monotone relaxation: S_A rises to 0.489 at t=37.5, sits at 0.420 at t=40, with 13 local maxima; never reaches ln2=0.693.
- Closure does not force Markovianity: BLP backflow Σ(ΔD)+ = 1.2035 > 0 (199 of 400 grid steps increase D); a CP-divisible reduced dynamics would give 0.
- Records: writing a CNOT record at t=3 changes A's later statistics (0.3385 vs 0.4491 at t=6); the register keeps p=0.35491 exactly because [H⊗I_R, Z_R]=0, a property of the supplied H, not of closure.
- Initial correlations: rho1=I/256 and rho2=I/256+Y_0X_1/512 have identical A and Abar marginals (deviation exactly 0) but rho_A(t) differ by trace distance 0.1025 at t=0.7; no map on the marginals reproduces both futures (Pechukas).
- Marginal compatibility: Werner AB and BC sharing B with identical overlaps admit a common 3-qubit state iff p ≤ 2/3 (SDP 0.66666667, closed form 2/3); above that, overlap agreement is necessary but not sufficient.
