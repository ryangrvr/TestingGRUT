# Operators track: numerical-analysis review of ff49e5d (causal CFS receiver)

Claude, 25 October 2026. Scope: CLAUDE_HANDOFF items 1, 2, 5 and 6 of
`research/reality_framework/grut_causal_receiver_20261024/` at ff49e5d. I read the packet only through
`git archive` into the scratchpad. Everything here is numerical analysis of one supplied-field working model.
Nothing here is evidence that the fields exist, and nothing here validates GRUT.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued. A number is **independent** (my code,
`code/*.py`) unless marked **REPRODUCTION** (Codex code run by me).

## Verdicts

| # | Codex claim | Codex value | Mine | Verdict |
|---|---|---|---|---|
| 1a | The direct CFS realization is equivalent to the stretched operator, starting from zero memory | equations in MATHEMATICS §1 | Rederived. The zero-memory system gives s²U−sU0−V0 = C_n D C_f G U − m²U exactly. Max relative difference 8.6e-14 (3 grids × 4 complex s) [D,C] | CONFIRMED |
| 1b | The C++ right-hand side implements those equations | — | The C++ matrix and my independently assembled matrix differ by at most 4.4e-16 (n=64) and 3.6e-15 (n=128) (REPRODUCTION export) [C] | CONFIRMED |
| 1c | Constant-coefficient Routh–Hurwitz analysis | Δ2=2(a³+σk²), Δ3=4σk²(am²+αa²+αk²) | Checked symbolically with sympy [D] | CONFIRMED |
| 1d | Initial charge fraction inside the layer | ≈1.63e-15 | 1.6315e-15 [C] | CONFIRMED |
| 1e | No growing linear modes | none at n=64 and 128 | Production geometry (R=360/L=450 at h=0.5 and 0.25; R=720/L=810 at h=0.5; both masses, state dimension up to 4317): largest Re λ is −2.3e-8 to −1.9e-7. RK4 amplification is below 1 at dt=0.03125 and at dt=0.015625 [C] | CONFIRMED (linear and spectral only) |
| 1f | Maximum sampled reflection | 1.82e-8 charged, 6.516e-4 neutral | 1.8222e-8 charged, 6.5161e-4 neutral, from a different method: a point source plus a two-wave fit, with no lead roots. Relative agreement ≤5.8e-7 wherever \|R\|² > 1e-14 [C] | CONFIRMED (for k ≥ 0.01) |
| 2a | G and D form a summation-by-parts pair with the odd/even closure; L = −GᵀG | stated | D+Gᵀ = 0 exactly; ‖L7−DG‖·h² = 4e-16; L is symmetric negative definite; the semi-discrete energy rate is 8.8e-16 relative [C] | CONFIRMED |
| 2b | Fourth-order accuracy | stated | Observed order 3.998 for the gradient and the Laplacian on regular U=rX, right up to r=0. Dispersion symbol matches to 1.3e-14, with leading error k−p_h = (3/640)k⁵h⁴ [D,C] | CONFIRMED, with a qualification (below) |
| 2c | Lead matching: 2 oscillatory roots + 2 left-decaying evanescent roots | stated | For h from 0.0625 to 2 and k from 1e-4 to 3, there are always 2 roots on the unit circle, 2 inside and 2 outside. The outside pair is complex with \|z\| ≥ 25.96, the limit of z²−26z+1=0. The selection rule \|z\|>1.01 is therefore safe. With no layer, \|R\|²=1 to within 1e-15 [C] | CONFIRMED |
| 5a | The extraction currents are exact discrete conservation laws | stated | Derived for any node-local potential and any U(1)-invariant force [D]. Residuals of my own implementation are ≤2.2e-14 relative (complex-step derivative, random states; cut and global ledgers) [C] | CONFIRMED for the semi-discrete scheme |
| 5b | Closure "with measured RK4 errors" | I025: 3.3e-7 (E), 6.28e-8 (Q) | In my run the fully discrete closure is RK4 truncation error. It falls 1.2e-4 → 3.9e-6 → 1.2e-7 → 4.1e-9 as dt halves from 0.1 to 0.0125 (observed order ≈4.9) [C]. I did not rerun I025 itself | CONFIRMED (wording: "exact" means semi-discrete) |
| 5c | Instrumented fields are bitwise identical to the plain run | stated | Final fields, memory, trace and history columns 0–7 are byte-identical on a short nonlinear run (REPRODUCTION) [C] | CONFIRMED |
| 6a | Bound-mode snapshot numbers | P025_time: one mode, Q=11182.655, ω=0.864632 | My own solver and my own re-implementation of the split give 11182.6551 and ω=0.8646321. My final Q60 fraction is 0.35020722689 (Codex 0.35020722691); my late mean is 0.351876650083 (Codex 0.351876650103) [C] | CONFIRMED |
| 6b | The split is a frozen single-field diagnostic, not the coupled fluctuation problem and not a subtraction | stated | Accurate as a description. The operator is the mean-field operator −L+S²+0.1\|X\|², not the charged fluctuation Hessian (see below) [D] | CONFIRMED (wording clarified) |
| 6c | (Implicit use) the snapshot split as a bound/free charge separation | Table: P025_time one mode/11183; P0125 two modes/13324 | The split oscillates in time: 11063–13954 in the late window while Q60 changes by only 1.1%. Errors reach 20–36% of the core charge in strongly excited cores (below) [C] | REFUTED as a separation of excited cores; usable only for relaxed, near-stationary cores |

## 1. Causal CFS realization (item 1)

**Derivation [D].** Stretch ∂_r → Γ⁻¹∂_r with Γ = 1+σ/(s+α), which gives C = Γ⁻¹ = (s+α)/(s+α+σ) = 1−σ/(s+α+σ).
Write C_f G U = G U + ψ with ψ̂ = −σ_f/(s+α+σ_f) G Û. This is ψ̇ = −(α+σ_f)ψ − σ_f G U, and it requires ψ(0)=0.
Similarly C_n W = W + q with W = LU + Dψ, so q̇ = −(α+σ_n)q − σ_n W with q(0)=0. Laplace-transforming the full system
gives s²Û − sU0 − V0 = C_n D C_f G Û − m²Û, for any initial U0 and V0, including data inside the layer. This matches
MATHEMATICS §1 and `radial_cfs.cpp` term by term. Zero memory is therefore the causal choice: initial data inside the
layer acts as a source in the stretched medium and is not set to zero. The prepared tail in the layer is 1.63e-15 of
the charge. The mass and the quartic force are local and are correctly left unstretched. That is the standard
nonlinear PML convention. It is reflectionless only for the linear part, as Codex states. The method follows
Berenger, Kuzuoglu–Mittra (CFS) and Roden–Gedney (convolutional/auxiliary realization).

**Reflection.** Mine is computed on the actual radial grid: a point source near r=8 is solved in the frequency
domain, and A e^{ikhj} + B e^{−ikhj} is fitted on 40 ≤ r ≤ 300. This method is independent of Codex's lead-root
matching. It agrees with Codex to 6e-7. At the band edge k=0.01, both equal the ideal continuum
round trip to the wall, exp(−σ_max W p ω/(α²+ω²)). That limit comes from the layer's thickness, not from
discretization. Values below about 1e-14 are solver noise in both codes.

**Scope limit [C,A].** Below the sampled band the reflection rises:

| k | Charged \|R\|² | Neutral \|R\|² |
|---|---|---|
| 0.005 | 1.35e-4 | 2.55e-2 |
| 0.002 | 2.83e-2 | 0.231 |

Such waves (v < 0.013) cannot reach R=360 from r ≤ 200 by τ=12000, let alone return. They therefore cannot
contaminate these runs [A]. The "< 2e-8" figure is a statement about k ≥ 0.01 only.

## 2. Staggered operators (item 2)

The SBP identity is exact (table). **Qualification [C].** The odd/even continuation has no closure error only when
the even derivatives vanish at the boundary. That holds for a regular U = rX at r=0, and at a Dirichlet wall with
U_rr = 0. If U_rr ≠ 0 at a boundary, the Laplacian error at the last three nodes is O(1): 0.17–0.33, and it does not
converge. The gradient error there is O(h). At the outer wall U_rr ≠ 0 generically, because the cubic σ is clipped
there, so C′ ≠ 0. This sits inside the absorber, and its effect is already included in the reflection numbers. It is
not a correction to any reported value.

## 3. Extraction currents (item 5)

**Derivation [D].** With V̇ = LU − F(U) + B, B = 0 near the cut, L symmetric and F node-local:

    d/dt Σ_{i≤c} e_i = −½ Σ_{i≤c<j} L_ij (U_i V_j − V_i U_j)

This holds for any local potential, including the quartic and cross-field terms. The charge version needs only
x_i F_y,i − y_i F_x,i = 0, which is U(1) invariance. Both match the `radial_flux.cpp` coefficients, 2π/h·c_ℓ and
4π/h·c_ℓ.

The energy split is the symmetric half-and-half node assignment. Another split changes J_E by the time derivative of
a bounded term at the cut, so the integrated flux is defined up to the field energy near the cut.

In the fully discrete scheme, quadratic and quartic invariants are not preserved by RK4. The closure therefore
measures time-step error. In the short reproduction at τ=600, the absolute closure (≈0.01) still equals the flux
(0.0103). The closure becomes informative only once the integrated flux is much larger than 1e-7 E0, as it is by
τ=12000.

## 4. Bound-mode diagnostic (item 6)

**Operator [D].** `bound_modes.py` diagonalizes H = −L + S² + 0.1|X|², identical for x and y. A stationary soliton
f e^{−iω0t} is an exact eigenvector of H with eigenvalue ω0², so the split is exact for it. My test gives a split
error of 1e-11 for my soliton (ω0=0.80701, Q=10435.3).

The linear fluctuation problem is different. In the co-rotating frame the a- and b-channels have mass terms
s²+0.3f² and s²+0.1f², they are coupled to δS through 2sf, and there is a gyroscopic term ±2ω0∂_t. For my soliton
that problem has:

- the expected zero modes;
- one internal mode Ω = 0.18896, just below the charged continuum at 1−ω0 = 0.19299, with lab-frame sidebands 0.618
  and 0.996;
- no growing mode (roundoff level, 4.7e-8) [C].

The frozen operator's second level, 0.9537, matches none of these frequencies. Its eigenvalues are not fluctuation
frequencies.

**Production time series [C].** My rerun of P025_time reproduces Codex to about 1e-10. Over 9600 ≤ τ ≤ 12000
(241 snapshots, every 10):

| Quantity | Value |
|---|---|
| Split bound charge Qb | 11063–13954 (mean 12869) |
| Q60 | 13647–13800; falls by only 152 over the window (dQ60/dτ = −0.062) |
| Qb/Q60 | 0.80–1.02 |
| Charge inside r<60 labelled continuum | −258 to 2731 |
| Number of subthreshold modes | 1 in 23% of snapshots, 2 in 49%, ≥3 in 28% |
| Dominant period | ≈35. Sampled every 10, so aliasing is not excluded. It is close to 2π/Ω_internal |

Codex's τ=12000 snapshot (one mode, 11183) lies near the minimum of this oscillation. The P025_time-versus-P0125
difference in REPORT (11183 versus 13324) is consistent with snapshot phase rather than resolution [E].

**Counterexamples (own soliton, own solver) [C].** The truth is defined as the charge that has escaped r<30 by the end
of the run.

- **Neutral S-velocity kick of size ε (adds no charge), ε = 0.02, 0.05, 0.1.**
  - The split's transient "free" charge peaks at 2.06×, 2.15× and 2.22× the charge that actually escapes (0.0234,
    0.143, 0.567).
  - It dips negative, to −0.83× the escaped charge.
  - Both quantities scale as ε², so near stationarity the error is of the same order as the free charge itself.
  - Once the radiation has left, the residual error is ≤2.5% of the emitted charge.
- **Soliton × 1.3, a strongly breathing core.**
  - Over 1000 ≤ t ≤ 3000 the split labels up to 3704 as free, but at most 163 escapes by t=4000. The core loses only
    0.028 per unit time.
  - The error reaches +3551 (20% of the core), with rms 964.
  - The number of modes flips between 1, 2 and 3.
- **Collapsing Gaussian.** The error ranges from −367 to +10 on a core of 1025, i.e. up to 36%, with Qb/Q30 up to 1.25.
- **Soliton × 0.75, relaxed.** The error is ≤ 2.4 (0.05%).

**Bound.** No small-parameter bound exists away from stationarity. Near stationarity the error is O(ε²), with a
coefficient of about 2 in these tests [C,E]. A defensible separation needs a time-averaged or Floquet projector over
the breathing period, or an asymptotic escape definition, or the coupled linear problem with the Klein–Gordon
symplectic norm [A].

## Uncertified / not checked

- Nonlinear matching inside the layer.
- Stability of the layer for nonlinear states, and at grids other than the five I computed.
- I025's own closure numbers (I did not rerun them), and formation retention apart from my P025_time rerun.
- Any asymptotic bound/free or particle spectrum.
- Hot-core leakage, cosmological formation and vacuum density.
- Credit: Q-ball formation and stability literature (Kusenko–Shaposhnikov; Enqvist–McDonald; Kasuya–Kawasaki);
  the Friedberg–Lee–Sirlin-type two-field structure is recalled [R].

## Reproduce

    cd research/reality_framework/grut_formation_review_20261025/tracks/operators/code
    python3 check_operators.py             # SBP, order, dispersion, roots        -> results/operators_sbp_order_roots.json
    python3 check_routh.py                 # symbolic RH, memory elimination      -> results/routh_and_memory.json
    python3 check_cfs.py [--codex-binary B]  # resolvent, (REPRODUCTION matrix), reflection, spectra -> results/cfs_checks.json
    python3 check_flux.py                  # discrete currents, RK4 closure       -> results/flux_checks.json
    python3 check_bound.py                 # soliton, neutral kicks, coupled spectrum -> results/bound_split_counterexample.json
    python3 check_bound_strong.py          # strongly excited cores               -> results/bound_split_strong.json
    python3 formation_bound.py             # own P025_time rerun + split time series (~30 min) -> results/formation_bound_P025_time_own.json
    sh repro_bitwise.sh                    # REPRODUCTION: Codex flux vs base binaries -> results/REPRODUCTION_bitwise_I025_vs_base.json

Here B is `operator_controls.cpp` from ff49e5d, compiled in the scratchpad. `results/reflection_comparison.json`
compares my reflection values with Codex's, row by row.
