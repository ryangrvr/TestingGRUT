# Track m1_vacuum: the M1 one-vacuum potential and the Q = 200000 object (Codex d4747cb)

1 November 2026 · Claude, issue #13 (BRIEF_M1.md) · same-session AI review, not outside replication. Codex's packet
`research/reality_framework/grut_matched_candidate_20261027/` was read at d4747cb through `git show` only
(CLAUDE_HANDOFF.md, REPORT.md, DERIVATION.md, DEFORMATION_PROTOCOL.md, CANDIDATE_M1.json, DEFORMATION_RESULTS.json,
STABILITY_RESULTS.json, COLD_MODE_RESULTS.json, VERIFICATION.json, PROFILE_M1.csv, code/deformation.py,
energy_stability.py, cold_modes.py); the bulk statement was read from ed73791 (BULK_RESULTS.json, bulk.py,
DERIVATION.md). **Codex's code was not run.** Every number below is from this track's own code
(`code/potential_m1.py`, `code/qball_m1.py`, `code/hessian_m1.py`, `code/coldfluid_m1.py`) unless marked [Q]
(quoted from Codex). The only Codex data file used is PROFILE_M1.csv, read once to difference it against my profile.
Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled. Nothing committed.

**Conventions (re-derived, agree with Codex's DERIVATION.md).** V_φ = λ/4(φ²−v²)² − εv(φ−v)³ with λ = 0.03,
ε = 0.002, v = 0.924089637 GeV; canonical complex X with V ⊃ g²φ²|X|² + λ_X|X|⁴ (g = 0.1, λ_X = 1e-3; this is
½g²φ²χ² + ¼λ_Xχ⁴ for |X| = χ/√2, the form in the brief). X = vF e^{−iωt}/√2, φ = vS, ρ = m_X r, m_X = gv =
0.09240896 GeV, Ω = ω/m_X; energy density in units g²v⁴ = 100 m_X⁴, U = (b/4)(S²−1)² − c(S−1)³ + ½S²F² + ¼aF⁴
with a = 0.1, b = 3, c = 0.2. Q = 2ω∫|X|²d³x gives Q̄ ≡ g²Q = Ω∫F²d³ρ (so Q = 200000 is Q̄ = 2000) and charge
density n/m_X³ = 100 ΩF². E/(Q m_X) = E_b/Q̄, mass = E_b v/g. These are exactly Codex's normalisations.

## 1. One-vacuum argument and the (φ, |X|) sector at fixed charge density [D] (`code/potential_m1.py`, 8 s)

| Claim (Codex) | Mine | Verdict |
|---|---|---|
| V'/v³ = (u−1)[λu² + (λ−3ε)u + 3ε], discriminant −0.000144 | factorisation exact (3.5e-16); (λ−3ε)² − 12λε = **−1.44e-4** | CONFIRMED |
| unique for ε > 0.00171573 | disc = 0 at ε = λ(1 ∓ 2√2/3) = **0.00171573** and 0.058284 | CONFIRMED |
| V(v) = V'(v) = 0, V''(v) = 2λv² | yes; V'''(v) = 6(λ−ε)v, cubic −6.67 % | CONFIRMED |
| V/v⁴ = (u−1)²[λ(u+1)²/4 − ε(u−1)] ≥ 0 for 0 < ε < 2λ | bracket discriminant ε(ε−2λ) = −1.16e-4 < 0; grid u ∈ [−50, 50]: min V = 0 at u = 1 only; V strictly decreasing for φ < v, increasing for φ > v | CONFIRMED |
| hidden tree bound: λ_Y,eff = 0.758516, large-Σ margin 2445.04 | min_x(λ_Xx⁴ − 2κyx³) = −27κ⁴y⁴/16λ_X³ (checked numerically, 1e-10) → **0.758516**; for s ≥ s_c = M_Y²/2μ = 62.854 H_i, V ≥ s²[(1−r²)²s_c² − μ²/λ_Y,eff], r² = h²/2s_c² = 0.2025, bracket **2445.04** | CONFIRMED |
| bulk (zero pressure): S = 0.15117211, F = 2.40647676, Ω = 0.77586471, n = 449.313 m_X³, "unique positive quartic root" [Q ed73791] | Ω² = min 2U/F² = min_S [S² + 2√(aU_φ)]; my own squaring gives the same quartic 2.1S⁴ + 3.76S³ + 2.864S² − 0.288S − 0.036 (one sign change → one positive root, Descartes); **S = 0.151172107845, F = 2.406476756778, Ω = 0.775864707460, n = 449.31336788 m_X³**, three-equation residuals ≤ 1e-16; a 164-start 2D search over S ∈ [−2.5, 1.5] finds no other minimum | CONFIRMED |

**Does the mixed term reopen a second minimum at large |X|?** No, in two senses. (i) At zero density the X-dependent
part ½S²F² + ¼aF⁴ is minimised at F = 0 for every S, so the static potential has the single minimum (S, F) =
(1, 0). (ii) At fixed charge density n̄ the homogeneous rotating states are the stationary points of e(S, F) =
n̄²/2F² + U(S, F). The φ equation U_φ'(S) + F²S = 0 has **exactly one root for every F > 0**, in (0, 1): F² =
b(1−S²) + 3c(S − 2 + 1/S) is strictly decreasing on (0, 1) (dF²/dS ≤ −4.66), and there is no negative root
because V_φ' < 0 for all φ < v (the one-vacuum result) and F²S < 0 for S < 0 (sign scan of U_φ' + F²S on S ∈
[−3, 3]: one root at every F from 0.1 to 10). The F equation n̄²/F⁴ = S(F)² + aF² then has exactly one solution
because G(F) = F⁴[S(F)² + aF²] is strictly increasing (checked on F ∈ (0, 15]); the 2×2 Hessian of e at fixed n̄
is positive at every sampled density (eigenvalues 0.94–14 in the dimensionless units). So **no second minimum, no
coexistence of two homogeneous charged phases at the same density**. What the tilt does change: the old φ = 0
charged interior (exact for F² ≥ b when ε = 0) becomes φ > 0 with S ≈ 3c/F² at large F; the homogeneous value at
the core density F = 2.631 is S = 0.1178 and the object centre has **S_c = 0.12577** (Codex 0.12577; ε = 0 gave
0.0060), so the interior "bag" energy density rises from U_φ(0) = 75.0 m_X⁴ to U_φ(S_c) = 86.0 m_X⁴ (83.8 in the
bulk). The only remaining instability is thermodynamic: along the family, Ω(n) decreases up to the spinodal
**n_sp = 233.2 m_X³ (F = 1.857, Ω_min = 0.6766)** and increases beyond it; the pressure P = nΩ − e is negative
for n < n_bulk = 449.3 m_X³ and zero there. That is the ordinary Q-ball phase-separation structure, now with
Ω_bulk = 0.77586 instead of the symmetric 0.74008 (Codex's "unique positive quartic root" is a statement about
the zero-pressure point, which I confirm; it is not a statement that e(n) is convex, and it is not).

## 2. Independent soliton and its stability [C] (`code/qball_m1.py`, 57 s; `code/hessian_m1.py`, 3.8 min)

**Method (own code, not Codex's).** Substitution u = ρF, w = ρ(S−1) removes the first-derivative terms:
u'' = (S² + aF² − Ω²)u, w'' = ρ[U_φ'(S) + F²S]; Numerov (4th-order) three-point discretisation on ρ_i = ih,
u = w = 0 at ρ = 0 and ρ = L; damped Newton with the exact sparse Jacobian; in fixed-Q̄ mode Ω is an extra
unknown with the extra equation Q̄ = 4πΩh Σu_i² (trapezoid, exponentially accurate here because the integrand
is even at 0 and vanishes at L). Energy from 4th-order differences, cross-checked by the virial identity
E = 2ΩQ̄ − 2E_U. The tanh initial guess converges in 6 Newton steps. Codex used scipy `solve_bvp` on
(F, F', S, S', ∫ρ²F²) seeded from the old symmetric profiles; nothing of that is reused.

**Primary Q = 200000** (L = 60, h = 0.05/0.025/0.0125, Richardson; convergence ratios 0.0625 = 1/16 as expected
for 4th order; L = 90 changes nothing at 1e-14; virial residual 2.9e-11):

| Quantity | Codex [Q] | Mine (extrapolated) | Difference |
|---|---|---|---|
| Ω = ω/m_X | 0.8817204919524214 | **0.8817204919524271** | 6e-15 |
| E/(Q m_X) | 0.9373999905262296 | **0.9373999905262296** | < 1e-16 |
| central S = φ_c/v | 0.1257717809745518 | **0.1257717809745547** | 3e-15 |
| central F | 2.6309605904405573 | 2.6309605904405690 | 1e-14 |
| mass (GeV) | 17324.832346421026 | **17324.832346421023** | −4e-12 GeV |
| profile F(ρ), S(ρ), 0 < ρ < 60 | PROFILE_M1.csv | max |ΔF| = 5.9e-11, max |ΔS| = 4.2e-11 | |

Verdict: **CONFIRMED** (mass 17324.83235 GeV, E/(Qm_X) = 0.93739999, Ω = 0.88172049, S = 0.12577; binding 6.260 %).
Radii: r₅₀ = 3.78, r₉₀ = 5.93, r₉₉ = 8.42 /m_X; E − ωQ = 11136 m_X = 1029 GeV (symmetric model: 1148 GeV).

**First law and dΩ/dQ.** Central differences at ΔQ = 500/1000/2000: dE/dQ / (Ω m_X) − 1 = 7.15e-8 / **2.87e-7** /
1.15e-6 — Codex's 2.9e-7 is reproduced at its own step ΔQ = 1000 and scales as (ΔQ)², i.e. it is the finite-difference
truncation, not a solver error; the first law holds to the solver's accuracy. dΩ/dQ = **−2.0500e-7** per unit charge
(Codex −2.0501e-7). Both CONFIRMED.

**Vakhitov–Kolokolov criterion** (Lee–Pang 1992; Friedberg–Lee–Sirlin 1976): dQ/dω < 0 ⇔ stable radial sector.
Here dQ/dΩ = 1/(−2.05e-7) = **−4.9e6 < 0** at the reference object; along my fixed-Ω branch (L = 100, h = 0.05,
Ω = 0.882 → 0.995) Q(Ω) decreases down to **Q_min = 5.67e4 at Ω = 0.9728** and rises again on the thick-wall branch, so
all objects with Q > 5.67e4 satisfy VK, including Q = 2e5. Consistent with Codex's positive fixed-Q Hessian, with
one correction to that statement:

**Fixed-charge radial Hessian** (own build, 3-point Laplacian like Codex's, on my re-solved profile; dense `eigh`):
L₊ alone has exactly one negative eigenvalue, **−0.09515 m_X²** (the usual single negative direction of a ground-state
Q-ball), and the charge-constraint rank-one term 4Ω²|ρF⟩⟨ρF|/∫ρ²F² (re-derived from δ²[Q̄²/2∫F²]) lifts it. The
lowest eigenvalue of L₊ + K is 0.23078, 0.22587, **0.22434**, **0.22332**, 0.22283, 0.22268 m_X² for box radii L = 40,
60, 80, 120, 200, 300 (h = 0.1; h = 0.05 at L = 80 gives the same 0.22434). The two bold values are Codex's 0.2243397
and 0.2233198 to five digits, so their build is reproduced; but the sequence converges to the **continuum edge
m_X² − ω² = 1 − Ω² = 0.22257** like 1 − Ω² + (1.04–1.3)(π/L)², and the eigenvector's weight inside ρ < 15 falls from
12 % (L = 40) to 0.0 % (L = 300): it is the lowest free-X box mode, not a mode of the object. **CORRECTED
(interpretation)**: the quoted "positive lowest eigenvalue 0.2233–0.2244" is box-dependent and carries no information
about the object; the correct statement is "no discrete radial fixed-charge mode below the continuum threshold
0.2226 m_X²", which is what VK requires and what I find. The L₋ phase zero mode (−2.77e-5 at h = 0.1, −6.9e-6 at
h = 0.05) and the ℓ = 1 translation zero mode (−1.25e-4, −3.1e-5) converge as h², numerically identical to Codex's:
CONFIRMED.

**E(Q)/Q and the approach to the bulk value** (fixed-Q continuation, h = 0.025, L = 60–110):

| Q | Ω | E/(Q m_X) | S_c | r₉₀ m_X | (E/Q − Ω_bulk) Q^{1/3} |
|---|---|---|---|---|---|
| 7.0e4 | 0.941722 | 1.001065 | 0.267 | 5.92 | 9.28 |
| 8.0e4 | 0.930994 | 0.992942 | 0.220 | 5.73 | 9.35 |
| 1e5 | 0.916103 | **0.978960** | 0.175 | 5.62 | 9.43 |
| 1.5e5 | 0.894402 | 0.953999 | 0.137 | 5.73 | 9.46 |
| 2e5 | 0.881720 | **0.937400** | 0.126 | 5.93 | 9.45 |
| 3e5 | 0.866569 | 0.916067 | 0.119 | 6.38 | 9.39 |
| 5e5 | 0.850935 | **0.892765** | 0.118 | 7.12 | 9.28 |
| 1e6 | 0.834297 | **0.867027** | 0.121 | 8.48 | 9.12 |
| 2e6 | 0.821559 | 0.846975 | 0.126 | 10.30 | 8.96 |
| 5e6 | 0.809039 | 0.827198 | 0.132 | 13.59 | 8.78 |
| 1e7 | 0.801973 | 0.816078 | 0.136 | 16.90 | 8.66 |
| 3e7 | 0.793787 | 0.803285 | 0.140 | 24.14 | 8.52 |

Codex's Q = 1e6 control (Ω 0.83429668, E/Qm 0.86702651, S 0.12122) [Q]: mine 0.83429668, 0.86702651, 0.12122 —
CONFIRMED. The ratio approaches Ω_bulk = 0.77586 from above with a surface term: a two-term fit to the last two rows
gives E/Q ≈ Ω_bulk + 8.2 Q^{−1/3} + 100 Q^{−2/3}, i.e. E/Q is still 0.803 at Q = 3e7 and reaches 0.78 only near
Q ~ 1e10. Codex's "E/Q → 0.77586 at infinite size; surface energy raises the reference object to 0.9374"
[Q ed73791]: CONFIRMED as a limit statement. **Binding threshold under M1: E/Q = m_X at Q = 7.13e4 (Ω = 0.9402)**;
objects exist down to Q_min = 5.67e4 (Ω = 0.9728), those with 5.7e4 < Q < 7.1e4 are stationary but heavier than
free X. Symmetric model (kin_thermo, ε = 0): 6.56e4 at Ω = 0.9312 and Q_min = 5.20e4 at 0.9675. The tilt therefore
raises every energy per charge: bulk 0.7401 → 0.7759, reference 0.9233 → 0.9374, threshold charge 6.6e4 → 7.1e4. [NEW]

## 3. Cold-fluid dispersion [D/E] (`code/coldfluid_m1.py`, 54 s)

Re-derived (NOTES convention X/v = (F + ψ₁ + iψ₂)e^{−iΩt}/√2, S = S₀ + s, modes e^{νt+ik·x}):
(ν² + k² + 2aF²)ψ₁ + 2Ωνψ₂ + 2FS s = 0, (ν² + k²)ψ₂ − 2Ωνψ₁ = 0, (ν² + k² + D)s + 2FSψ₁ = 0, D = U_φ''(S) + F²;
the determinant is Codex's cubic in z = ν² exactly. Own root finder on 30001 k values plus golden-section
refinement of the peak:

| supplied F | n (m_X³) | S, Ω | Γ_peak/m_X (mine) | k_peak/m_X (mine) | Codex [Q] | band k_max | verdict |
|---|---|---|---|---|---|---|---|
| 0.1 | 1.0 | 0.99833, 0.99883 | 0.0011667 | 0.04825 | 0.0011663, 0.0475 | 0.068 | CONFIRMED |
| 0.5 | 24.3 | 0.95762, 0.97059 | 0.029222 | 0.23465 | 0.029222, 0.235 | 0.334 | CONFIRMED |
| **1.0** | 88.0 | 0.82125, 0.88003 | **0.114280** | **0.41923** | 0.11428, 0.42 | 0.615 | CONFIRMED |
| 1.5 | 165.7 | 0.56356, 0.73661 | 0.183791 | 0.44656 | 0.18379, 0.4475 | 0.686 | CONFIRMED |
| 2.0 | 274.5 | 0.26631, 0.68624 | 0 | – | 0 | – | CONFIRMED |
| bulk 2.4065 | 449.3 | 0.15117, 0.77586 | 0 | – | – | – | NEW: bulk stable |
| core 2.631 | 581.6 | 0.11781, 0.84028 | 0 | – | 0 (F = 2.63) | – | CONFIRMED |

The growing branch is the k → 0 sound mode with z/k² → F²(2S² − aD)/[D(aF² + 2Ω²) − 2F²S²], which equals
**−c_s² = −n(dΩ/dn)/Ω** of the homogeneous equation of state (agreement 1e-4 at every density): instability ⇔
imaginary sound speed ⇔ dΩ/dn < 0, i.e. densities below the spinodal n_sp = 233.2 m_X³ (F < 1.8566) found in §1;
Codex's criterion aD < 2S² is the sign of the numerator (the denominator is positive throughout). **Meaning:** this is
the standard fragmentation of a dilute homogeneous charged condensate into dense objects (Kusenko–Shaposhnikov 1998;
Enqvist–McDonald 1998; Kasuya–Kawasaki 2000), with growth time 1/Γ = 8.8/m_X and fastest wavelength 2π/k = 15/m_X
at F = 1 (charge in one such cube 2.96e5, Codex 2.95e5 — a diagnostic, not a selected Q, as Codex also says). It is
not an instability of the bulk phase or of the finite objects: the bulk density, the object core density and every
n > 233 m_X³ have no growing mode, and §2 shows the finite objects satisfy VK. Codex's statement CONFIRMED with
this reading made explicit.

## 4. Walls and interfaces [A] with numbers

- **Energy of the former −v point:** V(−v) = 8εv⁴ = 0.016 v⁴ = **160 m_X⁴** (the brief's guess is right); V(0) =
  (λ/4 + ε)v⁴ = 95 m_X⁴ (was 75); V is monotone on φ < v with slope V'(−v) = −12εv³, inflection points at u = −0.385
  and 0.519. There is no barrier and no metastable point anywhere below v.
- **No static wall exists:** for a static 1D configuration with φ → v at infinity, ½φ'² = V(φ) ≥ 0, and a turning
  point needs V = 0, which happens only at φ = v; the only static solution is φ ≡ v. (The old ±v kink had tension
  σ₀ = (2√2/3)√λ v³ = 163.3 m_X³; its wall-network mechanism is gone, as Codex says.)
- **Transient negative-φ patches** (Codex's hot runs) are large-amplitude oscillations about +v with |δφ| > 2v. They
  **cost** energy to create — at least 160 m_X⁴ per unit volume of potential energy (1.9× the object interior's
  86 m_X⁴, 2.1× the φ = 0 barrier) supplied by hot kinetic/gradient energy — and **release** it on relaxation into φ
  oscillations (m_φ = √6 m_X) and, through g²φ²|X|², X pairs. Homogeneous roll from rest at −v to +v takes
  1.91/m_X = 4.67/m_φ; a patch of radius R has volume energy exceeding the old wall energy σ₀4πR² for R > 3σ₀/ΔV =
  3.1/m_X, so no patch is held up by surface tension. Implication: no wall formation, no wall collapse and no
  two-vacuum energy release; negative patches are a measure of local φ energy density > 0.016 v⁴, i.e. a hot-field
  diagnostic, and their relaxation is part of the ordinary φ-oscillation budget (whose decay Γ(φ→XX̄) = 1.9e-4 m_X,
  Codex's later packet, is slow). This does not exclude other defects (e.g. hidden-U(1) vortices in the gauged
  version), which Codex also does not claim.

## Limits

Radial, classical, tree-level, zero temperature; the ℓ ≥ 2 and non-spherical spectrum, gauge and quantum
corrections, finite-temperature profiles and the nonlinear fragmentation outcome are not computed (Codex claims none
of them). The hidden tree bound is checked at tree level with the heavy fields treated as in Codex's derivation;
nothing here bears on formation, the mass calibration chain or on whether the fields exist. No 3D runs were made.

## Reproduction

```
cd <scratch>; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; T=<track dir>
git show d4747cb:research/reality_framework/grut_matched_candidate_20261027/PROFILE_M1.csv > PROFILE_M1_codex.csv
python3 -I $T/code/potential_m1.py $T/results/potential.json                                   # 8 s
python3 -I $T/code/qball_m1.py PROFILE_M1_codex.csv $T/results/soliton.json primary_profile.npz  # 57 s
python3 -I $T/code/hessian_m1.py $T/code $T/results/hessian.json                                # 3.8 min
python3 -I $T/code/coldfluid_m1.py $T/results/coldfluid.json                                    # 54 s
```
Results: `results/potential.json`, `results/soliton.json`, `results/hessian.json`, `results/coldfluid.json`,
`results/claims.json` (verdict table).
