# Skeptic review of track m1_vacuum (Codex d4747cb / ed73791, M1 one-vacuum potential and the Q = 200000 object)

1 November 2026 · Claude (adversarial reviewer, key `m1_vacuum_skeptic`), issue #13 · same-session AI review, not an
outside replication. Target: `tracks/m1_vacuum/NOTES.md` and `results/claims.json` (16 claims M1V-1 … M1V-16).
Codex's packet was read through `git show d4747cb:…` and `git show ed73791:…` only; the hidden-sector potential form
was taken from our own `tracks/rp_reheating/code/rp.py`. **Nothing of Codex's code and nothing of the track's code was
run**; every number below is from the skeptic's own code in `skeptic/code/` (five scripts, 4 min total on one core),
chosen to differ in method from both Codex and the track. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate,
[A] argued, [R] recalled. Nothing committed. No validation, exclusion or new theory is claimed.

## Verdict on the track

**No claim of the track is refuted.** All sixteen are reproduced by different methods to the digits quoted (the
soliton to 3e-14 in Ω and 6e-14 in E/Q on a second-order cell-centred discretisation after two-level Richardson; the
algebra as exact rationals; the dispersion cubic as the 6×6 characteristic polynomial; the fixed-charge Hessian as a
direct second variation of a discrete energy with a secular-equation test). Two small corrections to the track's own
numbers/wording and a handful of sharpenings are listed under *New issues*. The track's one CORRECTED verdict (M1V-10,
box mode) is right in substance but slightly overstated as a correction of Codex (see N2).

## Method differences (independence)

| Item | Codex | Track m1_vacuum | This skeptic |
|---|---|---|---|
| algebra (potential, bulk quartic, bound) | float numpy | float numpy + grid scans | sympy exact rationals, Sturm root count, mpmath 30 digits (`sk_algebra.py`) |
| soliton | scipy `solve_bvp` on (F,F',S,S',∫ρ²F²), Ω parameter | Numerov 4th order on u=ρF, w=ρ(S−1), fixed-Q Newton | cell-centred conservative 2nd-order FD on (F, S−1), fixed-Ω Newton, Illinois on Ω for the charge, two-level Richardson; DOP853 shooting as a divergence-radius diagnostic (`sk_qball.py`) |
| Hessian | 3-point Laplacian on Dirichlet nodes, dense eigh of L₊+K | same operator, own profile, box scan | second variation of the discrete cell-centred energy as a generalised eigenproblem A v = λ M v; K validated by second differences of E_Q; secular identity f(0) = (1/I) dQ̄/dΩ; free-box comparison (`sk_hessian.py`) |
| cold fluid | cubic roots | cubic roots | eigenvalues of the 6×6 first-order matrix M(k) (`sk_coldfluid.py`); symbolic identity σ_lw ≡ −c_s² |
| surface term | – | two-point fit of E(Q)/Q | planar wall BVP at Ω_bulk, σ direct (`sk_wall.py`) [NEW] |

## Claim-by-claim (Codex value · track value · skeptic value · refuted?)

| id | label | claim | Codex [Q] | track | skeptic (own code) | refuted |
|---|---|---|---|---|---|---|
| M1V-1 | D | V′ factorisation, discriminant, threshold | −0.000144; 0.001715728752538096 | −1.44e-4; (0.00171573, 0.058284) | exact: V′/v³ = 3(u−1)(5u²+4u+1)/500, disc = **−9/62500**, thresholds **3/100 ∓ √2/50** = 0.001715728752538099, 0.058284271247461904 | no |
| M1V-2 | D | V(v)=V′(v)=0, V″(v)=2λv², cubic −6.67 % | 6.67 % | 6.67 % | exact: V‴(v)/v = 21/125 vs 9/50, reduction **1/15** | no |
| M1V-3 | D | V ≥ 0, bracket positive | bracket > 0 for 0<ε<2λ | disc −1.16e-4 | exact bracket 3u²/400 + 13u/1000 + 19/2000, disc **−29/250000**, leading coeff 3/400 > 0; V(−v) = 2/125 v⁴, V(0) = 19/2000 v⁴ | no |
| M1V-4 | D | hidden tree bound λ_Y,eff, margin | 0.758516376953125; 2445.0408755072076 | same | mpmath: **0.75851637695312502**; **2445.0408755072068**; ratios M_Y/H = 30.000, h/H = 40.000, μ/H = 7.1594562 from CANDIDATE_M1; brute-force Nelder–Mead on the actual rp.py potential at fixed |Σ|: V_min(s) ≥ s²·margin at s = s_c (9.93e6 vs 9.66e6), 70, 100 (8.45e7 vs 2.45e7), 200, 500; 2e5 random points: V_min = 4.5e4 > 0 | no |
| M1V-5 | D | bulk point, unique positive quartic root | 0.1511721078446304, 2.406476756777611, 0.7758647074595397, 449.3133678845478 | identical | quartic exact **21/10, 94/25, 358/125, −36/125, −9/250**; Sturm: exactly **1** positive real root (not just Descartes); S = **0.15117210784463039228**, Ω = **0.77586470745953986781**, n = **449.313367884548**; unsquared residual 5e-32; the other real root −0.0747 fails the unsquared equation (residual −0.297) | no |
| M1V-6 | D | no second minimum at fixed n; spinodal 233.2 | unique S(F) ∈ (0,1] | spinodal F=1.857, n=233.2, Ω_min 0.6766 | aD = 2S² solved at 30 digits: S = 0.336233013834, F = **1.85662560903**, Ω = **0.676578528186**, n = **233.220586942**; same point is the Ω minimum of the family | no |
| M1V-7 | C | Q=2e5 object | Ω 0.8817204919524214, E/Qm 0.9373999905262296, S₀ 0.1257717809745518, 17324.832346421026 GeV | 0.8817204919524271, 0.9373999905262296, 0.1257717809745547, 17324.832346421023 | **Ω = 0.8817204919524464 (Δ = +2.5e-14), E/Qm = 0.9373999905262875 (Δ = +5.8e-14), mass = 17324.832346405 GeV (rel −9e-13)**, F₀/S₀ to 4e-7/2e-6 rel (limited by my O(h²) centre extrapolation; error ratios 3.5); Ω/E/mass error ratios 4.000; L = 80 vs 60: 2e-13. DOP853 shooting from Codex's (F₀,S₀,Ω): divergence radius 16.6 vs 11.1 for S₀+1e-9 and 8.1 for S₀+1e-6 → Codex's centre values are good to < 1e-9 | no |
| M1V-8 | C | dE/dQ = ω to 2.9e-7; dΩ/dQ < 0 | 2.870671e-7; −2.0500727e-7 | 7.15e-8/2.87e-7/1.15e-6 | fixed-Q neighbours on my grid: **7.09e-8 / 2.871e-7 / 1.148e-6** at ΔQ = 500/1000/2000 (h = 0.02; 6.8e-8/2.85e-7/1.148e-6 at h = 0.01), pure (ΔQ)²; fixed-Ω branch: −3.4e-6 at δΩ = 5e-4, −1.37e-5 at 1e-3 (pure δΩ²); **dΩ/dQ = −2.05006e-7** | no |
| M1V-9 | C | VK dQ/dω < 0; Q_min | not stated | −4.9e6; Q_min 5.67e4 at 0.9728 | Q_min = **56713 at Ω = 0.9727**; secular identity f(0) = 1 + c wᵀA⁻¹w = **−21.5070** vs (1/I)dQ̄/dΩ = **−21.5073** (identity verified to 1e-5): the lowest fixed-charge eigenvalue is positive iff dQ̄/dΩ < 0 | no |
| M1V-10 | C | fixed-Q radial Hessian lowest 0.2233–0.2244; track: box mode | 0.2243399 (L=80), 0.2233198 (L=120) | 0.23078/0.22434/0.22332/0.22283 (L=40/80/120/200) | my discretisation: **0.23081/0.22436/0.22334/0.22285**; FREE operator −d²+1−Ω² on the same grids: **0.22876/0.22413/0.22327/0.22283** (gap 2e-3 → 2e-5); eigenvalues below the edge 0.22259: **0 in every box**; L₊ negative count **exactly 1 (−0.0950)**; mode weight inside ρ<15: 0.124 → 0.001; Hessian vs second differences of E_Q: 1e-10. Track's reading confirmed; see N2 on the label | no |
| M1V-11 | C | zero modes → 0 as h² | −1.25e-4 / −2.77e-5 at h=0.1 | same | L₋ zero mode is exact in my discretisation (|λ| ≤ 3e-13: F solves the same stencil); ℓ=1: +2.2e-5 (h=0.1), +5.4e-6 (h=0.05), ratio 4 | no |
| M1V-12 | C | E/Q → 0.77586; Q=1e6 control; fit 8.2 Q^{-1/3}+100 Q^{-2/3} | 0.8342966820, 0.8670265120, S₀ 0.12122 | 0.8342967, 0.8670265; (E/Q−Ω_b)Q^{1/3} 9.46→8.52 | Q=1e6: **0.83429668205, 0.86702651202, 80121.02151 GeV**; Q=1e5 **0.97895985**, 5e5 **0.89276485**, 7e4 **1.00106477**; planar wall σ = **99.037 m_X³** (Richardson ratio 4.00) → asymptotic coefficient **8.164** (track fit 8.2 ✓); residual × Q^{2/3} = 75, 95, 108, 111 at Q = 2e5, 1e6, 1e7, 3e7 (track's "100" is a mid value, see N3) | no |
| M1V-13 | C | E/Q = m_X at Q = 7.13e4; Q_min 5.67e4 | – | 71263 (0.9402); 56710 (0.9728) | **71263 at Ω = 0.94017; 56713 at 0.97266** | no |
| M1V-14 | D | cold-fluid peaks | F=1: 0.11428 at 0.42; F=0.5: 0.029222 at 0.235; F=1.5: 0.18379 at 0.4475; F=2, 2.63 stable | 0.114280 at 0.41923 (edge 0.615); … | 6×6 eigenvalues: F=1 **0.1142796 at k=0.419229, edge 0.61526**; F=0.5 **0.0292219 at 0.234653**; F=1.5 **0.1837907 at 0.446562**; F=0.1 0.0011668 at 0.048252; F = 1.8566, 2.0, 2.4065, 2.63, 2.631: max Re ν ≤ 1e-14 for **0 < k < 20** (not only the scanned band) | no |
| M1V-15 | D/E | growth = −c_s² below the spinodal | aD < 2S² criterion | z/k² = −c_s² (1e-4) | **symbolic**: det(M−ν) ≡ Codex's cubic (difference 0); σ_lw + n(dΩ/dn)/Ω ≡ **0** along the family (sympy simplify); numerically z/k² = −c_s² to 3e-6 at every density; charge per fastest cube 2.963e5 | no |
| M1V-16 | A | no static walls; patches cost 160 m_X⁴; roll time 1.91 | – | 160, 95, σ₀ 163.3, roll 1.9059 | V(−v) = 160, V(0) = 95 (exact 2/125, 19/2000 v⁴); σ₀ = 163.299316 (closed form 4√(λ/2)v³/3); roll time by tanh-sinh quadrature **1.9079/m_X** (track 1.9059: trapezoid misses the inverse-square-root endpoint singularity; 0.1 %, no consequence) | no (N1) |

Units audited and re-derived: m_X = gv = 92.408964 MeV, m_φ = √(2λ)v = 226.35481 MeV (Codex 226.355), Q = 2ω∫|X|² ⇒
Q̄ = g²Q = Ω∫F²d³ρ, n/m_X³ = 100ΩF², E/(Qm_X) = E_b/Q̄, mass = E_b v/g (my E_b → 17324.8323464 GeV), ε-term sign (ε > 0
raises the φ < v side: V(0) 75 → 95 m_X⁴, V(−v) 0 → 160 m_X⁴), H_i ratios of the heavy point (30, 40, 7.1594562), the
hidden potential λ_X|X|⁴ + M_Y²|Y|² + λ_Y|Y|⁴ + λ_Σ(|Σ|²−h²/2)² + 2Re[κY*X³ + μΣ*Y²] (rp.py) — all consistent with the
track's and Codex's conventions. Surface tension unit: g²v⁴/m_X² = 100 m_X³ (I first dropped this factor myself; the
check (36π)^{1/3}σ n_b^{−2/3} vs the fitted 8.2 caught it — anyone quoting σ should carry it).

## New issues (none changes a verdict)

- **N1 (track number, trivial).** Roll time from rest at −v to +v: 1.9079/m_X = 4.673/m_φ, not 1.9059 (0.1 %;
  trapezoid on the 1/√(φ+v) endpoint). Nothing depends on it.
- **N2 (track label).** M1V-10 "CORRECTED": the box-mode identification is right and useful (my independent operator
  gives the same L-sequence and the same free-box levels to 2e-5 at L = 200), but Codex's DERIVATION.md already says
  "Positive finite-grid eigenvalues are numerical evidence, not continuum coercivity" and STABILITY_RESULTS.json's scope
  says the same; the correction is to REPORT.md/CLAUDE_HANDOFF.md's wording ("positive lowest fixed-charge eigenvalues
  0.2233–0.2244"). A fairer label is "CONFIRMED (numbers), interpretation added". Worth adding: Codex's own table shows
  the L-dependence (0.22434 at L = 80 vs 0.22332 at L = 120, i.e. the (π/L)² shift) without comment.
- **N3 (track fit).** The two-term fit E/Q ≈ Ω_b + 8.2Q^{−1/3} + 100Q^{−2/3} is confirmed in its leading term by the
  direct surface tension (8.164 from σ = 99.04 m_X³, n_b = 449.31 m_X³), but the second coefficient is not a constant:
  (E/Q − Ω_b − 8.164Q^{−1/3})Q^{2/3} = 75, 95, 108, 111 at Q = 2e5, 1e6, 1e7, 3e7. The conclusion "0.78 only near
  Q ~ 1e10" stands (8e9 with 8.164/111).
- **N4 (track wording, §1).** "The only remaining instability is thermodynamic" is correct for the homogeneous family
  and k ≠ 0 (my σ_lw ≡ −c_s² proof); the track's own limits already exclude ℓ ≥ 2, gauge and quantum sectors — the
  sentence should be read inside those limits.
- **N5 (Codex grid, for the record).** The deformation grid straddles the one-vacuum threshold: at ε = 0.001 the
  quadratic has real roots u = −0.770 (local minimum, V = 0.00679v⁴ = 67.9 m_X⁴, metastable) and −0.130 (maximum); the
  ε = 0 point has the old two vacua. Codex states that the primary ε = 0.002 was chosen above the threshold; the grid
  objects at ε ≤ 0.0017 live in a two-minimum potential, which the packet does not say explicitly. No consequence for
  the primary.
- **N6 (hidden bound, limit restated).** The tree bound consumes the whole light quartic λ_X|X|⁴ to absorb the κY*X³
  term and uses λ_Y = λ_Σ = 1 at the matching point; it is a tree statement at the spot values (Codex's own one-loop
  pole at 5.8e15 GeV is the relevant caveat). My brute-force minimisation shows the bound is loose by a factor 1.03 at
  s = s_c and 3.5 at s = 100 H_i, i.e. there is slack.
- **N7 (what a fixed-charge Hessian can and cannot show).** The complete radial statement is: L₊ has exactly one
  negative eigenvalue, L₋ ≥ 0 with F its nodeless zero mode, and dQ̄/dΩ < 0 (f(0) = −21.5 < 0). Those three together are
  the Vakhitov–Kolokolov condition; the track's text has all three but spreads them over §2. Nothing in either packet
  addresses ℓ ≥ 2 deformations or the gauged (e = 1e-6) spectrum; both say so.

## Limits of this skeptic pass

Radial, classical, tree level, zero temperature, scalar sector only; no 3D runs; the mass-calibration chain, charge
loss, isocurvature and the later packets are not touched. My centre values F₀, S₀ are only O(h²)-accurate (2e-6);
Codex's and the track's are better. Comparisons are retrospective; nothing here validates the model or shows the
fields exist.

## Reproduction (skeptic code only)

```
cd <scratch>; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
T=research/reality_framework/grut_formation_review_20261025/tracks/m1_vacuum/skeptic
python3 -I $T/code/sk_algebra.py   $T/results/sk_algebra.json                      # 7 s
python3 -I $T/code/sk_qball.py     $T/results/sk_qball.json primary_cc.npz         # 20 s
python3 -I $T/code/sk_hessian.py   $T/code $T/results/sk_hessian.json              # 2.7 min
python3 -I $T/code/sk_coldfluid.py $T/results/sk_coldfluid.json                    # 2 s
python3 -I $T/code/sk_wall.py      $T/results/sk_wall.json                         # 8 s
```
`results/sk_firstlaw_neighbours.json` holds the Q = 199000/201000 reproduction of Codex's 2.9e-7; `results/sk_checks.json`
holds the per-claim table above in machine form.
