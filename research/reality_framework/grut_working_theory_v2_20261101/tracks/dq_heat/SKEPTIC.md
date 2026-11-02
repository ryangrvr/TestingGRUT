# SKEPTIC of track dq_heat: the hidden heat budget R(e) of the M1r release

2 November 2026 · Claude, adversarial skeptic in the decisive pass (BRIEF_DECISIVE.md; issue #13 context). This is
same-session AI work on one supplied model. It is not outside replication, and nothing here is evidence for the hidden
sector, a validation or an exclusion. Every input (e, the inflaton history, the portal condition, the static slice) is a
model choice. Labels: [D] derived, [C] computed with our own code, [E] estimate, [A] argued, [Q] quoted, [R] recalled.

Code is in `skeptic/code/`, results in `skeptic/results/` (total 0.3 MB). All of it is my own code, written fresh,
except `dq_model_copy.py`, an unchanged copy of `dq_heat/code/dq_model.py` used only for its stationary-point solver.

## Verdict

- **The heat budget at the static slice stands.** I recomputed every load-bearing number by a different method or code
  and found nothing wrong with R_total(e), N_A, K, the mediator rates, S/Q(R) or e_crit.
- **One claim is refuted.** The estimate that a Hubble rate growing with φ_I would "nearly cancel" the exit-state
  sensitivity (factor 1.06–1.15) is wrong. Re-solving the prepared state at the exit Hubble rate gives f = 0.77–0.40.
- **Consequence (new, [C] on the paths, [E] for the paths themselves).** On every rolling-inflaton path I tried, the
  declared e = 10⁻⁶ lands in the **Q-matter** branch, with S/Q = 3.11–3.37. The free-X verdict at 10⁻⁶ therefore holds
  only if φ_I is frozen over the last ~4 e-folds before the slice.

## 1. Recomputations

| claim (dq_heat) | their value | my value | my method | verdict |
|---|---|---|---|---|
| N_A(10⁻⁶), single-Higgs-direction Proca | 2.0249 | **2.0265** (+0.08%) | `sk_modes.py`: Stueckelberg angle θ in **cosmic** time; exact massless de Sitter Goldstone start at k/aH = 50 (theirs: adiabatic start at 1000, conformal time); splines in t, not η; two read-outs; Wronskian ≤ 10⁻⁸ | confirmed |
| N_A(10⁻⁵) | 0.6300 | **0.6304** | same | confirmed |
| per-mode n_k, k ≤ 5.3 | – | 0.05–0.9% agreement on all 61 modes | same | confirmed |
| frozen-angle factor (F_vac/F_inf)², K^{3/2}, e^{−1/2} | as stated | n_k(theirs)/n_k(frozen) = **0.9967–1.0001** for k = 2×10⁻³–0.05, in all 11 cases (full, noX, inst; 10⁻⁷–10⁻⁵) | `sk_frozen.py`: semi-analytic. The frozen Goldstone angle (a global rotation of the background solution) feeds a universal radiation-era equation, d/dx[x²/(q²+x²)θ_x] + x²θ = 0, giving N_A = (F_v/F_i)² K^{3/2} H² m^{−1/2} I/(8π²) with I = ∫qC²dq = 1.1815 [D/C] | confirmed. Modes k = 0.05–1 fall 1–6.5% below frozen (near-horizon modes); this is the 2.7% "full/inst ≠ K^{3/2}" |
| K = a²H | 1.98581 | **1.985805**, all digits | `sk_K.py`: inflaton-only background (φ released at rest, (3, 0.5)), independent of rptri_dq.c | confirmed |
| "0.70 of the GMR formula" | 0.70 | frozen solution with F_v = F_i is 0.938 of the GMR fit formula; 0.70 = 0.7565 (frozen angle) × 0.93 | `sk_frozen.py` | consistent. GMR quote a one-digit coefficient (6 μeV) [R, arXiv:1504.02102] |
| rates | Γ_L, Γ_H = 0.0350, Γ_σ = 0.398, Γ_A = e²m_A/48π, ⟨a⟩ = (πK/2Γ)^{1/2} | rederived from the Lagrangian: Y_L real ⇒ X³ and X̄³; √2μ G Y_L Y_H vertex; λ_Σ hσG²; Rayleigh mean in radiation | [D] | confirmed; R_Y,hom = 0.37887, R_cascade = 0.01855 |
| S/Q column | 4.95, 4.01, 3.51, 2.99, 2.55, 2.20 | 4.951, 4.010, 3.512, 2.988, 2.548, 2.200 | dq_entropy's independent closed form S/Q = 0.190255(W − 0.97494)^0.75013 − 0.00428 [Q] | confirmed |
| e_crit (bulk, Q = 2×10⁵) | 2.06, 3.20 ×10⁻⁶ | 2.05, 3.17 ×10⁻⁶ | direct solve with R_vec ∝ e^{−1.007}, R_c ∝ e^{−1/2} | confirmed (log-log interpolation error ≤ 1%) |
| M_P portals | S/Q 4.09; N_A × 0.945 | 4.094; (275.93/283.85)² = 0.9450 against the computed 0.9447 | closed form | confirmed |
| timing | a_on = 5.25×10¹², T_vis 81 GeV, ×895 | 5.253×10¹², 81.4 GeV, ×895 | n_Q = 151 m_X³ frame relation | confirmed |
| gate1 docstring wrong | – | gate1_vector.py line 13 indeed says "N_A ∝ e^{1/2} … R_A is e-independent". GMR give n/s ∝ m^{−1/2} | read | confirmed |
| fixed-H exit family | S/Q 4.24, 3.47, 2.72 at s = 1.25, 1.5, 2 | 4.235, 3.467, 2.721 | `sk_exit.py` | confirmed |

**Piece B against the multi-field modes [C].** dq_heat's own `ymodes` runs (rpmodes.c, e → 0, exact de Sitter BD start)
contain the Σ-phase eigenmode, which is A_L for k/a ≫ m_A. I compared it with the Proca band (`sk_pieceB_omitted.py`):
- On k ∈ [1, 30], XYS (τ = 40) gives **0.0498**, against the Proca "full" band's 0.0512 (ratio 0.973). YS gives
  0.019–0.022, against "noX" 0.0236. The band total stands.
- **The k = 4.45 "X-orbit resonance" is absent** from the multi-field Goldstone: n = 2.9×10⁻³ at k = 4.71, against
  1.8×10⁻² in the Proca band. The multi-field number sits at k = 1–3 instead. The resonance is a single-direction artefact;
  the integrated number agrees by coincidence of size. Numerical consequence: |ΔR| ≤ 0.003.
- The XYS Goldstone on k ∈ [15, 80] is 0.0104, which equals N_YH = 0.0104. The cascade products are resolved and counted
  once, in R_cascade.
- The M_P-portal background (c = −0.337) gives the same release band (0.0505). So inflaton → A_LA_L through the induced
  portal is not a hidden source.

## 2. Refuted: the "nearly cancelling" estimate for H growing with φ_I

dq_heat §5.1 ([E]; in the JSON, `inflaton_exit_relevance`) multiplies the fixed-H factor (F(1)/F(s))² by s². That keeps
F(s) from the fixed-H family. But the prepared X condensate is held by ξ_X R|X|² = 12ξ_X H²|X|², so F grows with H too.

**Recomputation (`sk_exit.py`).** I re-solved the conditional minimum at u = μ + s(10i − μ) and R = 12h_r² (H_exit = h_r H_i;
physical M_Y, h, μ, κ, λ unchanged). The piece-A factor is f = [h_r F(1,1)/F(s, h_r)]². It replaces the frozen angle
H/(2πF) at exit. All exit states remain heavy (lightest m²/H_exit² ≥ 10).

| s = h_r (H ∝ φ_I) | 1.25 | 1.5 | 2 | 2.5 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|---|
| dq_heat estimate s²(F₁/F_s)² | 1.150 | 1.055 | 0.911 | 0.833 | – | – | – |
| **re-solved f** | **0.774** | **0.623** | **0.507** | **0.461** | 0.436 | 0.410 | 0.397 |
| S/Q(10⁻⁶) | 4.34 | 3.92 | 3.58 | 3.44 | 3.37 | 3.29 | 3.25 |

- H growth alone (s = 1, h_r = 1.25–2.5) gives f = 1.04–1.07: F/H_exit is nearly invariant because |X| ∝ H.
- What sets the vector heat is the holder value at exit, s, not H.
- Along H ∝ φ_I, free X at 10⁻⁶ holds only for **s < 2.37**.

## 3. New issues

**N1 [C on the paths; E for the paths; A for the bounds]. The declared branch is the static-slice extreme.**
- Free X at e = 10⁻⁶ needs f_eff ≥ 0.470, i.e. (H/F)_exit ≥ 0.686 of its slice value.
- A rolling inflaton has |φ_I| decreasing (s ≥ 1), and H is non-increasing (h_r ≥ 1). On that region f ≤ 1.07, so S/Q(10⁻⁶)
  ≤ 5.14 for any history.
- Weighting dq_heat's per-mode spectrum by f(N(k)), with N(k) the exit e-fold (`sk_exit_path.py`):

  | path [E] | s, h_r at N = 3.8 | f_eff | S/Q(10⁻⁶) | e_crit (bulk) |
  |---|---|---|---|---|
  | static slice (M1r as specified) | 1, 1 | 1 | 4.95, free X | 2.06×10⁻⁶ |
  | fixed H, s = 1 + 0.2N | 1.76, 1 | 0.391 | 3.23, Q-matter | 0.86×10⁻⁶ |
  | quadratic large field, slice at ε = 1 | 2.93, 2.93 | 0.436 | 3.37, Q-matter | 0.95×10⁻⁶ |
  | plateau (Starobinsky-type, slice where m_I = 3H) | 2.02, 1.34 | 0.351 | 3.11, Q-matter | 0.81×10⁻⁶ |

- N_A is dominated by modes near k_*, so f_eff ≈ f(N ≈ 3.5–4).
- On these paths W_vec is 0.35–0.42 of wt_critic's broad transfer. WORKING_THEORY §8's failure condition (< 0.5) **is met**.
- s = 1 requires φ_I frozen for the last ~4 e-folds, i.e. a hybrid-type history in which φ_I is not the rolling inflaton.
  A field with m_I = 3H_i cannot sit still for e-folds by itself [A].
- The paths keep the M1r release (N_X, E_X, K) fixed. A rolling inflaton at the slice would also change these. Less
  post-slice acceleration means a smaller K, and R_vec ∝ K², which points the same way [A, not computed].

**N2 [C]. K sensitivity.** R_vec carries K² at fixed E_X. K = 1.75–2.97 over (m_I, Γ_I) = (2–4, 0.25–1) (`sk_K.py`), i.e. a
factor 0.78–2.24 in R_vec before E_X changes. Only (3, 0.5) has a computed release.

**N3 [C]. Piece-B attribution.** See §1. The band number is confirmed; the resonance explanation is not supported.

**N4 [E]. Channels missing from the ledger.** Each vector carries 1.78 E_X at 10⁻⁶, so I estimated every vector-producing
channel I could find (`sk_pieceB_omitted.py`). All are negligible:
- inverse-decay freeze-in XX̄ → A while T_h ≫ m_A (a ≲ 10⁵): N_A ≈ 10⁻⁴, R ≈ 1.8×10⁻⁴ at 10⁻⁶ and 5.6×10⁻⁴ at 10⁻⁵
  (∝ e^{1/2});
- Y_L → 3X + A_L through a virtual Y_H, via the unsuppressed √2μ vertex: BR ≤ 8×10⁻⁵, R ≤ 8×10⁻⁵;
- the condensate-induced two-body decay Y_L → XX (vertex 3κ⟨X⟩): 0.45% of Y_L decay early, so R_Y falls by ≤ 1.7×10⁻³;
- gravitational freeze-in from the SM bath: R ≲ 10⁻⁵;
- the Coulomb energy of the fragmenting charge: ~10⁻⁶.

**N5 (minor).** The M_P variant uses the coarse k > 1 grid (no dense band); the difference is 2.5% of the band, i.e. ΔR ≈ 0.002.

## 4. What I did not check

- R_Y,k (finite-k Y_L, ≤ 4×10⁻⁴) and the separate-universe piece B (3×10⁻⁴) were not recomputed. Both are negligible.
- The longitudinal mode is still the single-direction action for k ≤ 1. The frozen-angle symmetry argument makes this
  exact for superhorizon modes, and the X part of F² is 10⁻⁴ at re-entry.

## Reproduce (one or two cores)

```bash
cd research/reality_framework/grut_working_theory_v2_20261101/tracks/dq_heat/skeptic/code
export OMP_NUM_THREADS=1
python3 -B sk_frozen.py                 # 4 min
python3 -B sk_modes.py 1e-6 full dq 0 dq && python3 -B sk_modes.py 1e-5 full dq 0 dq   # 3 min each, 2 workers
python3 -B sk_K.py && python3 -B sk_exit.py && python3 -B sk_exit_path.py && python3 -B sk_pieceB_omitted.py
```

`sk_modes.py` reads dq_heat's scratch background npz (made by `dq_heat/code/background.py`).
