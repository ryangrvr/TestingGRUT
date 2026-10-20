# Adversarial review: `closed_whole_demo.py` (S1-S8) and `results/closed_whole_demo.json`

Reviewed files (read-only):
- `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/code/closed_whole_demo.py` (sha256 `ba50f60e…`, equals `provenance.script_sha256` in the JSON)
- `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/results/closed_whole_demo.json`
- `/home/user/TestingGRUT/research/reality_framework/grut_closed_whole_theory_20261019/CALCULATIONS.md` (the notes, for item (g))

Check scripts and raw outputs live in `/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wf2/review1/`
(`lib.py` execs the repo's own utility functions; `a_ptrace_check.py` … `g_notes_check.py`; outputs `a_out.txt` … `g_out.txt`;
`repro/` holds a scratchpad re-run of the unmodified script). Run any of them with `python3 <script>` from that directory.

## Verdict

**Sound with minor issues.** No recorded number is wrong, every S1-S8 claim is supported by the numbers it cites, and the
whole JSON reproduces bit-for-bit. One recorded check (S8 subadditivity / Araki-Lieb) is vacuous as coded and should be
re-run on the mixed state; several sentences in the notes are stronger or looser than what the numbers show. Nothing
blocking; no new physics is claimed and none is implied by the numbers.

## 0. Reproduction and provenance

- `sha256sum closed_whole_demo.py` = `ba50f60e58a0219a228ea5829b2116c4e1ca3f2c7008e4ff2eb8f51de7090ed2` = JSON `provenance.script_sha256`.
- Re-run of an unmodified copy (`repro/code/closed_whole_demo.py`, 53.7 s, numpy 2.5.3 / scipy 1.18.1): **worst absolute difference
  in any recorded number = 0.0** (every scalar and all 401-point arrays, provenance excluded).

## (a) `ptrace` vs brute force — PASS

`a_ptrace_check.py`: explicit index-loop partial trace on a random 4-qubit density matrix and a random non-Hermitian operator,
all 15 keep-sets (contiguous and non-contiguous: `[0,2]`, `[0,3]`, `[1,3]`, `[0,2,3]`, …): worst deviation 9.9e-16. Unsorted
input `[3,0]` returns the documented ascending order (dev 7e-18). Product-state ordering test `kron(a,b,c,d)` keep `[0,2]` =
`kron(a,c)` (1e-16). 8-qubit cases used by the script (`[0]`, `[1..7]`, `[0,4,7]`) vs einsum: ≤ 1.7e-16.

## (b) S3 — inversion well conditioned; −0.159 is robust — PASS

`b_s3_check.py`:
- cond(Φ_t) at sample times recomputed = recorded (1.1137, 1.3533, 1.8787, 1.9077, 2.3495, 1.4123). Worst over the full grid
  t∈[0,40]: **2.858 at t=17.3**, min singular value 0.389. Inversion is never a concern.
- TCL generator vs finite-difference step: ‖K_h=1e-2 − K_h=1e-3‖ ≤ 7.6e-5, ‖K_h=1e-4 − K_h=1e-3‖ ≤ 7.6e-7 (O(h²) as expected).
  The recorded table recomputed with h=1e-4 agrees to 6 decimals.
- CP-divisibility (no finite difference involved): `inv` route vs `solve` route vs least squares from 10 inputs agree to 1e-15.
  Min Choi eigenvalue at t=5: τ=0.1 → −0.0351, τ=0.25 → −0.0858, **τ=0.5 → −0.1590 (recorded)**, τ=1.0 → −0.2614, τ=2.0 → −0.3896.
  Negative at every sampled t for τ ≤ 0.5; for τ ≥ 1 the first interval (t=0.5) is CP (+0.0018, +0.0166), which does not
  affect the conclusion (one non-CP interval suffices).
- τ-free cross-check: canonical (Kossakowski) rates of K(t) — one rate is negative at every sampled t (−0.0070 at t=0.5 …
  **−0.3523 at t=5** …, two negative at t=20); the τ=0.1 Choi eigenvalue ≈ τ·(most negative rate) (−0.0351 vs −0.0352).
  Sanity: rates of −i[Z,·] are 0, of (ZρZ−ρ)/2 are (0,0,1). This matches the inversion-free S2 BLP backflow (1.2035).

## (c) S4 — identical marginals, positive ρ₂, analytic inhomogeneity — PASS

`c_s4_check.py`, all three cases: Tr_A χ = Tr_Ā χ = 0 exactly; marginal deviations 0 / 0 / 2.2e-16; ρ₂ Hermitian, trace 1,
min eigenvalue 1.953e-3 / 1.953e-3 / 1.68e-9 (>0); D(ρ_A¹(0), ρ_A²(0)) = 1.6e-18 / 4.4e-16 / 6.9e-18 and later 0.1025 (t=0.7) /
0.0097 (t=14.1) / 0.0544 (t=0.7) — all equal to the recorded values.
Analytic inhomogeneity for χ = Y₀X₁: [Z₀, Y₀X₁] = −2i X₀X₁ (traces to 0 over Ā), [H_Ā, Y₀X₁] = Y₀⊗[H_Ā, X₁] (trace of a
commutator → 0), [gX₀X₁, Y₀X₁] = 2ig Z₀⊗1, so −i Tr_Ā[H, εχ] = 2gε·2⁷ Z₀ = 2·0.5·(1/512)·128 Z₀ = ¼ Z₀, Frobenius norm
√2/4 = **0.353553** = recorded 0.353553; numeric − analytic = 0.0; central difference of d/dt[ρ_A²−ρ_A¹](0) agrees to 3e-9.
Analytic ε_max = 1/256 (inf-T) and 1/(2(1+|⟨X₁⟩|)) = 0.5 (⟨X₁⟩ = 3e-17 because [H_Ā, ΣZ] = 0) match the bisection.

## (d) S6 — CNOT, identities, composition — PASS (with wording items below)

`d_s6_check.py`: ‖C†C − 1‖ = 0; C|a,mid,r⟩ = |a,mid,r⊕a⟩ on 20 random basis states; [H′, Z_R] = [H′, X_R] = 0.
Born p(Z_A=+1,t₁) = register pop₀ = 0.354909508264 (identity by construction: R starts in |0⟩).
**Why dephasing = record traced out:** Tr_R[CNOT(ρ⊗|0⟩⟨0|)CNOT†] = P₀ρP₀ + P₁ρP₁ is an exact operator identity (the CNOT with a
fresh |0⟩ ancilla is a Stinespring dilation of Z-dephasing on A) — verified 1.6e-15; and because H′ = H⊗1_R, U′(t) = U(t)⊗1_R
commutes with Tr_R — verified 8.7e-15. Hence 0.338531123755 = 0.338531123755 to all digits is a consistency check, not evidence.
Joint distribution sums to 1; its z_R marginal (0.35490951) equals the t₁ Born probability (perfect record).
Additional check not in the script: Φ_{t₂−t₁}∘D∘Φ_{t₁}(ρ_A(0)) (thermal bath re-prepared at t₁) gives **0.288331 ≠ 0.338531**
(diff 0.050), so the notes' claim "multi-time statistics are not functions of one-time maps" is true in this model but was
not recorded.

## (e) S7 — PASS (with wording item)

`e_s7_check.py`: max|Im H| = 0 and H = Hᵀ; U(t)* = U(−t) exactly (0.0); ρ_T* Hermitian (6e-17), min eig −2e-16, trace 1, same
spectrum as ρ₀ (3e-16); ρ_T* ≠ ρ_T (trace distance 0.588, so the demo is not vacuous); ρ_rev(t) = ρ(T−t)* to 5e-15, so
S_A(rev,t) = S_A(T−t) exactly and the recorded six values equal the S1 array read backwards. S_A along the reversed run is
**not monotone**: 201 of 400 grid steps increase, it peaks at 0.489 (t=2.5), 26 sign changes; it ends at exactly 0.

## (f) S8 — S_A = S_Ā and purity PASS; inequality check VACUOUS; closest-return statement fair

`f_s8_check.py`: for the pure global state S_tot ≤ 7e-15 and S_A = S_Ā (Schmidt), so subadditivity reads 0 ≤ 2S_A and
Araki-Lieb reads 0 ≥ 0: **neither can fail**; the recorded "0 violations" tests nothing beyond S_A = S_Ā. Non-trivial version
on the recorded mixed S1 trajectory (S_tot = 1.260931): min I(A:Ā) = 1.8e-15 (t=0, product state), min Araki-Lieb slack
S_tot − |S_A − S_Ā| = 0.000 at t=0 (saturated: S_A(0)=0), 0.2115 at t=3, 0.2439 at t=20; 0 violations.
Closest return: 0.401232 at t=11.1 reproduced; refined (step 0.001) 0.401223 at t=11.095; fidelity 0.839; long-time average
fidelity Σ|c_k|⁴ = 0.309 (typical D ≈ 0.83), 19 eigenstates with weight > 1e-3. The notes say "no recurrence time is
claimed" — fair; 0.40 is a partial return, not a recurrence.

## (g) Notes vs numbers

`g_notes_check.py`: every entry of the S1 table in CALCULATIONS.md matches the JSON to the printed digits; "first exceeds 0.4
at t = 8.8" true; Tr ρ² = 0.489575 true; E_A + E_Ā + E_int − E_tot ≤ 9e-11 (12-digit JSON rounding); S2 sign changes 26
reproduced; S5 b′ Bohr frequencies {0, 0.5, 1.5, 2.0, 2.5} follow from eigenvalues ±1.25, ±0.75 → period 4π.

## Defects (location → deciding check → minimal fix)

| # | Sev. | Location | Problem | Deciding check | Minimal fix |
|---|---|---|---|---|---|
| D1 | major | `closed_whole_demo.py` L556-565; JSON `S8_what_closure_forces.{subadditivity,araki_lieb}_violations`; CALCULATIONS.md S8 "subadditivity and Araki-Lieb never violated" | Both inequalities are checked only on a pure global state, where they hold identically (S_tot = 0, S_A = S_Ā); the recorded 0 carries no information. | `f_s8_check.py`: AL slack = S_tot − 0 ≤ 7e-15 at every sampled t; cannot be violated. | Evaluate both on the mixed S1 trajectory (arrays already recorded): min MI = 1.8e-15 ≥ 0, min AL slack 0 (t=0, saturated) then > 0.2; state that S8's pure run tests S_A = S_Ā and purity only. |
| D2 | minor | JSON `S7_time_reversal.note` ("its forward evolution lowers S_A back"); CALCULATIONS.md S7 "lowers the local entropy" | Endpoint statement reads as monotone; S_A rises to 0.489 (t=2.5) and increases on 201/400 steps before ending at 0. The reversed run is identically the forward run read backwards (ρ_rev(t) = ρ(T−t)*). | `e_s7_check.py` | Write "non-monotonically returns S_A to 0 at t = 40; ρ_rev(t) = ρ(T−t)* exactly". |
| D3 | minor | `closed_whole_demo.py` L503-514; JSON `S6.note_persistence`; CALCULATIONS.md S6 "makes it drift" | With κX_R the register is still uncoupled; the "drift" is free single-qubit precession (1-qubit prediction 0.250487 = recorded to 1e-15). Nothing about the record is lost — a rotated readout recovers it exactly. The claim is supported only for fixed-basis Z_R readout. | `d_s6_check.py` last block | Say "rotates the record out of the Z_R basis", or couple R to the bath (e.g. κX_R X_7) if genuine degradation is meant. |
| D4 | minor | JSON `S6.note_two_time`; CALCULATIONS.md "Reading" item 4 ("multi-time statistics are not functions of one-time maps") | Script compares only to the uninterrupted Φ_{t₂}(ρ_A(0)) = 0.4491; the composed prediction Φ_{t₂−t₁}∘D∘Φ_{t₁} was never computed. | `d_s6_check.py`: composed = 0.288331 vs exact 0.338531 (claim true, unrecorded). | Add the composed number to S6. |
| D5 | minor | JSON `S6.born_pZ_plus_at_t1` / `register_population_0_after_write` / `pZ_plus_A_at_t2_with_measure_and_forget_dephasing` and `note_two_time` | Register-population = Born probability and dephased = record-traced-out are exact identities (Stinespring dilation + idle register), not independent numerical evidence; the JSON note does not say so (CALCULATIONS.md does, partly). | `d_s6_check.py`: identities hold to 1.6e-15 / 8.7e-15. | Label them "identity (consistency check)". |
| D6 | minor | CALCULATIONS.md S4 "The first-order term is the Nakajima-Zwanzig inhomogeneity"; JSON `norm_of_inhomogeneous_term_at_t0` | The NZ inhomogeneity I(t) is time dependent; the recorded number is its t = 0 value. For χ = Z₀Z₁ both the first- and second-order terms vanish at t = 0 (−Tr_Ā[H,[H,εχ]] = 0 computed), the difference starts at O(t³) and stays ≤ 0.0097. | `c_s4_check.py` | Say "its t = 0 value"; note the classical case is O(t³). |
| D7 | minor | `DATUM["H_Abar"]` "non-integrable"; CALCULATIONS.md "The complement is nonintegrable … nothing below depends on a fine-tuned spectrum" | Asserted, not tested; 7 sites are too few (sector ⟨r⟩ = 0.41-0.56, inconclusive). Also unstated exact symmetries: [H_Ā, ΣZ] = 0 (why ⟨X₁⟩ = 0 and ε_max = 0.5 exactly) and [H, ∏Z] = 0 for the full H. | `g_notes_check.py` | Write "disordered XXZ (no known integrability)"; mention the U(1)/parity symmetries. |
| D8 | minor | JSON `S2.{S_A_max,t_of_S_A_max,S_A_min_after_max}`; CALCULATIONS.md S2 | "max 0.489 at t=37.5, falling to 0.404 by 39.5" covers only the last 2.5 time units; the claim is carried by the 26 sign changes and 13 local maxima (first fall 0.338 at t=2.1 → 0.235 at t=3.6). BLP sum of increments is a 0.1-grid quantity (lower bound). | `g_notes_check.py` | Quote the first local max/min; label the backflow sum "on the 0.1 grid". |
| D9 | minor | `closed_whole_demo.py` L138-141 (`embed_bath_state`, unused), L478 (`CNOT = kron_all(...)` overwritten) | Dead code. | inspection | Delete. |

## Confirmed (numbers support the claim as stated)

S1 invariants and Born local/global identity; reduced equation check (H_Ā term traces to 3e-16); S2 non-monotone S_A and BLP
backflow; S3 CPTP for product preparations, time- and bath-state-dependent TCL generator, CP-divisibility failure (robust to
τ and inversion route, corroborated by negative canonical rates and by BLP); S4 identical marginals / different futures with
the hand-derived inhomogeneity; S5 non-unique local law and exact 4π recurrence of the commensurate complement; S6 record as an
included subsystem, post-intervention statistics differ from the uninterrupted ones; S7 exact reversibility and legal
conjugated preparation; S8 S_A = S_Ā and constant purity for the pure run; closest-return statement is fair.
