# Track rp_preparation: handoff item 2 of Codex e160356 (the preparation obstruction)

26 October 2026 · Claude (independent review track, issue #13) · packet e160356
(`origin/codex/grut-curvature-boundary-20261026`, `research/reality_framework/grut_reheating_preparation_20261026/`),
read only through `git archive e160356 … | tar -x` into the scratchpad. Nothing was committed. Own code in `code/`,
outputs in `results/`. No code of Codex's and no code of our earlier tracks is reused (the finite quantum model has no
counterpart in earlier tracks). Codex's `boundary.py` was run once, in the scratchpad copy, as a separately labelled
reproduction (§3). Files read: CLAUDE_HANDOFF.md, REPORT.md, MATHEMATICS.md, BOUNDARY.md, PORTALS.md, PROTOCOL.md,
SOURCES.md/json, REPRODUCE.md, BOUNDARY_RESULTS.json, VERIFICATION_RESULTS.json, and the heads of REHEATING_RESULTS,
REHEATING_VERIFICATION, INDEPENDENT_REHEATING and PORTAL_RESULTS; items 1, 3 and 4 of the handoff are outside this track.

**One finite toy model with supplied frequencies, couplings and coherent amplitudes. Nothing here is evidence that X,
Y, Σ or the boundary state exist. This is not a validation, an exclusion or a GRUT prediction.** All agents are
same-session AI, so this is not outside replication. Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate,
[A] argued, [R] recalled (not re-read this session).

## 1. What Codex claims [Q]

Three bosonic modes a, b, c with charges +1, +3, −6; cutoffs (18,8,6) and (24,12,9); the exactly neutral projection
P₀ (n_a + 3n_b − 6n_c = 0) of the coherent state |√6, √2, √2 e^{iδ}⟩; neutral dimensions 31 and 58;
H = n_a + 1.7 n_b + 2.3 n_c + 0.025 (b†a³ + h.c.) + 0.04 (cb² + h.c.); light charge L = n_a + 3n_b; torque
d⟨L⟩/dt = i⟨[H, L]⟩. Claims: (i) opposite phases δ = ±π/2 have identical energy, occupations, single-mode entropies
and zero hidden/SM mutual information (SM factor diag(.2,.8)) but opposite torque ∓1.3563 (cutoff 1), ∓1.3576
(cutoff 2); (ii) a doubled conjugate copy has zero initial signed light charge and signed torque ±2.7153, cutoff change
0.102%, energy change 0.00374%; (iii) a five-reading open clock with one time-independent constraint C satisfies CV = 0
and the physical effects E_{a|t} = 5P(|t⟩⟨t|⊗Π_a)P commute with C and reproduce the conditional Born probability
(n_a ≥ 6 at reading 2); (iv) therefore "neutrality, these entropy conditions and an initially uncorrelated SM factor
do not select the charge-producing phase"; (v) homogeneous FLRW has zero Weyl curvature for every phase, so a Weyl
penalty selects nothing; (vi) a real diagonal Lorentzian boundary phase preserves configuration probabilities; (vii)
scope: not a continuum constraint solution, not a derived cosmic phase, the projected state need not be Gaussian.

## 2. Method (own code) [C]

`code/neutral_boundary.py`: ladder operators on the full truncated Fock space as Kronecker products (dimensions 1197,
3250, 6448, 11248 for cutoffs (18,8,6), (24,12,9), (30,15,12), (36,18,15)), H assembled from them, then restricted to
the neutral block; the off-block part of H is checked to vanish (charge conservation). Codex instead enumerates the
neutral basis and writes matrix elements by hand. Computed per phase: projection probability (relative to the
truncated and to the untruncated coherent state), occupations, gauge charge, energy, light charge, torque, Shannon
entropies of the occupation marginals and von Neumann entropies of the single-mode reduced states (partial trace in
the full space), hidden/SM mutual information. Added checks: complex conjugation in the number basis maps ψ_δ to
ψ_{−δ}; the rotation identity e^{iθL}ψ_δ = ψ_{δ+6θ}; a 49-point δ scan of E(δ), T(δ) and 6 dE/dδ (finite differences);
torques of all H eigenstates and of Gibbs states; the explicit doubled product state (961 and 3364 dimensions) and its
swap-symmetrized version; the clock (isometry, CV, spectrum of C, invariance under clock reversal times conjugation,
effect commutators and Born matching for three events at all five readings); a complex cb² coupling 0.04 e^{iχ}
(χ = 0.3, 1.0); an SM-correlated preparation. `code/projection_identity.py`: the exact identity
⟨cb²⟩_projected = α_c α_b² at machine precision.

## 3. Reproduction

| Quantity (cutoff (18,8,6) / (24,12,9)) | Codex BOUNDARY_RESULTS | Own code | Codex code rerun (reproduction) |
|---|---|---|---|
| Neutral dimension | 31 / 58 | 31 / 58 (88, 124 at (30,15,12), (36,18,15)) | identical |
| Neutral projection probability | 0.0410786 / 0.0408848 | 0.0410786 / 0.0408848 (relative to the truncated coherent state; untruncated 0.0408829) | identical |
| Energy δ = 0, ±π/2, π | 14.310439, 14.084395, 13.858350 / 14.311195, 14.084921, 13.858647 | same to all digits | max abs difference 3.4e-15 over the whole JSON |
| Mean occupations (±π/2) | 5.833525, 1.784986, 1.864747 | same | identical |
| Mode entropies (nats) | 1.154788, 1.595187, 1.049856 | same; von Neumann of the reduced single-mode states equal to 4e-16 | identical |
| Torque at δ = +π/2 | −1.3562654 / −1.3576447 | −1.3562654 / −1.3576447; converged −1.35764502 at (30,15,12) and (36,18,15) | identical |
| Doubled signed torque | −2.7125308 / −2.7152894 | −2.7125308 / −2.7152894 (explicit product state); converged −2.7152900 | identical |
| Hidden/SM MI | ≤ 4.5e-16 | ≤ 7e-16 | identical |
| Clock: isometry, CV, effect commutator, Born mismatch | 2.5e-14, 6.2e-15, 4.9e-15, 2.5e-15 | 2.3e-14, 5.7e-15, 5.1e-15, 4.3e-15 (3 events × 5 readings) | identical |
| Event n_a ≥ 6 at reading 2, δ = +π/2 | 0.7184909 | 0.7184909 | identical |
| Cutoff change (energy, light charge, torque) | 3.74e-5, 3.98e-5, 1.016e-3 | 3.74e-5, 3.98e-5, 1.016e-3 relative to the largest cutoff; (24,12,9) is within 2.4e-7 of convergence, (30,15,12) within 1e-10 | identical |

Every number in BOUNDARY_RESULTS.json is reproduced by independent code. Codex's "refined torque ±2.71528939" is the
(24,12,9) value; the converged value is ±2.7152900 (difference 2.4e-7 relative, inside their stated cutoff change).

## 4. The logic: is "these conditions do not select the phase" a theorem?

**Yes, for the listed conditions, and it does not depend on the finite example [D].** Let K be complex conjugation in
the occupation-number basis (antiunitary, K² = 1). Facts, all checked numerically: (a) H is a real symmetric matrix in
that basis (imaginary part 0, asymmetry 0), so KHK = H; (b) Kψ_δ = ψ_{−δ} up to phase (overlap 1 − 2e-16), because
the coherent amplitudes are real except α_c = √2 e^{iδ}, and P₀ is real; (c) L, the n_k and Q are real diagonal, so
K-even; (d) the torque operator i[H, L] is K-odd: K(i[H,L])K = −i[H,L]. Hence for every state ψ: ⟨f(H)⟩, all
occupation distributions, the spectra of all reduced states (Kρ_AK has the spectrum of ρ_A), hence all entropies and
mutual informations, and the total charge are equal for ψ and Kψ, while T(Kψ) = −T(ψ). If every selection condition
is K-even and H is K-invariant, the set of admitted states is closed under K and its torques come in ± pairs; a rule
that returns a unique state returns a K-invariant one, whose torque is zero. Codex's four rows are one instance
(energy difference 0, torque sum 0, entropy difference 0 at both cutoffs; the δ scan gives E(δ) − E(−δ) ≤ 7e-15 and
T(δ) + T(−δ) ≤ 1.2e-15 on 49 points).

**Label correction.** Codex calls the flip "CP"/"charge-conjugate copy". In the three-mode model no mode has an
antiparticle partner, so charge conjugation is not an operator on this Hilbert space; the symmetry actually used is
antiunitary conjugation in the number basis, i.e. the time-reversal-type transformation (torque = a velocity, odd;
energies, occupations, entropies = even). In the field theory both C (or CP) and T flip dQ/dt, with different
even/odd assignments (C flips Q and swaps particle/antiparticle occupations; T keeps Q and occupations). The theorem
holds for either, under invariance of the dynamics. gp_phase already states the classical C version (z ↔ z* with
real couplings) and its quantum corollary ⟨Q⟩ = 0 for C-invariant states [Q gp_phase §2]. CORRECTED (label only).

**Two sharper statements the example contains but does not state [D, C].**
1. *The torque is the phase gradient of the energy.* L is diagonal, so [L, P₀] = 0, and e^{iθL} rotates the
   gauge-invariant phase: e^{iθL}ψ_δ = ψ_{δ+6θ} exactly (overlap 1 − 2e-16). Therefore
   T(δ) = i⟨[H,L]⟩ = d/dθ ⟨ψ_δ|e^{−iθL}He^{iθL}|ψ_δ⟩ = 6 dE/dδ (finite-difference check: |T − 6E′| ≤ 3e-9 at all
   cutoffs). Consequence: **any energy-extremal member of the phase family has zero torque.** Energy does not only fail
   to select the sign; used as a minimization rule it selects the torque-free phase δ = π (E minimum, T = 1.7e-16),
   the finite quantum analogue of hw_gates gate 2 ("the phase sits at the torque-free minimum δ = 0, where ε = 0")
   and of gp_phase's "own-potential alignment gives zero charge".
2. *Exact projection identity.* With P₀ = (2π)⁻¹∫dθ e^{iθQ} and cb² a neutral annihilation monomial,
   cb²|α_θ⟩ = α_cα_b²|α_θ⟩ for every θ, so ⟨cb²⟩_projected = α_cα_b² = 2√2 e^{iδ} exactly (error 7e-14 at (36,18,15);
   6.8e-7 at (24,12,9) is truncation). Hence E(δ) = E₀ + 2·0.04·2√2 cos δ = E₀ + 0.226274 cos δ and
   T(δ) = −12·0.04·2√2 sin δ = −1.357645 sin δ, independent of the projection; the mean occupations are not
   preserved (6, 2, 2 → 5.8335, 1.7850, 1.8647) because n_k is a mixed monomial. The ±π/2 comparison is the
   δ → −δ symmetry of these two exact formulas.

**Stationary selections give zero torque regardless of symmetry [D, C].** Tr(ρ[H,L]) = Tr([ρ,H]L) = 0 whenever
[ρ, H] = 0. Every eigenstate of H (max torque 0.0) and Gibbs states at β = 0.1, 1, 10 (torque 0.0) confirm it. Ground
state, thermal state, the maximum-entropy state on the physical clock subspace (RECORDS.md's P/8) and any unitarily
equivariant rule F(H) (calc 1 S9) are of this type. This is the finite version of the third Sakharov condition.

**What breaks the theorem (each tested or argued):**
- *A K-odd (CP-violating) term in H.* With 0.04 → 0.04 e^{iχ}: E(+π/2) ≠ E(−π/2) (14.0176 vs 14.1512 at χ = 0.3), so
  energy now distinguishes the two phases, as the theorem's hypothesis says. But by statement 1 the energy minimum
  moves to δ = π − χ (2.845 vs 2.842 on the grid) and the torque there is 0.004, zero up to the grid spacing. CP
  violation in H alone does not make an energy-selected state transfer charge; the state must be displaced from the
  extremum, dynamically (a moving minimum: gp_phase found no resolved net charge, |⟨ε⟩| ≤ 0.034) or by hand.
- *An asymmetric initial state.* This is what Codex does: δ is supplied. The doubled state ψ_+⊗ψ_− breaks the copy-swap
  symmetry; its swap-symmetrized version has torque 0.0 (own check). The "zero initial signed light charge" of the
  doubled example is just ⟨L⟩_+ = ⟨L⟩_− (K-evenness of L) and the signed torque is twice the single-copy value
  (−2.7125308 = 2 × (−1.3562654)); it adds no information beyond the single-copy rows.
- *A K-odd observer condition*, e.g. a condition on the sign of the torque itself (circular), or a correlated SM factor
  carrying a CP-odd SM asymmetry. Own illustration: hidden sign recorded in an SM bit with weights (p, 1−p): mutual
  information 0.69 nats at p = 0.5 with hidden torque 0; torque −0.27 at p = 0.6 and −1.085 at p = 0.9, i.e.
  T_net = (2p−1) T_+. Correlation is not what selects; the SM-side asymmetry is. Codex's MI = 0 is imposed by the
  product form ρ_h ⊗ diag(.2,.8), not a tested condition: CORRECTED (status of the MI check).

**The clock embedding selects nothing [D, C].** V is a fixed linear isometry; every ψ gives a physical history with
CV = 0 (constraint residual ≤ 1.2e-15). The spectrum of C is 1 − cos(πk/5), k = 0…4, each d-fold: kernel dimension
31/58 (= d), first gap 0.190983 = 1 − cos(π/5) (Codex's 8-cycle in RECORDS.md had 0.292893). C is invariant under
clock reversal combined with complex conjugation (5.6e-16), so the K-argument extends to the history formulation,
with the reversed history starting from K of the final state. The effects are positive, sum to P, commute with C
(≤ 8.6e-15) and give the direct conditional probabilities (≤ 4.3e-15 for three events at five readings), as claimed.
Codex's statements that this is "not a periodic clock" and "not a cosmic time endpoint" are accurate.

## 5. Relation to GRUT: what is new, what is restated

- FORMAL_CORE DEF 1.1 admits every ρ_Ω(0); DEF 1.4 "restricting does not select the state"; §4 "a preferred initial
  state: BY DEFINITION"; calc 1 S9: a unitarily equivariant F(H) is stationary, so one-time marginals do not evolve;
  §7 state-selection row: a selection must break equivariance; COND 5.3 excludes embedding-covariant rules in locally
  covariant QFT. THM 2.1 (global invariants) is the reason energy, spectrum and entropy of the whole are conserved and
  therefore cannot encode a time-asymmetric sign in a stationary whole. Codex's example is an instance of these rows
  with a *discrete* symmetry: K-equivariance of the rule already forces zero torque, a weaker hypothesis than S9's full
  unitary equivariance but the same conclusion ("stationary ⇒ no one-time evolution" becomes "K-symmetric ⇒ no
  charge-transfer sign"). Substance: a restatement in a concrete gauge-neutral quantum state.
- RECORDS.md (parent 04c53f5) already showed "closure leaves the state free" in the 8-cycle clock (any record
  probability p; P/8 gives no changing records; a supplied filter F₀ selects). The new packet moves the same
  construction to an open 5-reading chain with charged modes. Nothing about selection changes.
- The closure fork (CLOSURE_SELECTION_RESPONSE.md: P1 adds no state-selection condition) and the vacuum pass
  (SYNTHESIS §5 row 11: "neither is selected by closure"; no state-selection rule found) are consistent with it.
- The field-theory tracks already contain the physical content: gp_phase (real couplings ⇒ C symmetry ⇒ ⟨Q⟩ = 0 for
  C-invariant states; zero physical phases at dimension ≤ 4 with Z₂_φ; the real subspace is unstable, giving
  random-sign local charge with zero mean), hw_gates gate 2 (P(π/2)/P(0) = e^{−43,277} in the de Sitter equilibrium;
  a holder with a supplied CP phase is required), sp_phaseprep (symmetric equilibrium distribution; nil probability
  for the observable universe). The finite example is the 0-dimensional quantum version of these statements.
- **New here (Codex):** an explicit exactly neutral, pure, normalizable boundary state (the "positive normalizable
  state" the §7 row demands) and its open-chain Page–Wootters embedding with checked physical effects. **New here
  (this track):** the identities T = 6 dE/dδ and ⟨cb²⟩_proj = α_cα_b², which say precisely why energy-type selection
  fails: it selects torque-free phases, with or without CP violation in H; and the clock-constraint spectrum and
  reversal invariance. **Not new:** the conclusion that the listed conditions do not choose the cosmic phase. Codex's
  own scope sentences (no continuum constraint solution, no derived phase, no new closure law) are accurate.

## 6. Verdicts

| # | Codex e160356 statement | Ours | Verdict |
|---|---|---|---|
| 1 | Neutral dimensions 31 and 58; exactly neutral, normalizable, positive; every transition conserves charge | 31, 58 (88, 124 at larger cutoffs); off-block H = 0; ⟨Q⟩ = 0 [C] | CONFIRMED |
| 2 | Projection probabilities 0.04108, 0.04088 | same; they are relative to the truncated coherent state; untruncated 0.0408829 [C] | CONFIRMED (normalization stated) |
| 3 | Equal energy, occupations, mode entropies for δ = ±π/2 (differences 0) | reproduced to all digits; E(δ) = E₀ + 0.226274 cos δ exactly, so energy is even in δ but selects δ = π if minimized [C, D] | CONFIRMED; sharpened |
| 4 | Mode entropies 1.1548, 1.5952, 1.0499 nats | same; equal to the von Neumann entropies of the single-mode reduced states (superselection), 4e-16 [C] | CONFIRMED |
| 5 | Torque ∓1.3563 / ∓1.3576; sum 0; cutoff change 0.102% | same; converged ∓1.35764502 = 12·0.04·2√2 (exact identity); (18,8,6) is 0.1016% low [C, D] | CONFIRMED |
| 6 | Doubled example: initial signed light charge 0, signed torque ∓2.7125 / ∓2.7153, reversed by phase exchange | explicit product state: 5.6e-16, ∓2.7125308 / ∓2.7152894 = 2 × single copy; converged 2.7152900; swap-symmetrized state torque 0 [C] | CONFIRMED (adds nothing beyond the single-copy rows) |
| 7 | Hidden/SM mutual information zero for both phases; "this condition also fails to choose a phase" | ≤ 7e-16, but imposed by the product form, not tested; a correlated SM bit gives T_net = (2p−1)T_+ [C] | CORRECTED (status: assumption, not a test) |
| 8 | Clock: V isometry, CV = 0, effects positive, sum to P, commute with C, Born match (0.7184909) | 2.3e-14, 5.7e-15, 5.1e-15, 4.3e-15 for 3 events × 5 readings; 0.7184909; spectrum 1 − cos(πk/5), kernel = d, gap 0.190983; C invariant under reversal × conjugation [C, D] | CONFIRMED; extended |
| 9 | "Neutrality, these entropy conditions and an uncorrelated SM factor do not select the charge-producing phase" | a theorem for K-even conditions under K-invariant H (not example-dependent); the transformation is conjugation in the number basis (time-reversal type), not charge conjugation; energy extremization selects a torque-free phase (T = 6E′); stationary rules give T = 0 unconditionally [D, C] | CONFIRMED as theorem; label CORRECTED |
| 10 | CP-violating action or asymmetric state needed (implicit) | χ ≠ 0 lifts E(+π/2) = E(−π/2) but the energy minimum stays torque-free (T = 0.004, grid-limited); displacement or a K-odd rule is what produces a sign [C] | CONFIRMED; sharpened |
| 11 | FLRW has zero Weyl curvature for all phases; a Weyl penalty selects nothing in the homogeneous class | Weyl ≡ 0 for every FLRW metric [D]; credit Penrose [R] | CONFIRMED (trivial) |
| 12 | A real diagonal Lorentzian boundary phase preserves configuration probabilities | |e^{iφ(x)}ψ|² = |ψ|² [D]; no boundary-term model exists to test more | CONFIRMED (trivial); rest NOT CHECKED |
| 13 | Scope: not a continuum constraint solution, not a derived cosmic phase, projected state need not be Gaussian | agreed; a Q-eigenstate with nontrivial number statistics is not Gaussian [A] | CONFIRMED |
| 14 | Credits: Page–Wootters named in RECORDS.md text; SOURCES.json pins only Ashtekar–Gupt for this item | no pin for Page–Wootters, conditional-probability clocks, coherent-state/charge-projection or Gaussian-state literature | CORRECTED (credit gap) |

No numerical error was found in the boundary part of e160356.

## 7. Credits and read levels

Page & Wootters, Phys. Rev. D 27, 2885 (1983); Wootters, Int. J. Theor. Phys. 23, 701 (1984) [R]. Conditional
probabilities in constrained systems: Gambini, Porto & Pullin (2004); Dolby (2004); Giovannetti, Lloyd & Maccone,
Phys. Rev. D 92, 045033 (2015); Smith & Ahmadi, Quantum 3, 160 (2019); Höhn, Smith & Lock, Phys. Rev. D 104, 066001
(2021) [R]. Coherent states and projections onto charge sectors: Glauber (1963); Perelomov, *Generalized Coherent
States* (1986); Klauder & Skagerstam (1985); Sanders, Bartlett, Rudolph & Knight, Phys. Rev. A 68, 042329 (2003)
(number superselection and projected coherent states); Bartlett, Rudolph & Spekkens, Rev. Mod. Phys. 79, 555 (2007)
[R]. Gaussian states: Weedbrook et al., Rev. Mod. Phys. 84, 621 (2012) [R]. Antiunitary symmetries: Wigner (1932,
1959) [R]. Charge generation conditions: Sakharov, JETP Lett. 5, 24 (1967); Dine & Kusenko, Rev. Mod. Phys. 76, 1
(2004) [R]. Weyl curvature hypothesis: Penrose (1979) [R]. None of these was re-read in this session; Codex's packet
pins none of them for the boundary item.

## 8. Limits

- A three-mode toy with supplied frequencies, couplings and amplitudes; no gravitational constraint, no gauge field,
  no continuum limit, no antiparticle modes (so no C operator). The K-argument needs only a real H in some basis and
  K-even conditions; whether the cosmological preparation conditions GRUT might adopt are K-even (or C-even) is a
  question about those conditions, not settled here.
- "Energy extremization selects a torque-free phase" is exact for the one-parameter family ψ_δ; a selection over a
  larger family (amplitudes, several phases) is extremal in all directions and the identity T = 6E′ then holds along
  the L-rotation direction only, which is still enough for zero torque of the extremal state.
- The CP-violating test is a grid scan (721 points); the torque at the minimum is zero up to grid spacing, not to
  machine precision.
- The clock is the model's; nothing about physical cosmological time, dressing or clock accuracy (FORMAL_CORE COND
  2.20) is addressed.

## 9. Reproduction (own code, 45 s; the Codex rerun 4 s)

```
cd research/reality_framework/grut_formation_review_20261025/tracks/rp_preparation/code
python3 -B neutral_boundary.py ../results/neutral_boundary.json
python3 -B projection_identity.py ../results/projection_identity.json
# labelled reproduction of Codex's code, in a scratchpad copy of the packet (git archive e160356):
# GRUT_REHEATING_WORKDIR=<scratch> python3 -B <scratch>/.../code/boundary.py   -> max |diff| vs committed JSON 3.4e-15
```
