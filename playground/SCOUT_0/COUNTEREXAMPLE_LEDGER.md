# SCOUT_0 COUNTEREXAMPLE LEDGER

One decisive counterexample per broad claim attempt. Smallest first.

## CE-01 (P-06): "CM-ness of a kernel is isolated in its deformation space"

- **Claim attacked:** that the S5-1 parent's kernel might be CM-isolated (the charter's null).
- **Smallest counterexample:** the family S_p(x) ∝ (x²+κ²)^{-p} at p=2: kernel f(t)=e^{-κt}(t+1/κ), f''(t)=e^{-κt}(κ²t−κ)<0 on (0,1/κ). Exact, artifact-free, CM fails on an open interval.
- **Nearest counterexample:** p=1.05 — non-CM detected at t=0.08 (order 6), nothing at t≥0.5. Boundary location between 1.0 and 2.0 is bracketed, not proven.
- **Same assumptions, opposite result:** p<1 members are all CM at tested resolution — the opposite side of the same family.
- **Same result, fewer assumptions:** none needed beyond the recorded κ and GR transform table.

## CE-02 (P-06b audit): "if g is CM and α ∈ (0,1], then z^α g(z) is CM" (lemma used in `81ad661`)

- **Claim attacked:** the lemma carrying the `1/2 < p ≤ 1` leg of the banked P-06b proof.
- **Smallest counterexample:** g ≡ 1 (CM, Bernstein measure δ₀); z^α is increasing, not CM.
- **Nearest counterexample:** g = e^{-z}, α = 1: (z e^{-z})′ = e^{-z}(1−z) > 0 on (0,1). The script's integral form of the lemma is also false (RHS/LHS = 1.5 at α=1/2, s=1, z=1; exact RHS = e^{-sz}[s z^{α−1} + (1−α) z^{α−2}]).
- **Same assumptions, opposite result:** none — the lemma is simply false; the specific conclusion it was used for (z^{2ν}·ℒ[(s²−1)^{ν−1/2}] is CM for 0<ν≤1/2) is nevertheless true, because that function equals z^νK_ν(z), which is CM by the unified Bernstein representation.
- **Same result, fewer assumptions:** the unified representation z^νK_ν(z) = √π2^ν/Γ(1−p)·∫_1^∞e^{-zs}(s²−1)^{-p}ds proves CM on all of 0<p<1 with no lemma at all (DLMF 10.32.8 at index 1/2−p). See `probes/P06B_CORRECTION_01.md`.
- **Lesson for the scout:** every numerical check in `81ad661` tested the (true) conclusion, never the lemma. Lemmas introduced to bridge a gap need their own smallest-counterexample pass before banking.
