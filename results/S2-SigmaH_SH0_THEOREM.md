# S2-ΣH ΣH-0 — theorem-first: state independence alone selects no TPS

**Status:** TRUE DERIVATION (elementary), plus a dimension count. Scope: a finite-dimensional ℋ ≅ ℂ^N with N = Π d_i.

## Proposition 1 (productness no-go)

For any unit vector |ψ⟩ ∈ ℂ^N and any factor dimensions with Π d_i = N, there is a unitary W with
W|ψ⟩ = e^{iφ}|0⟩⊗…⊗|0⟩. In the TPS defined by W, ψ is a product state: every factor entropy and every total correlation is
zero.

**Two sets, distinguished (REPAIR 06, Y-13).**

1. **Fixed target ray.** For a *fixed* product ray |0…0⟩, the set {W : W|ψ⟩ ∝ |0…0⟩} is the coset {V W₀ : V|0…0⟩ ∝ |0…0⟩}
   ≅ U(1) × U(N−1), of real dimension **(N−1)² + 1**.
2. **The full product-making family.** F_ψ = {W : W|ψ⟩ is *some* product state} also ranges over the manifold of pure
   product rays, of real dimension 2 Σ_i (d_i − 1). Hence

       dim F_ψ = (N−1)² + 1 + 2 Σ_i (d_i − 1),
       codim F_ψ in U(N) = 2N − 2 − 2 Σ_i (d_i − 1),

   which is the codimension used in Proposition 2.

*Proof.* Take a Householder reflection W₀ sending ψ to |0…0⟩ up to a phase, and compose with the stabilizer: this gives
set 1. Set 2 is the union over product rays r of the analogous cosets, fibred over the product-ray manifold. ∎

**Consequence.** Minimizing state correlation over unrestricted TPSs has minimum 0, attained on a **high-dimensional,
non-unique family** of frames. → **STATE INDEPENDENCE ALONE → NO TPS SELECTOR.** The same holds for "initial low entropy", which is
zero for every factor in such a frame.

## Proposition 2 (where the information has to come from)

**Joint compatibility.** Fix a locality class (k, d, n). Let F_H be the set of frames in which H is k-local, and F_ψ the
set of frames in which ψ is product.

- By CPR, when H is k-local at all, F_H is generically a finite union of LU × Sym orbits.
- F_ψ is an LU-orbit-saturated set of codimension (2N − 2) − Σ_i (2d_i − 2) inside the frame manifold, which is the
  codimension of product states among pure states.
- So **F_H ∩ F_ψ ≠ ∅ is a non-generic condition on the pair (H, ψ).** For instance, for Haar ψ it has probability 0 given
  H, measure stated.

**Dichotomy (classification of what a joint selector can do):**

| case | condition | what a joint selector gives |
|---|---|---|
| **(a) compatible pair** | F_H ∩ F_ψ ≠ ∅ | the joint condition is satisfied exactly, and the selected class is CPR's (unique up to LU × Sym, **given k, d, n**). The joint principle adds no selection power beyond H-locality. Its content is the *statement that ψ happens to be product there*: a boundary datum, H_corr |
| **(b) incompatible pair** (the generic case) | F_H ∩ F_ψ = ∅ | there is a genuine trade-off between H-locality and ψ-independence: a Pareto family. Choosing in it needs a preference / threshold / weight |

Either way, the local dimension, the number of factors and the locality class enter as in S2-Σ / CPR.

The numerics below test whether finite candidate families follow this dichotomy, and whether anything (epoch, arrow,
records, MDL) collapses case (b).
