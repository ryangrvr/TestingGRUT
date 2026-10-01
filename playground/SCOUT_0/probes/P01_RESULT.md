# SCOUT_0 W3 P-01 RESULT — do the E-1 floor edges port off the C1-a chain?

**Charter:** `PROBE_CHARTERS.md` P-01 (frozen): test E-1's gap → P_memory and passivity → P_positivity
on different coupling families with explicit spectral arguments. **Record read:**
`L0_1A_CHARTER_01.md` (frozen gates G-1…G-3, P-1…P-3; the L0-1a comparator),
`L0_1A_OWNER_RULING_01.md` (the scope: "within the declared C1-a construction class"; the exact
factorization `k_pin = e^{−pin·τ}k₀`). **Script:** `p01_portability.py` (log `p01_portability.log`).
NEW HYPOTHESIS CLASS — DOES NOT ALTER OLD TERMINAL. E-1 stands at its recorded scope.

## 0. Verdict

> **Neither edge ports as stated. Each ports as a statement about the retained site's local spectral
> measure `μ_r`, and the C1-a chain is a family where the stated form and the spectral form coincide.**
>
> - **Gap → memory.** The portable theorem is **`inf supp μ_r > 0 ⇒ exponential-grade memory`** (for the
>   symmetric spectral representation `k(τ) = ∫e^{−λτ}dμ_r`). What is amenability-specific is only
>   **whether the supplied pin is load-bearing**: on the amenable chain, removing the pin lets `μ_r`
>   reach 0 (algebraic memory); on the non-amenable infinite graph/tree class a positive spectral
>   bottom exists without a pin, so the pin is not load-bearing. On **non-amenable** bath graphs (random d-regular, d ≥ 3, converging to the
>   d-regular tree) memory is exponential **with no pin at all**: the soft mass seen by the retained
>   site vanishes like `≈ 0.36/n` (d = 3), and the rest decays at a rate at or above the Kesten–McKay
>   bottom `d − 2√(d−1)`. On the amenable chain the soft mass is n-independent (0.013) and the tail is
>   algebraic. **The theorem ports everywhere; the pin's necessity is family-specific** (amenable vs
>   non-amenable). *Wording repair 01 (external audit): earlier summaries said "gap → memory ports only
>   to amenable baths"; that conflated the portable theorem with the load-bearing status of the pin.*
> - **Passivity → positivity:** **sufficiency** ports to every symmetric family (all `λ ≥ 0` ⇒ CM).
>   **Necessity** ports only where the retained site sees every mode (the chain does: Jacobi
>   eigenvectors never vanish at the end). A symmetric two-branch bath with an unstable
>   **antisymmetric** mode is non-passive yet has an exactly CM retained kernel. Under **asymmetry**,
>   sufficiency fails too: a directed ring is accretive (`min eig(K+Kᵀ) = 0.15`) but its kernel is
>   non-monotone (complex modes; the record's E-4).

Charter hypothesis ("gap → memory transfers; passivity → positivity transfers only for symmetric
Laplacians"): **gap half: the theorem `inf supp μ_r > 0 ⇒ exponential memory` transfers to every symmetric family; the pin is load-bearing only on amenable families**; **passivity half
TRUE for sufficiency, FALSE for necessity** (symmetric families with invisible modes). Charter null
("neither transfers"): false as stated; both transfer in their spectral forms.

## 1. Gap edge

**Exact identity (L0-1a ruling):** `K_b = pin·I + L_b`, so `k_pin(τ) = e^{−pin·τ}k₀(τ)` for **any**
bath graph. Deleting the pin therefore matters exactly when `k₀` is not already exponential, i.e.
when `μ_r` (the local spectral measure of `L_b` at the coupling vertex) has mass reaching down to 0.

| Bath family (pin = 0) | `n` | soft mass below threshold | raw window grade | soft-removed kernel |
|---|---|---|---|---|
| chain (C1-a, amenable) | 23 / 95 / 383 | **0.013 / 0.013 / 0.013** (n-independent) | ALGEBRAIC | — (the tail is real) |
| random 3-regular | 200 / 800 / 3200 | `1.9e-3 / 4.5e-4 / 1.1e-4` (`n·mass` = 0.374, 0.363, 0.361) | ALGEBRAIC (finite-n plateau) | **EXPONENTIAL**, rate 0.26–0.28, min visible `λ` → 0.175 (bottom 0.1716) |
| random 4-regular | 200 / 800 / 3200 | `n·mass` = 0.53 | ALGEBRAIC (finite-n plateau) | **EXPONENTIAL**, rate 0.64–0.66, min visible `λ` → 0.543 (bottom 0.536) |

**Limit argument.** Random d-regular graphs converge locally to the d-regular tree (Benjamini–Schramm;
McKay). The tree Laplacian has spectrum `[d − 2√(d−1), d + 2√(d−1)]` (Kesten–McKay). The coupling
spring adds `+1` at the root. That is a Loewner-positive perturbation, so it cannot create spectrum
below the bottom. Hence in the limit `k₀(τ) ≤ e^{−(d−2√(d−1))τ}` — exponential-grade memory with no
pin. The finite-n plateau is a single near-constant mode of weight `O(1/n)`. Note the order of
limits: the raw window grade at any finite n is ALGEBRAIC; the exponential statement is the
`n → ∞` local limit.

**General form (theorem-grade, standard).** For a bounded-degree vertex-transitive bath graph, the
pin-free bottom of spectrum is > 0 **iff the graph is non-amenable** (Kesten 1959; Dodziuk 1984 via
the Cheeger constant). So, in that class: *gap-deletion destroys exponential memory iff the bath is
amenable.* The C1-a chain is amenable, which is why E-1 certified the pin as load-bearing.

## 2. Passivity edge

**Exact statement (symmetric `K_b`):** `k(τ) = Σ_k u_k² e^{−λ_k τ}` with distinct exponentials, so the
Bernstein measure is unique. Therefore **k is CM ⇔ every λ_k with `u_k ≠ 0` satisfies `λ_k ≥ 0`.**

- **Sufficiency** (passive ⇒ CM): every symmetric family.
- **Necessity** holds iff every unstable mode is visible (`u_k ≠ 0`). On a Jacobi (chain) matrix every
  eigenvector has a nonzero end component (min retained weight `3.8e-4 > 0`), so necessity holds there.
  This is why E-1 certified it.
- **Counterexample off the chain:** bath `{h, a, b}`, springs `h–a = h–b = 1`, `a–b = −1.5`, pin 0.3.
  Eigenvalues `{−1.7, 0.568, 4.03}`: not passive. The unstable mode is antisymmetric in `a ↔ b`, with
  retained weight `0.0`. The kernel is strictly decreasing and CM, with `k(40) = 2.9e-11`. Breaking the
  symmetry by 0.05 makes the mode visible (weight `2.9e-5`), and the kernel grows (`k(40)/k(1) = 2.7e25`).
  So necessity is **generic but not structural**: it fails on the measure-zero set where a symmetry
  fixing the retained site hides the unstable mode.

**Asymmetric (directed) Laplacian:** a directed ring with pin 0.05 has `min eig(K+Kᵀ) = 0.150`
(accretive, i.e. passive in the dissipative sense) and complex eigenvalues. The kernel stays positive
but is **not monotone**, hence not CM. Passivity in the accretive sense no longer implies positivity
of the CM kind. This is the record's E-4 (cycle affinity ⇒ ¬CM), reached from the portability side.

## 3. What carries the floor edges

Both edges factor through one object: the retained site's local spectral measure `μ_r`, plus symmetry
of `K_b`.

| Edge | Carrier | C1-a chain | Off-chain |
|---|---|---|---|
| gap → memory | `inf supp μ_r > 0` | supplied by the pin (amenable graph) | can be supplied by the graph (non-amenable) |
| passivity → positivity | sign of the `μ_r`-visible spectrum; symmetry | every mode visible (Jacobi) | invisible unstable modes; asymmetry breaks it |

**Classification:** REDISCOVERED-KNOWN as mathematics (Kesten–McKay; Kesten/Dodziuk amenability;
Bernstein uniqueness; Jacobi endpoint non-vanishing). KNOWN-BUT-NEW-IN-GRUT: E-1's two certificates
are properties of the chain's spectral measure, not of "gap" and "passivity" as such.

## 4. Feed to the zoom-out

This is the **quotient principle again**. The retained kernel factors through `μ_r` (consistent with
E-7/C-7: geometry is underdetermined at single-site access), so E-1's edges are statements about
`μ_r`. P-02 found the same pattern on the sector side: the IR law class factors through the
occupation-edge velocity `v_b^∞`, a **local spectral datum at an edge**. Two independent probes now
say that **the law-level invariants live in local spectral data at an edge** — the bottom of `μ_r`, and
the dispersion at `k_b^∞` — and that whatever supplies that edge (pin vs non-amenability; exclusion
vs interaction) is fibre data.

**Status: P-01 COMPLETE (one pass). `inf supp μ_r > 0 ⇒ exponential memory` ports to every symmetric
family; the pin is load-bearing only on amenable baths (non-amenable baths give pin-free exponential
memory). Passivity edge: sufficiency ports to symmetric families; necessity only
where the retained site sees every mode; asymmetry breaks sufficiency (E-4). Both edges factor
through the retained local spectral measure. Next: P-24 lock census.**
