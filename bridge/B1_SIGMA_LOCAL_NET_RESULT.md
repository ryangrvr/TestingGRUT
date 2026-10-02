# B1 — Σ / LOCAL-NET RECONSTRUCTION (result)

**Question (owner).** If the GRUT local net {𝒜_x} is **deleted from the input**, can the remaining GRUT structure
(generator, response structure, noise, drift) reconstruct it?

**Scripts:** `bridge/b1/b1_net_reconstruction.py` (+ `.log`) and `bridge/b1/b1_rigidity.py` (+ `.log`), with shared helpers in
`bridge/b1/b1_lib.py`.
- The canonical models are **re-implemented** from the canonical charters at `b935099`, not imported.
  - C1-a: K = 0.3·I + L(path), N = 24.
  - L0-1b: K = 0.3·I + L(w), w ≥ 0, with its geometry predicate (Q ≥ N/2 and R₁ⱼ monotone).
  - L0-1c: on-site convex quartic.
  - L0-1e: Q = 2·diag(Tᵢ).
- **"Net deleted"** means the operator is handed over in a Haar-random orthogonal frame. Only basis-free data survive.

**Firewall (no circular criteria).**
- No criterion that names sites is allowed. "k-local", "on-site", "nearest neighbour **in the canonical labelling**" and
  "retained site 0" are all excluded.
- Admitted criteria are frame-free: sparsity, bandwidth-1 / connectedness, autonomy, response-kernel factorization,
  invariant subspaces, passivity (pins ≥ 0), and the declared GRUT class constraints (uniform pin, unit springs).
- Every admitted criterion is recorded as a **price** where it does the work.

## Results by route

| route | test | finding | compression class |
|---|---|---|---|
| **sparsity / autonomy / response factorization** | sparsest frame of K | the normal-mode frame makes K **diagonal** (nnz 24 < canonical 70). This is the optimum of every factorization criterion, and it has **no edges and no geometry** | **NO BRIDGE** (the canonical net is not the optimum) |
| **connected sparsest (chain)** | Lanczos from a unit start vector q | **every** q gives a tridiagonal (nearest-neighbour chain) frame, so the chain frames form a **continuum** (S^{N−1}/sign) | **NO BRIDGE** without a further price |
| **+ passivity** (spring-network reading, all pins ≥ 0) | q = canonical end site + ε·noise | passive chains: 50/50, 38/50, 11/50 and 1/50 at ε = 0.02, 0.05, 0.1 and 0.2. **Every one passes the L0-1b geometry predicate and none is signed-permutation-equivalent to the canonical chain.** Random q: 0/400 (the admissible set is a neighbourhood, not the sphere) | **NO BRIDGE**: a continuum of inequivalent passive nets |
| **dimension** | weighted 4×4 grid generator (canonical net 2D; line test FAILS, Q = 3.85) | a **passive** 1D chain frame of the *same operator* exists (min pin +0.1475; chain Q = 40.5, monotone, **SURVIVES**) | **NO BRIDGE**: not even the dimension of the net-derived geometry is fixed by K |
| spectral multiplicity | uniform 5×5 grid | 11 degenerate pairs, so it is **never** a connected chain (a Jacobi matrix has a simple spectrum) | a **partial** constraint only |
| **L0-1b class** (uniform pin 0.3, weights ≥ 0) | rotations fixing the uniform vector, minimizing positive off-diagonals | **N = 12:** an exactly feasible, **inequivalent**, isospectral weighted Laplacian (1 hit in 24 restarts). Details below. **N = 24:** not found in 6 restarts (**not established**) | **NO BRIDGE** at N = 12 (class does not pin the net); open at N = 24 |
| **narrow C1-a class** (uniform pin, **unit** springs) | known theorem: the path is determined by its Laplacian spectrum (DLS) | the net is unique up to automorphism and the K-commutant (H-relative gauge) | **CONDITIONAL COMPRESSION, CRITERION-PRICED**: unit springs plus "graph Laplacian" *are* the supplied net class |
| **invariant sectors** | eigenspaces of a nondegenerate K | global modes only; none localized | **NO BRIDGE** |
| **memory kernel / retained partition** | k(τ) of the canonical site vs a random-q chain "site 0" | both are well-formed, decaying kernels (3.9e-1 → 5e-7 vs 2.9e-1 → 2e-7). **Every** unit vector is the retained site of some chain net | **NO BRIDGE**: the retained/bath partition is not fixed by K |
| **noise layer** (L0-1e) | eigenframe of Q′ | uniform Tᵢ: Q′ ∝ I, so **no** information. Non-uniform Tᵢ: recovers the site axes exactly (min overlap 1.000000) | **RELOCATION** (the net is carried by the supplied environment layer) |
| **nonlinear drift** (L0-1c on-site quartic) | odeco tensor decomposition (power iteration + deflation) | recovers all 24 site axes exactly (overlap 1.00000000; residual 4e-13). Odeco decompositions are unique, so the drift **determines** the net. **Control:** a quartic written on-site in another frame recovers **that** frame (overlap with the canonical sites 0.587) | **RELOCATION / REDUNDANT SUPPLY** + **CONSISTENCY** (drift and net agree) |

**The N = 12 L0-1b counterexample** (`b1_rigidity.log` iii–iv):
- K′ = 0.3·I + L(w′) is orthogonally conjugate to the canonical path K, with the uniform vector fixed.
- w′ ≥ 0, and all **66/66** pairs carry positive weight, in the range [1.29e-5, 1.52].
- Row sums are 1e-15 and the spectrum difference is 2e-14.
- It is **not** signed-permutation-equivalent to the path.
- Its geometry predicate **FAILS** (Q = 34.1, but R₁ⱼ is not monotone).
- So inside the stated L0-1b class, the generator does not determine the net, and the net-derived geometry differs.
- **Caveat:** the smallest weights are tiny (1e-5). The example is exact but near the boundary of the cone.
- A dimension count makes a large family plausible: the weight cone has dimension N(N−1)/2, against N−1 spectral
  constraints. This count is **not** a proof.

## Verdict

> **Σ (the GRUT local net) REMAINS SUPPLIED in GRUT's earned linear core.**
> - It is **CONSISTENCY-CHECKED / OVERDETERMINED** where GRUT also supplies site-local nonlinear drift (L0-1c) or
>   non-uniform site noise (L0-1e). There the net is recoverable from those layers. That is **RELOCATION, not derivation**:
>   the information was written into them in site-local form, and the control shows the drift selects whatever net it was
>   written in.
> - It is **uniquely reconstructable only in the narrow C1-a class** (unit springs, uniform pin). That is a
>   **CONDITIONAL COMPRESSION** whose condition *is* a supplied net class (**CRITERION-PRICED**).

**Not Σ REDUNDANT IN GRUT.** The owner's elimination test fails: the net cannot be reconstructed without supplying it,
or supplying a layer that encodes it.

**Firewall check.** No reconstruction used site labels.
- The two positive recoveries (drift, noise) use frame-free procedures: tensor decomposition and an eigenframe.
- Their *inputs* were supplied site-locally. That is exactly why they are classed RELOCATION.

**Relation to SCOUT.** This matches SCOUT's reviewed Σ wall, with a **category difference** (crosswalk row 1):
- GRUT's earned substrate is classical and direct-sum, so the "net" is a frame plus edge weights.
- SCOUT's Σ is a tensor-product factorization.
- The GRUT-side reconstruction problem is **weaker** (fewer invariants to recover) and **still fails**. That is a negative
  for "GRUT-specific structure does extra work".

## Proposed canonical note (proposal only)

> **CANONICAL UPDATE CANDIDATE (B1-CUC-1):** "The Level-0 local net is not determined by the linear generator K, even
> within the L0-1b class (explicit isospectral inequivalent counterexample at N = 12, `grut-bridge-1`
> `bridge/b1/b1_rigidity.py`). Where the on-site quartic drift (L0-1c) or non-uniform site temperatures (L0-1e) are
> supplied, the net is **overdetermined** by them (recoverable by odeco decomposition / noise eigenframe). The net remains a
> supplied layer; the drift and noise layers carry redundant copies of it."
>
> This is consistent with, and sharpens, EA0_OWNER_RULING_02 §3 ("existence or uniqueness of the local net itself" not
> derived). **Not applied; canonical GRUT untouched.**

## Scope and limits

- The numerics are a single code path by the same agent: an **independent code path, not an independent reviewer**.
- N = 12 / 24 / 16 / 25 only. The N = 24 L0-1b case is open.
- The quantum (tensor) version of the question is SCOUT's S2-Σ / S2-1b, and was **not** redone here.
- The power-law L0-1b canonical weights (α family) were not separately searched. The path is used as the α → ∞ member of
  the class.
