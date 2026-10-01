# SCOUT-1 W5-ES RESULT — emergent symmetry (the C3 hostile against TC-4′)

**Charter:** `PROBE_CHARTERS.md` §W5-ES (pre-registered).

**Files:** `w5es_emergent_symmetry.py`, with log `w5es_emergent_symmetry.log`.

**Labels:**
- **A. C3-PARTLY-ERASABLE in principle**, together with **B. C3-SUPPLIED-IN-GRUT-SCOPE**;
- KNOWN RESULT IMPORT: Aharony cubic-anisotropy RG; emergent rotation invariance at 2D Ising criticality — STANDARD ✓.
  The N_c ≈ 2.9 estimate is SECONDARY.

## 0. Verdict

> **C3 is partly erasable: RG can forget symmetry-*breaking detail* and produce `S_IR ⊋ S_UV`.**
>
> - **O(N) + cubic anisotropy:** the cubic coupling v has `y_v = ε(N−4)/(N+8)`, which is −0.20 / −0.09 / −0.04 at
>   N = 2 / 3 / 3.5. Flows from different UV values of v (+0.15, +0.30, and for N = 2 also −0.10) all converge to
>   **exactly** (u*, 0): O(N) emerges from a theory with only cubic symmetry.
> - **2D Ising lattice:** the microscopic symmetry is C₄ (the square lattice). At T_c the axis/diagonal correlation
>   ratio at equal |r| is 1.037 → 1.010 → 1.003 → 1.0018 → 1.0007 for |r| = 1.4 … 22.6. This decays ≈ r⁻², as
>   expected for the irrelevant lattice-anisotropy operator (y = −2). **SO(2) emerges.**
>
> **What is forgotten:** the *magnitudes and signs* of the symmetry-breaking couplings (v; the lattice-anisotropy
> coefficients). This is C1-type information attached to a C3 deformation.
>
> **What is NOT forgotten, i.e. still supplied:**
> - **N, the field content (C4).** It decides whether enhancement occurs at all (`N < N_c`) and *which* larger group
>   emerges (O(N)).
> - **The dimension (C5):** the sign of `y_v` depends on ε.
> - **The basin.**
>   - For N ≥ 3 a UV start at v = −0.10 **runs away** (a fluctuation-induced first-order transition): no fixed point,
>     no emergent symmetry.
>   - Criticality itself is a tuning (W1-R).
> - **Contested marginality.** One loop gives N_c = 4. Best estimates put N_c ≈ 2.9 (secondary). So at N = 3 the cubic
>   perturbation is in fact weakly relevant, and enhancement is **not robust**.
>
> **Hostile control: redundant ≠ erased.** With uniaxial anisotropy `J_y/J_x = 0.5` at its own T_c, `C_x/C_y` at equal
> r stays at 1.125 → 1.111 → 1.098 → 1.087 (r = 2 … 16). That is far from the isotropic convergence (0.07 % at r ≈ 23).
> The anisotropy is a **redundant** operator (a coordinate rescaling). Rotation symmetry emerges only in rescaled
> coordinates, and the aspect ratio survives as geometry. The slow drift is a finite-size effect (r up to L/8), not
> extrapolated.
>
> **Not computed (cited, SECONDARY):** a dangerously irrelevant anisotropy (3D XY with Z₆: an emergent U(1) at
> criticality, broken again at a longer scale in the ordered phase). It is a further case where emergence is
> scale-limited.

## 1. Information accounting (W5-3)

| Item | UV specification | IR specification |
|---|---|---|
| symmetric couplings (u) | value | forgotten (fixed point u*) |
| symmetry-breaking coupling (v; lattice anisotropy) | value + sign | **forgotten** (→ 0) *inside the basin* |
| microscopic symmetry group | cubic / C₄ | **replaced** by O(N) / SO(2) |
| N (field content), d | supplied | **still supplied**; it decides whether and which enhancement |
| basin membership (sign of v for N ≥ 3; critical tuning) | — | **still supplied** (the state / tuning) |
| redundant anisotropy (uniaxial) | ratio | **kept** as geometry (not erased) |

Classification: **SYMMETRY-BREAKING DETAIL FORGOTTEN** (scoped). This is **not** full symmetry selection: the enlarged
group is fixed by C4/C5 content and the basin.

## 2. W5-5 GRUT connection

**Do any admitted GRUT variables wash out like irrelevant symmetry-breaking deformations?**

| GRUT datum | Behaviour | Classification |
|---|---|---|
| coupling anisotropy in K (W1-S: `J_y/J_x = 0.2`; `1 : 0.5 : 0.1`) | γ unchanged; amplitudes change | **redundant** (Gaussian: absorbed into a coordinate rescaling), **not erased**. Same as the uniaxial control |
| lattice-shape detail of K (NNN with a quadratic minimum, W3-TC3) | IR γ unchanged | **irrelevant** analytic dispersion data forgotten. Gaussian-class erasure of C1/C2 data, not of a symmetry group |
| site-dependent rates / random couplings | not shown to wash out. Spectral data at the retained site depend on them (EDA-01) | not erased |
| sector distinctions (SF-1) | define *different* IR classes (z = 1 vs 2) | not erased |
| lift differences (P-08/P-09) | invisible in observables | **quotient invisibility** (a representation datum not observed), not an RG flow to enlarged symmetry |
| noise-profile asymmetry | **generates** operators once a nonlinearity exists (W4) | the opposite of erasure |
| boundary / readout | changes γ by an integer (W1-S) | not erased |

**Answer:** no GRUT-supplied **C3** datum is known to wash out under the earned flow. The earned flow is Gaussian. It
erases irrelevant analytic C1/C2 data and absorbs anisotropy redundantly, but it has **no interacting fixed point
at which a symmetry-breaking operator could become irrelevant**. Hence **B within GRUT scope**. Lift invisibility is a
quotient effect and was already banked in SCOUT-0.

## 3. TC-4′ wording (W5-7)

> **RG cannot violate exact microscopic symmetry constraints (absent anomalies), but it can erase symmetry-breaking
> deformations and produce an enhanced IR symmetry within a basin. What remains upstream is the field /
> representation / anomaly content (which decides whether, and to which group, enhancement occurs), the dimension,
> and the data specifying the fixed-point basin (tuning, state, sign of the breaking coupling). Redundant deformations
> are not erased but absorbed into geometry.**

**Status: W5-ES COMPLETE — A (C3-PARTLY-ERASABLE, in principle: breaking detail forgotten; the enlarged group set by
C4/C5 + basin) with B (no GRUT-internal mechanism: Gaussian floor, no interacting fixed point).**
