# S2-G2 / G3 / G4 RESULT — the dimension origin campaign

**Charter:** `probes/PROBE_CHARTERS.md` §S2-G2/G3/G4. Pre-registered at `e401bbe`, before any run.

**Files** (`probes/S2-G/`):

| script (+ `.log`) | covers |
|---|---|
| `s2_g2_graph.py` | graphs, node-averaged estimators |
| `s2_g2_exact.py` | exact product / circulant spectra; local estimators |
| `s2_g3_causal.py` | causal sets, lattice spacetimes, propagation cones |
| `s2_g4_capacity.py` | capacity: LPs, noiseless subsystems, access, noise, composition |

**G1** (local Hilbert d) was already adjudicated: **Σ PRIMITIVE**.

## Headline

> **"Dimension" is not one primitive.** Each notion tested is a **derived observable of a different parent structure**:
>
>     d_g, d_s, d_w  = f(interaction graph, probe dynamics, scale, site/averaging)
>     d_causal       = g(causal order or propagation rule, region shape, estimator)
>     N (capacity)   = h(state space, effect/access set, dynamics/noise)
>     d (local)      = Σ (supplied)
>
> **None determines another.** Every cross-hostile the owner asked for was found:
> - G2 same / G3 different;
> - G3 same / G4 different;
> - G4 same / G1 and representation different.

## G2 — graph / spectral dimension

**G2-1 controls.**
- Exact product spectra, L = 4000: d_s(t = 1000) = **1.000 / 2.000 / 3.000 / 4.001** for ℤ¹ … ℤ⁴.
- On the 41³ torus: d_g → 2.89 (r = 15, lattice correction) and d_w = 1.97 – 2.07.
- On the graph-based ring and 2D torus:
  - ring: d_g 0.96, d_s 1.02, d_w 2.00;
  - 45² torus: d_g 1.92, d_s 2.12 – 2.19, d_w 2.00.

**G2-2 same order of N, different connectivity:**

| graph | d_g | d_s | d_w | notes |
|---|---|---|---|---|
| random 3-regular | 1.55 → 3.41 (r = 2 → 5), growing | 1.9 → 14 (t = 2 → 100) | 1.5 → 15 | no finite dimension (exponential growth) |
| binary tree | 0.95 → 1.88, growing | 1.1 – 1.7 | 2.4 → 3.3 | no stable value |
| small world p = 0.05 | ≈ 1 | ≈ 1.05 up to t = 100 | — | p = 0.1: 1.02 → 1.25 by t = 300 (crossover) |
| **Sierpinski gasket** | **1.584** (corner, r = 2⁷; theory 1.585) | **1.32 – 1.43** (local; theory 1.365, log-periodic) | **2.26 – 2.32** (theory 2.322) | **three different numbers on one graph** |
| **comb, backbone site** | **1.98** (theory 2) | **1.51** (theory 3/2) | — | same growth as the square lattice, different diffusion |
| comb, mid-tooth site | 0.995 | 1.000 | — | |
| comb, node-averaged | — | ≈ 1.00 | — | **site / averaging choice changes the answer** |

**G2-3 crossovers** (exact; ε = 0.01):

| system | d_s(t) |
|---|---|
| 2D bundled chains | 1.17 (t = 3) → 1.86 (t = 30) → **2.00** (t ≥ 1000) |
| 3D layered | 2.22 (t = 3) → 2.86 (t = 30) → **3.00** (t ≥ 1000) |

→ **SCALE-PRICED.**

**G2-4 diffusion-operator hostile** (same vertex set):

| operator | d_s | theory |
|---|---|---|
| random weights U(0.5, 1.5) | 2.12 – 2.18 | unchanged from standard 2.12 – 2.19 |
| anisotropic | 1.17 → 2.16 at fixed t = 5 → 100 | — |
| fractional Δ^{a/2}, a = 1 | **4.02 – 4.03** | 4 |
| fractional Δ^{a/2}, a = 0.5 | **8.13** at t = 10 | 8; finite-size beyond |
| 1D Lévy, α = 1.5 | 1.36 – 1.43 | 1.33 |
| 1D Lévy, α = 1 | 1.94 – 2.02 | 2 |
| 1D Lévy, α = 0.5 | 3.44 at t = 10 | 4; finite-size beyond |

→ **G2 SPECTRAL DIMENSION PROBE-DYNAMICS- (D-) PRICED.**

**G2-5 structural.**
- d_s(t) depends **only** on the normalized-Laplacian spectrum.
- Among all 853 connected 7-node graphs there are **12 cospectral non-isomorphic classes, all 12 with different ball-growth
  profiles**. Example: the star K₁,₆ vs K₂,₅. Their d_s(t) is identical at every t while their growth differs.
- The converse also holds: comb vs square lattice have the same growth and different d_s.
- → d_s and d_g each see structure the other cannot. **G2-NONUNIQUE as a geometry selector.**

**G2 verdict:** **G2-DERIVED-FROM-D** (from the graph plus a probe walk), and **SCALE-, PROBE-DYNAMICS- and
SITE / AVERAGING-PRICED**. The graph connectivity itself remains a D / Σ input (S2-1b / S2-Σ).

## G3 — spacetime / causal dimension

**G3-1 / G3-2: sprinkled Alexandrov intervals** (4 samples each):

| true d | N | Myrheim–Meyer | midpoint scaling | longest-chain scaling (500 → 2000) |
|---|---|---|---|---|
| 2 | 500 / 2000 | 2.009 / **1.998** | 2.027 / 2.000 | **1.948** |
| 3 | 500 / 2000 | 2.984 / **2.996** | 3.065 / 3.044 | **2.678** |
| 4 | 500 / 2000 | 3.984 / **4.014** | 4.147 / 4.129 | **2.974** |

- **CAUSAL DIMENSION DERIVED FROM ORDER**, given an interval-shaped region. MM converges fastest.
- Chain scaling is strongly biased at these N. → **ESTIMATOR-PRICED.**
- **Shape hostile (a slab instead of an interval):** MM gives **3.47 / 5.45 / 7.37** for true 2 / 3 / 4. The region shape is a
  supplied input.

**G3-3 same spatial graph, different causal structure** (lattice spacetimes, interval between tips):

| spatial graph | causal rule | MM d |
|---|---|---|
| ℤ | L1 cone, c = 1 | 1.958 |
| ℤ² | L1 cone, c = 1 | 2.832 (T = 10), **2.902** (T = 16) |
| ℤ² | L∞ cone | 2.854 (T = 10), 2.897 (T = 14) |
| ℤ² | L1 cone, c = 2 | 2.889 |
| **ℤ² (same graph)** | **absolute time (instantaneous propagation)** | **1.142** (midpoint 1.220) |

**The same spatial graph (d_g = d_s = 2) yields causal dimension ≈ 3 (finite-speed cones) or ≈ 1 (absolute time).**
→ **GRAPH DIMENSION ≠ SPACETIME DIMENSION.** The extra input is a **finite-speed propagation rule** plus the region.

**G3-4 / G3-5: propagation cones** from hopping Hamiltonians (tail threshold 10⁻¹⁰; the threshold-defined speeds exceed the
group velocity):

| lattice | cone speed | cone-volume exponent |
|---|---|---|
| chain | 3.00 | 0.95 |
| chain + NNN | 5.00 | 0.93 |
| square | 2.70 axis vs 3.75 diagonal (**anisotropy 1.39**, ≈ √2) | 1.91 |
| triangular | 3.90 / 5.52 / 7.57 (**anisotropy 1.94**) | 1.97 |

- The cone volume gives back the **graph** dimension (≈ 1, 2). It adds a causal structure, but a **direction-dependent** one
  with a preferred lattice frame.
- Different lattices with the same dimension have different speeds and anisotropies. The low-energy dispersion is
  quadratic (non-relativistic).
- → **DIMENSION SELECTED (= graph) but LORENTZ STRUCTURE NOT DERIVED.**
- In sprinkled causal sets, Lorentz invariance is **inherited from the Minkowski parent of the sprinkling** (priced; KNOWN
  RESULT IMPORT: Bombelli–Henson–Sorkin, SECONDARY).

## G4 — operational / information capacity

**G4-1 matched capacity, different state spaces** (LP over vertex subsets for polytopes):

| system | N | K | Hilbert d |
|---|---|---|---|
| classical bit | 2 | 1 | — |
| rebit | 2 | 2 | 2 (real) |
| gbit (square) and 5- … 8-gons | **2** | 2 | — |
| qubit | 2 | 3 | 2 |
| spin factor in ℝ⁴ | 2 | 4 | — |
| triangle = classical trit | 3 | 2 | — |
| qutrit | 3 | 8 | 3 |
| real d = 3 | 3 | 5 | 3 |

→ **CAPACITY DOES NOT FIX REPRESENTATION.**

**G4-2 capacity from dynamics.**
- A Markov chain with 3 closed classes: asymptotic capacity **3** (from 7).
- Collective SU(2) noise on 3 qubits: commutant dimension 5, centre dimension 2, block sizes [4, 4], multiplicities
  **(2, 1)**. The noiseless capacity is **3** (from 8).
- → **CAPACITY DERIVED FROM D** (the noise symmetry).

**G4-3 access hostile.**

| system | readout | N |
|---|---|---|
| 2 qubits | all effects | 4 |
| 2 qubits | collective S_z only | **3** |
| qutrit | readout {\|0⟩⟨0\|, \|1⟩⟨1\| + \|2⟩⟨2\|} | **2** (same N as a qubit, different d) |
| qubit | unsharp effects (spectrum ⊂ [0.1, 0.9]) | **1** |

→ **CAPACITY A-PRICED.**

**G4-4 noise** (depolarizing qubit; ε = 0.05):

| p | exact N | ε-capacity | Holevo (bits) |
|---|---|---|---|
| 0 | 2 | 2 | 1.0000 |
| 0.01 | **1** | 2 | 0.9546 |
| 0.1 | 1 | 2 | 0.7136 |
| 0.3 | 1 | **1** | 0.3902 |

The three notions disagree. They are not merged.

**G4-5 composition.**
- Classical N_AB = N_A N_B.
- Complex QM: N = d_A d_B and K + 1 multiplicative.
- **Real QM: N = 4 for two rebits, but K = 9 ≠ 8**: capacity is multiplicative while state-space dimension is not.
- Two gbits: the min tensor (local polytope, 16 vertices) and the max tensor (boxworld, 24 vertices, NS check True) both have
  K = 8 and both have **N = 4** (no distinguishable set of size 5 in an exhaustive LP search).
- → capacity alone does not see the composition difference. **COMPOSITION-PRICED for K, not for N, in these cases.**

## Cross-hostiles

| pattern | instance |
|---|---|
| **G2 same, G3 different** | the spatial ℤ² graph gives causal d ≈ 2.9 (finite-speed cones) vs 1.14 (absolute time) |
| **G3 same, G4 different** | causal order constrains no local state space: the same causal set can carry a qubit (N = 2) or a qutrit (N = 3) at each element. Construction-level, but the information accounting is exact: no G3 datum enters N |
| **G4 same, G1 / representation different** | N = 2 for K = 1 … 4 and for the access-restricted qutrit (d = 3) |

**Status: S2-G2/G3/G4 COMPLETE.**
- **G2:** derived from (graph + probe walk); scale-, probe- and site-priced; non-unique as a geometry selector.
- **G3:** derived from causal order + region (estimator-priced); not determined by G2; Lorentz not derived.
- **G4:** derived from (state space + effects + dynamics); A-priced; does not fix the representation.
- **G1:** Σ.
