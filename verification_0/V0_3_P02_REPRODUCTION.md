# V0-3 — P-02 blind reproduction (free-boson collapse; HCB ↔ SF-1; occupation-edge theorem)

**Firewall.** The only project file read was `verification_0/specs/V0_3_P02_TARGET.md`. Everything below was derived
from that spec, standard results (Jordan–Wigner, Fock-space diagonalisation of quadratic Hamiltonians, free-particle form
factors), and fresh code in `verification_0/code/v0_3_*.py`. Logs: `v0_3_*.log`, next to the scripts.

**Conventions (fixed once).** Sites j = 0..L−1. H = Σ_j Σ_r h_r a†_j a_{j+r} with h_{−r} = h_r*.
c_k = L^{−1/2} Σ_j e^{−ikj} c_j, so H = Σ_k ε(k) n_k with ε(k) = Σ_r h_r e^{ikr}. The SF-1 band is h_{±1} = −1,
ε = −2cos k. The readout is

  ρ_q = Σ_j e^{−iqj} n_j = Σ_k a†_{k−q} a_k   (for q > 0 a particle moves from k to k − q).

The spec writes Σ_k a†_{k+q} a_k. That is the same operator in the opposite Fourier convention, which also sends ε(k) to
ε(−k). The two are equivalent for reflection-symmetric bands. For asymmetric bands the choice decides which edge P sees
(§3, CE-4), so it has to be stated with ε.

---

## TARGET 1 — free-boson collapse

**Derivation.** H_B = Σ_k ε(k) n_k. The Fock states |{n_k}⟩ are an eigenbasis with E = Σ_k n_k ε(k) ≥ N·min ε. Equality
holds iff every boson sits in argmin ε. For ε = −2cos k, argmin = {0} on every L, so:

1. The sector-N ground state is unique: (b₀†)^N|0⟩/√N!, with E₀ = Nε(0) = −2N.
2. The gap is the cheapest single transfer, ε(±2π/L) − ε(0) = 4 sin²(π/L). It does not depend on N: 1.000 at L = 6 and
   0.381966 at L = 10.
3. ρ_q acting on |N₀⟩ keeps only the k = 0 term: b†_{−q} b₀ |N₀⟩ = √N |(N−1)₀, 1_{−q}⟩ for q ≠ 0. This is one particle
   transferred out of the condensate mode (claim 2).
4. The transfer energy is ε(−q) − ε(0) = 4 sin²(q/2) (claim 3), with weight N/N = 1.
5. So S(q, ω) = δ(ω − 4 sin²(q/2)) for every q ≠ 0, and this does not depend on N (claim 4).
   **Qualification:** at q = 0, ρ₀ = N̂ gives S(0, ω) = N δ(ω), which does depend on N. P never uses q = 0, so the
   classes are unaffected.

**Classes.**
- **Under P:** ω⁻(q) = 4 sin²(q/2) for every family, whatever the N rule. This gives z_P = 2 and Q_soft = {0}, so
  I-q = 1.
- **Under P′:** ω₁ = 4 sin²(π/L), so z_P′ = 2.
- Every family, C-6 included, collapses to **(2, 1)** under both prescriptions (claim 5).
- **Status:** established for the cosine band.

**Hostile checks** (code: `v0_3_counterexamples.py`, CE-1 and CE-2):

| hostile case | result |
|---|---|
| more than one degenerate minimum, ε = 2cos k + cos 2k (minima ±2π/3) | The GS manifold of sector N has dimension N + 1 (ED: 3 at L = 6, N = 2; 4 at L = 6, N = 3 and L = 12, N = 3). The "sector GS" is undefined unless a state is declared. |
| non-unique sector GS | Same model. The support depends on the declared state. At q = 2π/3 (≡ 2k₀ under the quotient) there is ω = 0 weight **0** for "all bosons in +k₀", **1** for "all in −k₀" or n₊ = n₋ = 1, and **0.5** for (\|2,0⟩ − \|0,2⟩)/√2. The class is therefore (2, 1) or (2, 2) depending on the declaration. Even for fermions with N = 1 on this band the GS is degenerate. |
| q connecting the degenerate minima | Soft at q* = k_i − k_j exactly when the declared state has a non-zero amplitude to move a boson from minimum k_i to minimum k_j. Weights do not simply add for superpositions; they can interfere. |
| flat (non-quadratic) minimum, ε = −4cos k + cos 2k, so ε − ε(0) = 8 sin⁴(k/2) | GS unique and support a single line 8 sin⁴(q/2) (ED at L = 6, 10). z_P = z_P′ = 4, so the free-boson class is **(4, 1), not (2, 1)**. |

**Corrected statement.**
- Free bosons with a **unique** minimum k* of order r* (the first non-zero derivative) have class (r*, 1) for every N
  rule, under both P and P′.
- "Collapse to (2, 1)" therefore needs a unique, quadratic minimum, which the cosine band has.
- With degenerate minima, the object is undefined until a ground state is declared.

---

## TARGET 2 — HCB ↔ SF-1 via Jordan–Wigner (1D)

**Mapping.** b_j = exp(iπ Σ_{l<j} n_l) c_j. This maps n_j ≤ 1 hard-core bosons onto fermions, unitarily on the 2^L
space and within each N sector.

**Bulk bonds.** b†_j b_{j+1} = c†_j c_{j+1}: the strings cancel.

**Boundary bond.** b†_{L−1} b₀ = c†_{L−1} (−1)^{Σ_{l<L−1} n_l} c₀. Once c₀ has acted and site L−1 is empty, the string
counts the other N − 1 particles. So b†_{L−1} b₀ = (−1)^{N−1} c†_{L−1} c₀, which gives:
- the boundary twist (−1)^{N−1} (claim 2);
- **periodic fermions for odd N, antiperiodic for even N.**

**Density readout.** n_j = b†_j b_j = c†_j c_j, so ρ_q is JW-invariant (claim 3).

**Consequence.** For odd N the HCB and SF-1 sector Hamiltonians are unitarily equivalent, and ρ_q is fixed under the map.
Spectra and every matrix element ⟨m|ρ_q|0⟩ coincide, so S is the same measure. Supports are therefore equal at every
finite L, and every P or P′ limit along any odd-N sequence is equal too. All registered families have odd N
(D: L ≡ 2 mod 4; D-¼ and D-¾: L ≡ 4 mod 8; E-3; Ē; C-6 by definition). So the class data are identical (claim 4).

**Scope repairs** (code: `v0_3_A2_hcb_vs_fermions.py`, CE-7):

| condition | evidence |
|---|---|
| **Odd N only** | HCB with even N equals antiperiodic fermions. Example L = 12, N = 6: E₀(HCB) = −7.7274066 equals the antiperiodic value, while periodic SF-1 has −7.4641016 and is 2-fold degenerate. Same pattern at (6,2), (10,4), (12,4). |
| **Nearest-neighbour hopping only** | With an extra next-nearest hop, the string (−1)^{n_{j+1}} survives and HCB becomes an interacting fermion model. ED at L = 10, N = 3 with t₂ = 0.5: E₀ = −7.389 (HCB) vs −6.854 (free F), and the supports differ. |
| **Periodic ring** | With open boundaries there is no twist at all. |
| **d = 1 only** | 4×4 torus, N = 5: E₀(HCB) = −15.3004 vs −12 for free fermions. The HCB support at q = (π/2, 0) has a dominant line at 2.888 plus a continuum; free fermions give a single point at 2.000. The JW strings do not cancel in d > 1. |

"Exclusion suffices, exchange antisymmetry not required" holds exactly in 1D for nearest-neighbour hopping on a periodic
ring with the odd-N convention. Under those conditions the hard-core constraint plus the JW string *is* the
antisymmetry.

---

## TARGET 3 — occupation-edge theorem

### (i) Support from ρ_q

**Fock-state reference.** A unique sector GS of a quadratic H is a Fock state |n⟩: the Fock states diagonalise H, and a
non-degenerate eigenvector must be one of them.

**Action of ρ_q.** For q ≠ 0, a†_{k−q} a_k|n⟩ is a Fock state with energy shift ε(k−q) − ε(k) and amplitude:
- fermions: ±1 if n_k = 1 and n_{k−q} = 0, else 0;
- bosons: √(n_k (n_{k−q} + 1)).

**No interference.** Different k give orthogonal final Fock states, and grouping by energy adds non-negative weights. So
the support is exactly

  {ε(k−q) − ε(k) : n_k > 0, k−q admissible}

where "admissible" means unoccupied for fermions and anything for bosons. The spec writes the same set as
ε(k+q) − ε(k) in its convention.

**Status:** correct, under two conditions:
- **the reference must be a single Fock state.** It fails for superpositions in a degenerate manifold, where
  interference occurs (CE-2: weight 0.5).
- **the Fourier convention must be stated together with ε.**

**Numerical checks.** Verified against real-space ED for:
- the cosine band;
- the cubic band (L = 10, 14);
- the asymmetric band and its mirror (L = 8, 12). The ED matches the k → k−q rule and **not** the k → k+q rule in the
  stated convention.
- the three-pocket band (L = 12);
- the ladder (L = 6);
- the boson "Fermi pattern" with ω < 0 (L = 6, 10).

### (ii) Exclusion, cosine band

**Occupied set.** With L even and N odd, the levels pair as ±k and the shell is closed. O_N = {2πj/L : |j| ≤ (N−1)/2} is
unique.

**Edge momentum and velocity.**
- k_b(L) = π(N−1)/L → **k_b^∞ = πν**.
- **v_b = ε′(πν) = 2 sin(πν)**: 2 for D, √2 for D-¼ and D-¾, 0 for E, E-3 and Ē.
- Using the midpoint πN/L instead of π(N−1)/L for k_b(L) makes no difference at fixed ν.

**Closed-form lower edge (derived; checked against brute force up to L ≈ 2×10⁶).**
- For ν ≤ ½ and 0 < q ≤ 2k_F, the admissible k lie in [−k_F, −k_F + q). There
  ε(k−q) − ε(k) = 4 sin(q/2) sin(q/2 − k), and concavity puts the minimum at the endpoint.
- For q ≥ 2k_F all of O is admissible.
- Particle–hole symmetry (c_j → (−1)^j c_j†) maps ν to 1 − ν. This makes D-¼ and D-¾ **identical even at finite L**,
  which the log confirms.

Altogether:

  **ω⁻(q) = 4 sin(q/2) · |sin(π·min(ν, 1−ν) − q/2)|,   ν ∈ [0, 1].**

**P′, exactly:** ω₁(L) = 4 sin(π/L) sin(πN/L).

**Classes:**
- For ν ∈ (0, 1): ω⁻ ≈ v_b q, so z_P = z_P′ = 1. Q_soft = {0, 2π·min(ν, 1−ν)}, so I-q = 2. The second soft point is π
  for D and π/2 for D-¼ and D-¾.
- For ν ∈ {0, 1}: ω⁻ = 4 sin²(q/2), so z = 2 and I-q = 1.

**Leading behaviour of ε(k_b + q) − ε(k_b).** Let r be the first order with ε^{(r)}(k_b) ≠ 0. Then
ε(k_b + q) − ε(k_b) = ε^{(r)}(k_b) q^r / r! + O(q^{r+1}). It is linear iff ε′ ≠ 0, quadratic iff ε′ = 0 and ε″ ≠ 0, and
of order r otherwise.

**The lower edge is not this single difference.** It is a minimum over the window of pairs that straddle the edge.

**Fermi edge (μ crossing, r odd).** With ε − μ ≈ c (k − k_b)^r:

  ω⁻(q) = |c| q^r · min_{t∈[0,1]} [t^r + (1−t)^r] = |c| 2^{1−r} q^r,

so z = r. Checks:
- r = 1 gives |ε′| q.
- r = 3 gives q³/12 for −cos³(k)/3 and q³/16 for the asymmetric band. Both are verified to 4 digits.

**Degenerate sea (ν → 0 or 1, fixed N), or a condensate.** At an extremum of even order r,
ω⁻(q) ≈ |ε^{(r)}| q^r / r!, so z = r.

**Hostile r > 2.**
- A quartic minimum gives z = 4 (CE-1, z_P′ = 3.9999).
- A stationary-inflection edge gives z = 3 (CE-3, z_P = 3.000, z_P′ = 3.0000).
- So z = r, not 2.

**Audit of "z = 1 iff ε′(k_b^∞) ≠ 0; z = 2 iff k_b^∞ is a band extremum."**
- The first half is correct for a single edge type.
- The second half is **false as stated**:
  - An extremum of order r ≥ 4 gives z = r.
  - A non-extremal stationary edge (an odd-order inflection) gives z = r ≥ 3, not 2.
- **Correct statement:** **z = r(k_b)**, the order of the first non-vanishing derivative at the edge seen by P.
  z = 2 iff the edge is a non-degenerate (ε″ ≠ 0) extremum. For a fermionic sea that happens only in the degenerate
  ν → 0 or 1 limits, or at tuned Lifshitz points.
- For the cosine band, which has only quadratic extrema and ε′ ≠ 0 on (0, π), the original wording is true.

### (iii) Free bosons

O_N = argmin ε, independent of N. "Collapse" holds iff the minimum is unique. The class is then (r*, 1), which is (2, 1)
only when the minimum is quadratic (see TARGET 1).

### "Hence": does the class factor through v_b?

**No, in general.** It does for the cosine band, and in fact only through the Boolean "v_b ≠ 0". The counterexamples are
in §3.

---

## 3. Critical generality audit

| # | test | model (code CE-#) | result |
|---|---|---|---|
| 1 | disconnected pockets | ε = −2cos 3k − ½cos k, μ = −1 (CE-5); three pockets, 6 Fermi points, ν = 0.3299. Odd-N closed shells are unique, e.g. L = 6000, N = 1981 | z_P = 1. Q_soft = all pairwise Fermi-point differences under the quotient: **I-q = 10** (predicted set = numerical zeros). There is no single k_b or v_b. |
| 2 | inequivalent edges | asymmetric ε = ½sin k − ¼sin 2k, ν = ½ (CE-4): edges at 0 (r = 3, v = 0) and π (r = 1, \|v\| = 1) | **z_P = 1 for this band and z_P = 3 for its mirror ε(−k)**, with the same edge data and both I-q = 2. P (q ∈ (0, π]) probes only edges whose occupied side lies at larger k (in our convention). The quotient q ~ −q hides the other orientation. |
| 3 | degenerate minima | 2cos k + cos 2k (CE-2) | Bosons: GS degenerate, so the class depends on the declared state ((2,1) vs (2,2)). Fermions: the E family has a degenerate GS. |
| 4 | asymmetric bands | CE-4 | O = (−π, 0), not [−πν, πν], so **k_b^∞ = πν fails**. Chirality decides z_P. The Fourier convention must be stated. |
| 5 | flat / higher-order extrema | 8 sin⁴(k/2) (CE-1); −cos³k/3 (CE-3) | Bosons (4, 1); fermion E (4, ·). For D at ν = ½ in the cubic band, **v_b = 0 but the class is (3, 2)**. Cosine-band E, E-3 and Ē also have v_b = 0 but class (2, 1). |
| 6 | two bands | ladder, t⊥ = 1, μ = ½ (CE-6); ED-verified intraband-only support for the leg-summed ρ_q | Two Fermi momenta k₋ = 2.4189 and k₊ = 1.3181, with velocities 1.323 and 1.936. Q_soft = {0, 1.4455, 2.6362}, so **I-q = 3**, z = 1. |
| 7 | L-dependent edge | C-6 (see below) | P and P′ disagree: z_P = 2, z_P′ = 3/2. |
| 8 | soft count from v_b? | all of the above | **No.** v_b = 0 occurs with I-q ∈ {1, 2}, and v_b ≠ 0 occurs with I-q ∈ {2, 3, 10}. For a single interval the soft set is {0, 2πν}, which depends on ν (the interval length), not on v_b. |

### Exclusion map

k_b^∞ = πν and v_b = 2 sin(πν) are derived above for the cosine band with odd N and periodic boundaries. They hold for
a symmetric band whose ground-state occupied set is a single interval centred at 0. They fail for asymmetric bands
(CE-4) and for multi-pocket seas (CE-5, CE-6).

### Free-boson "occupation edge": a notational unification, not the same object

**Fermionic k_b.** A Fermi edge is a boundary point ∂O_N of an occupied set of measure ν. The occupation n(k) jumps
1 → 0 there. It moves with ν, it has an orientation, and low-energy weight comes from pairs on either side of the edge.

**Bosonic "edge".** The condensate momentum is argmin ε, a band property that does not depend on N or ν. It is the
support of an occupation that is a δ of weight N. There is no unoccupied side and no orientation: a transfer from a single
mode replaces the window minimisation.

**Where the two agree, and where they do not.**
- The two agree, with identical normalised supports, only with the N = 1 sector, i.e. the ν → 0 degenerate-sea limit
  of the fermion problem.
- Their exponents agree with each other only in that limit: z = order of the minimum, which is even.
- At fixed ν ∈ (0, 1), "statistics act only through ν ↦ k_b^∞" is false for bosons, because ν plays no role at all.
- What statistics really change is (a) the occupation rule (sublevel set vs argmin) and (b) admissibility.

### Soft-point count, single interval

Take a single occupied interval of length 2πν, with ε = μ only at its two endpoints (no tuned Lifshitz touching).
- Q_soft = {0, ±2πν}, so I-q = 2 iff ν ∉ {0, 1}, and I-q = 1 otherwise. This holds even for asymmetric bands, where
  2k_b is replaced by the interval length 2πν.
- With several pockets or bands: Q_soft = {k_e − k_f : e, f Fermi points} under the quotient (CE-5, CE-6).
- For bosons, Q_soft = {differences of the occupied minima}. For superpositions it further requires non-zero
  transfer amplitudes (CE-2).

### C-6 two-scale analysis (N = nearest odd √L)

**Exact finite-L lower edge** (verified against brute force with 0 mismatches at L = 1000, 4096 and 10000):
- ω⁻_L(m) = 4 sin(πm/L) sin(π(N+1−m)/L) for m ≤ N;
- ω⁻_L(m) = 4 sin(πm/L) sin(π(m−N+1)/L) for m > N.

**Scaling form.** For small q:

  **ω⁻ ≈ q·|q − 2k_F(L)|**,  with k_F(L) ≈ π L^{−1/2}.

There are two scales: q and k_F(L), equivalently v_b(L) = 2 sin k_F ≈ 2π L^{−1/2}. On top of them there is a soft point
at 2k_F(L) that moves into q = 0.

| path | q | ω behaviour | exponent | class / remark |
|---|---|---|---|---|
| **P** (fixed q) | fixed | → 4 sin²(q/2) | z_P = 2 | (2, 1) |
| **P′** | 2π/L | ω₁ = 4 sin(π/L) sin(πN/L) ≈ 4π² L^{−3/2} | z_P′ = 3/2 (fit 1.4995) | — |
| general path | c L^{−a} | — | z(a) = 2 for a ≤ ½; z(a) = 1 + 1/(2a) for a > ½ | verified: a = 0.6 → 1.833, 0.75 → 1.666, 0.9 → 1.557 |
| along the moving soft line | q ≈ 2k_F(L) (a = ½, c = 1) | — | 3 | — |
| N ~ L^α | 2π/L | — | z_P′ = 2 − α | — |

**Comment on the original.** The original says the path dependence is the two-scale approach of v_b to 0. That is
correct but incomplete. The soft count is also limit-order dependent: at finite L there are two soft momenta,
{0, 2k_F(L)}, which merge under P. C-6 is genuinely **MULTISCALE / PATH-DEPENDENT** and has no single (I-z, I-q).

### Corrected theorem

**Setting.** A free, quadratic, number-conserving, translation-invariant 1D ring. ε is smooth. The sector GS is unique
along the sequence, so it is a Fock state. Fixed ν, or fixed N.

**Fermions at fixed ν ∈ (0, 1).** Let F be the set of Fermi points (ε = μ, with no tuned touchings). Each e ∈ F has a
local order r_e (odd) and an orientation σ_e.
- **I-z = z_P = z_P′ = max{r_e : e ∈ F with σ_e = the orientation probed by q > 0}.**
- For a reflection-symmetric band this is max_e r_e.
- ω⁻(q) ≈ min_e |c_e| 2^{1−r_e} q^{r_e}.
- **Q_soft = {k_e − k_f : e, f ∈ F} under the quotient q ~ −q ~ q + 2π**, and I-q = |Q_soft|.

**Degenerate seas** (fixed N, or ν → 0 or 1) **and bosons** with a unique extremum of order r*: class (r*, 1).

**Degenerate extrema:** the class is undefined until a reference state is declared.

**Invariant needed.** The class depends on the **set of edge data** {(k_e, σ_e, r_e)}: positions modulo translation, for
the soft count; orientation and local jet order, for the exponent. The exponent and the soft count are separate
invariants. **They do not factor through the scalar v_b.**

**Single-scalar corollary (the original's "hence", restricted).** Assume:
- ε is even;
- ε′ ≠ 0 on (0, π);
- ε″(0) ≠ 0 and ε″(π) ≠ 0 (the cosine band satisfies all three);
- fermions with fixed ν;
- unique GS.

Then (I-z, I-q) = (1, 2) iff v_b ≠ 0 iff ν ∈ (0, 1), and (2, 1) iff v_b = 0 iff ν ∈ {0, 1}. Only the zero / non-zero
status of v_b matters.

**L-dependent edges (C-6)** are outside the theorem: P and P′ differ.

---

## 4. Acceptance tests

| test | result | key numbers | verdict |
|---|---|---|---|
| **A1** bosons, L ∈ {6, 10}, N ∈ {1, 3, 5} | GS unique (degeneracy 1) in all 6 sectors; overlap with the condensate = 1 to 1e-15 | E₀ = −2N exactly. Gap 1.000000000000 (L = 6) and 0.381966011250 (L = 10), = 4 sin²(π/L) and independent of N. Every m ≠ 0: one line ω = 4 sin²(q/2), weight 1. m = 0: ω = 0 with weight N. | **PASS** (qualification: S(0, ·) depends on N) |
| **A2** HCB vs SF-1, all odd N at L ∈ {6, 10, 12} (14 sectors) | HCB ED = fermion ED (canonical signs) = free-fermion form factors | Full spectra equal; grouped supports *and weights* equal; GS unique in every sector. Even N fails: L = 12, N = 6 has E₀(HCB) = −7.72741 vs −7.46410 (periodic SF-1, 2-fold degenerate). | **PASS** (1D, NN, odd N only) |
| **A3** classes | Fermions: D (1, 2), v_b = 2; D-¼ and D-¾ (1, 2), v_b = √2 (identical by particle–hole); E, E-3, Ē (2, 1), v_b = 0 | z_P = z_P′ for all six: 1.0/1.0 and 2.0/2.0. Bosons: (2, 1) for all seven families including C-6, z_P′ = 2.00000. **C-6 fermions:** z_P = 2, I-q = 1 but z_P′ = 1.4995 → 3/2. | **PASS** for the stated items. C-6 confirmed MULTISCALE (P ≠ P′). |
| **A4** "Fermi pattern" in bosons | Exact eigenstate, not the GS (ED, L = 6 and 10). At L = 1026 **only N = 513 = L/2 (family D)** gives E − E_gs within 0.05 of 372.8. | E − E_gs = 2N − 2/sin(π/L) = **372.8271**. min ω = **−1.993876** (q = π/2 → m = 256 by the tie rule), **−0.763484** (π/8 → m = 64), **−0.195623** (π/32 → m = 16). These round to the quoted −1.99, −0.76, −0.20. The thermodynamic limit is −2 sin q; the finite-L formula is −4 sin(πm/L) cos(π(m+1)/L). | **PASS**: the sector is D (N = L/2) |

## 5. Corrections to the spec and the original claims

1. **T1, claim 4** holds only for q ≠ 0, because S(0, ω) = N δ(ω).
2. **T1, claim 5 and (iii)** need a unique, quadratic minimum. In general free bosons have class (r*, 1), and the
   object is undefined for degenerate minima.
3. **(i)** needs a Fock-state reference and a stated Fourier convention. The spec's a†_{k+q} a_k and the real-space
   e^{−iqj} readout correspond to opposite conventions for ε(k), which matters for asymmetric bands.
4. **"z = 2 iff extremum"** is false. The correct rule is z = r, the order of the first non-zero derivative at the edge.
5. **k_b^∞ = πν** holds only for a single symmetric interval.
6. **{0, 2k_b}** holds only for a single interval without tuned touchings; it is really {0, 2πν}.
7. **The class does not factor through v_b**, and statistics do not act only through ν ↦ k_b^∞. The corrected invariant
   is the set of edge data (above).
8. **T2** requires nearest-neighbour hopping, odd N, a periodic ring and d = 1.
9. **C-6:** z_P = 2 but z_P′ = 3/2. The soft count is also path-dependent; the moving 2k_F(L) is the second scale.
10. **A4:** the sector is confirmed as D (N = 513). The quoted numbers are reproduced to rounding.

No spec number was found to be wrong.
