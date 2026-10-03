# V0-3 TARGET SPEC — P-02 (free-boson collapse; hard-core-boson ↔ SF-1 match; occupation-edge theorem)

**Definitions, assumptions, target statements and acceptance tests ONLY. No proof.**

**Sources** (frozen `scout-0 @ ab2da47`):

| file | blob | what was read |
|---|---|---|
| `playground/SCOUT_0/probes/P02_RESULT.md` | `f766a727…` | heading grep; lines 1 – 104 (verdict, model and checks, zero-mode firewall, minimal-carrier table, theorem with its proof lines) |
| `playground/SCOUT_0/probes/P02_SF1_CONTROL_MAP.md` | `fe73d236…` | all 29 lines (definitions) |
| `SF1_FORMATION_CHARTER_01.md` | `2f05f461…` | grep; lines 96 – 135 (object, support, prescriptions, invariants) |

**Not read:** P-02 scripts and logs; P02B; any other SF-1 file.

**Exposure note.** The orchestrator saw the original's short proof lines (the "Proof of the collapse" paragraph and the
§4 proof sentence). They are **not** reproduced here.

**After this spec is committed, all of the above plus `p02_boson_sectors.*` and the P02B material are SEALED until the
reproduction is committed.**

## 1. SF-1 definitions (frozen from the record)

**Ring and parent (fermionic SF-1).**
- Ring ℤ_L, periodic, L even, spinless, hopping 1.
- H = −Σ_j (c_j†c_{j+1} + h.c.).
- One-particle levels ε(k) = −2cos k, with k = 2πj/L.
- No chemical potential and no sector term.

**Sector selection.**
- Conserved N; the reference state is the **lowest-energy eigenstate within sector N** (a state declaration).
- **Odd-N convention.**

**Readout and object.** ρ_q = Σ_j e^{−iqj} n_j, with q = 2πm/L, and

  S_{L,N}(q, ω) = (1/N) Σ_m |⟨m|ρ_q|0_N⟩|² δ(ω − (E_m − E₀)).

- **Finite-L support:** Σ_{L,N} = {(q, E_m − E₀) : weight > 10⁻¹²}, with the weight summed over each degenerate
  eigenspace.
- **Lower edge:** ω⁻_{L,N}(q) = min{ω : (q, ω) ∈ Σ_{L,N}}.

**Prescription P.**
1. Exact finite-L support.
2. L → ∞ along the declared sequence at fixed q ∈ (0, π], with q_L the nearest grid momentum (ties go to the lower m):
   ω⁻(q) := lim_L ω⁻_{L,N_L}(q_L).
3. q → 0⁺, then z_P := lim log ω⁻(q)/log q.

**Prescription P′.** ω₁(L) := ω⁻_{L,N_L}(2π/L), and z_{P′} := −lim log ω₁(L)/log L.

**Invariants.**
- I-z := z_P.
- I-q := |Q_soft|, where Q_soft = {q* : lim_{q→q*} ω⁻(q) = 0}, taken **modulo q ~ −q ~ q + 2π**, with representatives in
  [0, π] and q* = 0 included via the one-sided limit.
- **Class** := (I-z, I-q).

**Sector families (N rules; L sequences):**

| family | N | L |
|---|---|---|
| D | N = L/2 | L ≡ 2 (mod 4) |
| E | N = 1 | |
| D-¼ | N = L/4 | L ≡ 4 (mod 8) |
| D-¾ | N = 3L/4 | L ≡ 4 (mod 8) |
| E-3 | N = 3 | |
| Ē | N = L − 1 | |
| C-6 | N = nearest odd √L | report-only crossover |

**Recorded SF-1 fermion classes:** D-type (D, D-¼, D-¾) = (1, 2); E-type (E, E-3, Ē) = (2, 1). C-6 reports
k_F(L) = π(N_L − 1)/L and may be **MULTISCALE / PATH-DEPENDENT**.

## 2. Target claims (establish or reject independently)

### TARGET 1 — free-boson collapse

**Parent:** H_B = −Σ_j (b_j†b_{j+1} + h.c.) on the same ring, full bosonic Fock space, with conserved N. Same readout,
object, prescriptions, invariants and families (applied literally as N rules).

**Claims:**
1. The sector ground state is the condensate (b₀†)^N|0⟩/√N! in the **unique** single-particle minimum.
2. ρ_q on it creates only the one-particle transfer out of the condensate mode.
3. The supported excitation energy is ω(q) = ε(q) − ε(0) = 4sin²(q/2).
4. After normalisation, S is **independent of N**.
5. Hence every registered family collapses to the same class (I-z, I-q) = (2, 1), under both P and P′.

**Required:** analytic derivation; finite diagonalisation only as a check.

**Hostile checks:**
- more than one degenerate minimum;
- a flat (non-quadratic) minimum;
- a non-unique sector ground state;
- q connecting degenerate condensate minima.

### TARGET 2 — hard-core boson (HCB) ↔ SF-1 match

**HCB:** n_j ≤ 1, with commuting bosonic operators on distinct sites. **This is an interaction / supplied exclusion rule
— not "free bosons".**

**Claims:**
1. Jordan–Wigner maps HCB to free fermions.
2. On a periodic ring the boundary twist is (−1)^{N−1}, so the odd-N SF-1 convention gives periodic fermions.
3. The density readout is JW-invariant.
4. Hence the support object and the class data equal SF-1's exactly.

**Conclusion to test:** occupancy exclusion suffices; exchange antisymmetry is not required, **in 1D only** (do not
extend to d > 1).

### TARGET 3 — occupation-edge theorem (free, quadratic, number-conserving, single-band, 1D, smooth ε)

**Original statement.** Let O_N be the ground-state occupation set, k_b(L) its edge momenta, and k_b^∞ the limit.
- **(i)** The support is {ε(k + q) − ε(k) : k ∈ O_N, k + q admissible}. "Admissible" means unoccupied under exclusion,
  and anything for free bosons.
- **(ii) Exclusion:** O_N is the N lowest levels, so k_b^∞ = πν for ε = −2cos k (ν = N/L, under the SF-1 momentum
  convention).
  - At a regular edge ω⁻(q) ≈ |ε′(k_b^∞)| q, so z = 1.
  - At an extremal edge, z = 2.
  - I-q counts the soft momenta {0, 2k_b^∞} modulo the quotient.
- **(iii) Free bosons:** O_N = {argmin ε} for every N, so every sector has the band-extremum class: collapse.
- **Original "hence":** "the law class is a function of the single sector invariant k_b^∞ (equivalently the edge
  velocity v_b = ε′(k_b^∞)). Particle statistics act only through ν ↦ k_b^∞."
- **Original wording of (ii):** "z = 1 iff ε′(k_b^∞) ≠ 0; z = 2 iff k_b^∞ is a band extremum."

**Required tasks:**
- Re-derive (i) from ρ_q = Σ_k a†_{k+q} a_k.
- Derive when ε(k_b + q) − ε(k_b) is linear, quadratic or higher order in q.
- **Hostile case:** ε′(k_b) = 0 with the first non-zero derivative of order r > 2 may give z = r, not 2. Address it
  explicitly.
- Audit whether "z = 2 iff extremum" is correct as stated.

## 3. Critical generality audit (required; any needed scope repair is grade-bearing)

Does the law class really factor through the **single scalar** v_b = ε′(k_b^∞)? Test at least:
1. multiple disconnected occupied pockets;
2. multiple inequivalent occupation edges;
3. degenerate single-particle minima;
4. asymmetric bands;
5. flat or higher-order extrema;
6. more than one band;
7. an edge approaching an extremum at an L-dependent rate;
8. whether the soft-point count is determined by v_b alone.

**Possible outcomes:**
- the single-edge theorem holds under extra assumptions;
- a vector / set of edge data (with local dispersion jets and connectivity) is required;
- the exponent and the soft count need separate invariants.

**Further items:**
- **Exclusion map:** derive k_b^∞ = πν and v_b = ε′(πν) for the cosine band, under the exact convention. Do **not**
  assume this for dispersions whose ground-state occupied set is not one contiguous Fermi sea.
- **Free-boson map:** is calling the condensate momentum an "occupation edge" the same mathematical object as a fermionic
  Fermi edge, or only a unifying notation? Audit it.
- **Soft-point count:** for the single-interval exclusion case, test the claimed set {0, 2k_b} modulo the quotient. State
  the assumptions under which it holds. Do not generalise to multiple pockets, condensates or arbitrary bands without
  derivation.
- **C-6:** with N = nearest odd √L, the original claims k_b(L) ~ L^{−1/2}, so v_b(L) → 0 at a rate competing with
  q ~ 1/L, and "the path dependence is the two-scale approach of v_b to 0". Analyse this separately from the fixed-density
  and fixed-N limits.

## 4. Acceptance tests (finite exact; fresh code)

| test | content |
|---|---|
| **A1 (bosons, exact sector blocks)** | L ∈ {6, 10}, N ∈ {1, 3, 5}. **Original:** the sector ground state is unique, with gap 1.000 (L = 6) and 0.382 (L = 10); E₀ = Nε(0); the grouped ρ_q support is the single line ω = 4sin²(q/2) for every q |
| **A2 (HCB vs SF-1)** | for **all odd-N sectors** at L ∈ {6, 10, 12}, the HCB grouped ρ_q support equals the free-fermion (SF-1) support exactly. Even N gives the twisted ring and is not used (the original notes L = 12, N = 6) |
| **A3 (classes)** | the SF-1 fermion classes D, D-¼, D-¾ = (1, 2) and E, E-3, Ē = (2, 1), with v_b = 2, √2, √2 (D, D-¼, D-¾) and v_b = 0 (E, E-3, Ē), and z_P = z_{P′} for each. Free bosons give (2, 1) for every family, including C-6 |
| **A4 (optional diagnostic, "Fermi pattern in bosons")** | the free-boson eigenstate with one boson in each of the N lowest levels is an exact eigenstate but not the ground state. **Original:** E − E_gs = 372.8 at L = 1026, and the support contains ω < 0, with min ω = −1.99, −0.76, −0.20 at q = π/2, π/8, π/32. The sector is not stated in the extracted lines; D-type N = L/2 is presumed and must be checked |

## 5. Firewalls

- P-02b (interacting / Luttinger physics) is **not** part of V0-3 and may not be used to rescue or modify the free
  theorem.
- Classical exactness only: no μ, no regulator.
