# V0-4 — EDA-01 SECOND READER (blind, classification-first)

**Reader:** independent second reader. It is context-isolated from the original audit: it never opened
`EDGE_DATA_AUDIT_01.md`, the `eda_witnesses` files, or any other deny-listed file (the full read list is in §8).
**Spec:** `verification_0/specs/V0_4_EDA01_TARGET.md`, applied as written.
**Fresh code:** `verification_0/code/v0_4_sr_w{1..4}_*.py`. Each script has a matching `.log`.
**Status:** nothing committed or pushed. No record file was modified.

---

## 0. Conventions

### How the rubric is applied

**FIX** means the earned predicate *as banked* pins the value of an edge datum. The test is the spec's decisive
test:

- Any two models that satisfy the earned predicate must agree on the datum.
- Equivalently, the value must not change when a supplied primitive (Table 3, S-1 … S-13) is varied inside the
  predicate's scope.

**Rejected literal reading.** On a literal reading, "with no supplied primitive varied" lets any member-level
computation count as FIX. A theorem about one declared parent computes numbers for that parent once every input
is fixed. On that reading E-12, E-13, E-14 and E-19 would be "FIX (member-scoped)". I reject this reading because
it makes FIX the same as "computed from supplied data", which is the definition of READ. It is recorded as an
interpretation flag in §6.

### Edge-data classes used

The five classes in spec §1:

- (L) edge location / spectral bottom;
- (X) edge exponent / local DOS exponent;
- (V) edge velocity at an occupation edge;
- (S) visible spectral sign;
- (P) number and location of soft points.

**Added explicitly, as spec §1 allows:**

- (J) **edge jet / first non-zero derivative order r_e.** This is what V0-3's corrected quotient actually uses.
- (M) **low or inverse spectral moments.** Examples are K₁₁ = ∫λ dμ₁₁, a resistance distance (a Laplacian
  pseudo-inverse element) and r = (K⁻¹)₁₁.

  These are not edge data in the strict sense. They are listed because three rows can only be labelled READ once
  (M) is admitted. I give both labels for those rows.

### Upstream primitives

Upstream primitives use Table 3 labels:

| label | primitive |
|---|---|
| S-1 | substrate K, net and pin |
| S-2 | generator and orientation |
| S-3 | lift and statistics |
| S-4 | ħ |
| S-6 | gravitational quadruple |
| S-7 | sector / boundary / conservative realization |
| S-8 | environment / noise |
| S-9 | L0-1c drift |
| S-10 | access / readout |
| S-11 | preparation / partition / initial class |
| S-12 | coarse-graining |
| S-13 | cosmological transport inputs |

---

## 1. Row table E-1 … E-22

| ID | Label | Datum involved | Upstream supplied primitive | Justification and FIX test |
|---|---|---|---|---|
| **E-1** gap → P_memory; passivity → P_positivity | **READ / CON** | (L) spectral bottom = gap = pin; (S) sign of bath spectrum (passivity ⇔ λ_min ≥ 0); (X) ungapped edge exponent (τ^{-3/2} envelope) | S-1 (C1-a substrate including pin); S-2 (relaxational G-D kernel e^{−Kτ}) | The record's own identity is k_pin(τ) = e^{−pin·τ}k₀(τ), so the memory rate *is* the supplied spectral bottom (READ). The necessity direction restricts the datum to gap > 0 and spectrum ≥ 0, and leaves the value free (CON). **Two models:** pin = 0.3 and pin = 0.6 are both EXPONENTIAL-grade. Their decay rates are 0.30 and 0.60, and k_pin/(e^{−pin τ}k₀) = 1 exactly (W2a). The record has its own pair too: the L0-1b members keep P_memory with λ_min ranging from 0.304 to 0.383. |
| **E-2** locality → P_geometry | **IND** (strict classes) / **READ** if (M) is admitted | Strict classes: none. Extended: the resistance metric is a Laplacian pseudo-inverse functional, i.e. inverse/low spectral data (M). | S-1 (weights w_ij ∝ \|i−j\|^{−α}); S-10 (the geometry criterion and access) | P_geometry is a geometry criterion, not an effective law's edge datum. **Two models:** in the record, the α = 4 member and the anchor both carry the line (Q = 18.0 vs 23.0) at different λ_min, so nothing spectral is pinned. |
| **E-3** linearity → P_exact-reduction | **IND** | none | — | The result is definitional (the ruling says so). It makes the "linear spectral machinery" *available*. It gives no spectral value, so the FIX test cannot even be posed. |
| **E-4** cycle affinity → ¬P^resp_CM | **CON** (dual **READ**) | (S) visible spectral sign: CM ⇔ a positive measure on [0, ∞) (Bernstein); complex eigen-pairs | S-2 (the generator's affinity / orientation) on the S-1 ring | **Two models:** the banked predicate is an implication. A zero-affinity ring satisfies it vacuously and is CM (L-1). A non-zero-affinity ring satisfies it and is ¬CM. So the sign bit is a function of the supplied affinity (READ). Among ¬CM models the complex-pair locations move with the rates. Analytic witness: a directed ring with forward/backward rates a, b has λ_k = (a+b)(1−cos θ_k) + i(a−b) sin θ_k, so Im λ_k ∝ (a−b). Only "not CM" is fixed, and only once affinity ≠ 0 is supplied (CON). |
| **E-5** seven floor non-implications | **IND** (dual **READ** on two entries) | Five entries: none. "spectral stability ⇏ accretivity ⇏ monotone response": (S), shown *insufficient*. "emergent dissipation ⇏ strict ordering": (X) band-edge t^{-3} ripples of the O-6 chain. | S-1, S-2, S-11 | A non-implication is witnessed by a two-model pair. It is an anti-FIX by construction: it says the datum does *not* determine the property. The O-6 witness reads the supplied chain's edge exponent (READ for that entry). |
| **E-6** no single deep organizing pair | **IND** | none | — | This is a structural/ontological terminal. No spectral value is stated. |
| **E-7** TRIVIAL/IDENTITY observability | **IND** | none. det DΦ = ∏w_k^{N−k} is in coupling weights, not spectral data. | S-10 (the h = x₁ readout); S-1 | The observability rank is constant for *any* weights. It has no spectral content. |
| **E-8** P-1/X2 intrinsic partition selection | **IND** | none. "gap 0.678" is a selection-score margin, not a spectral gap. | S-1 (the supplied impurity structure of the X2 world) | This is a selection terminal. Its content is "which subsystem", not a spectral value. |
| **E-9** spine intact; P-2 cone 𝔠_Gauss = {J ≥ 0, ν ≥ ħJ/2} | **CON** | (S) the sign of the influence spectral density J(ω) ≥ 0, plus the noise floor | S-4 (ħ, "located, not generated"); S-3 (quantum class) | The cone restricts sign and floor. The values stay free: the record's own PC-D says λJ is admissible for every λ > 0. **Two models (W4):** J₁ = ωe^{−ω} and J₂ = 7ω³e^{−ω/2} are both in the cone, with spectral-bottom exponents 1 vs 3. |
| **E-10** D-1 ontology demotion at fixed influence data | **IND / SUP** | The (K, N) influence data are *held fixed as input* | S-1, S-8 (via (K, N)) | The result quotients microscopic realizations at fixed spectral input. It transports the input and never changes it. |
| **E-11** exact FRW kernel transport computable | **SUP** (dual **READ**) | K(t, t′) spectral content: KMS/T, the Matsubara pole ladder 2πnT, H(t) | S-13 (FRW state, KMS/T, stationary reduction, dS→FRW map) | "Computable from already-admitted inputs" is the definition of transported input. **Two models:** different T or H(t) histories give different kernels, and the same "computable" predicate holds for both. |
| **E-12** GR2-L6 counting asymptote | **READ** | (X) the low-frequency mode-counting exponent: 2d = 6 from d = 3 and acoustic (linear) edge dispersion; +0.2142 window curvature; O(1/s) open-boundary term | S-1 (lattice class (b), d, open-chain BC, dispersion) plus the declared window | **Two models:** the same mechanism on a d = 2 cube gives 2d = 4. A different window gives a different curvature offset; the record's own halved window gives 6.0454. Note that 6.2142 is a slope on [0.9, 1.8], *away* from the band bottom, so only the "2d" part is an edge exponent (flag in §6). |
| **E-13** SF-1 FORMATION-OF-LAW-CLASS | **READ** (dual **CON**) | the full edge-data set (L, V, J, P); z = the dispersion order at the relevant edges | S-1 (one-particle dispersion ε(k) of the ring hopping); S-3/S-7 (CAR/exclusion statistics; conservative realization); S-7 (the conserved sector ν or N) | See §4. **Two models (W1):** parents ε_A = −2cos k, ε_B = −2cos k + ½cos 2k and ε_C = −2cos k + 2cos 2k all satisfy the predicate. Their dilute classes are (2, 1), (4, 1) and (2, 2), and their dense edge sets differ (two vs four Fermi points). |
| **E-14** S5-1 MARKOV-LIMIT-OTHER-CLASS | **READ** | (L) band [0.3, 4.3] and absence of bound states; (X) square-root edges → t^{-3/2} tails; the interior density at ω_s → κ = 1/(2√2.3) | S-1 (K, pin); S-2 (the declared conservative parent q̈ = −Kq); the admitted L-vH deformation (a declared limit) | **Two models (W2):** varying pin gives κ = 1/(2√(2+pin)) (0.3297 → 0.3101 at pin = 0.6), and M-1 YES / M-2 NO still hold. A Gegenbauer-graded nearest-neighbour chain with the same band edges gives tail exponent 2.48 instead of 1.49, while "native → non-Markovian dissipation only" still holds. If ω_s² is moved outside the band (a supplied choice of K₁₁), a bound state appears and the Markov damping is lost. So the result also *uses* the location of the edge relative to ω_s (CON aspect). |
| **E-15** S2-1 FULL-DISCRIMINATOR | **IND** (strict) / **READ** if (M) is admitted | Strict classes: none (Δc₂ = −24βT₁a is local noise × drift curvature). Extended: Δc₃ contains K₁₁, the first moment of the site-1 spectral measure (M). | S-9 (β, a); S-8 (T₁); S-1 (K₁₁) | The discriminator's leading term has no spectral datum. The O(t³) coefficient reads a supplied moment. |
| **E-16** S2-HB UNRESTRICTED-REALIZATION-ONLY | **IND** | none | S-8 | This is a realization/representation terminal (path-space skew product). |
| **E-17** S3-0 FORMULABLE-ONLY-WITH-CHANGE | **NF / IND** | The crossed cell is not formulable at fixed invariants. In the report-only R₂ note, C = T·k makes correlation CM ⇔ response CM, i.e. it reads the response's (S). | S-1, S-2, S-8 | The question "does affinity fix the correlation spectral sign at fixed invariants" cannot be posed in the declared theory (NF). The R₂ note is an identity that reads a supplied response (READ, report-only, non-adjudicative). |
| **E-18** REV-0 NOT SUPPORTED | **IND** | none | — | A graph-governance terminal. |
| **E-19** S6-0/S6-1 NET-ARROW-CONFIRMED; LS-1/2/3 | **READ** | (L) purely a.c. spectrum, no bound state, band edges; (M) r = (K_∞⁻¹)₁₁ = (c − √(c²−4))/2, the resolvent at 0, a function of both band edges; (X) band-edge t^{-6} entropy tails | S-1 (K, pin); S-11 (initial product state, T_s, T_b, partition); S-12 (declared observables / J bookkeeping) | **Two models (W2):** at pin = 0.3 the code reproduces the record exactly: r = 0.582109; X_J = 1.16943 / 9.16943 / 0.33057 / 0.73057; κ check. At pin = 0.1, 0.6 or 1.0, NET-ARROW still holds on the declared pairs, but X_J(∞) moves. The temperature ratio where X_J(∞) changes sign, 1 − r²/2, moves from 0.73 to 0.93 (0.83 at pin = 0.3). For the undeclared near-equilibrium pair (0.9, 1), X_J is +0.069 at pin = 0.3 and −0.027 at pin = 1.0. |
| **E-20** SYN-0 architecture | **IND** | none | all nine layers (described, not valued) | An architecture/ontology statement. |
| **E-21** sector-selection S-1, PROVISIONAL (IN CLASS) | **CON** | (X) the low-frequency (spectral-bottom) exponent *class* of the influence data; per-branch increments (+2 per vertex power, +(d−1) phase space, +4 per order-2 cancellation); min-dominance over branches | S-1 (locality, d); S-6 (coupling class / symmetry, spin-2, tracelessness, when applied to gravity); S-7/S-11 (amplitude, state, boundary: "sector-supplied") | The record says it is "a class/compatibility selection statement, not a unique point-selection law", and "stop saying 'counting selects the exponent'". **Two models (W3):** both obey the same per-branch law (exponents s, s+4). With base amplitude a > 0 the ω→0 exponent is s; with a = 0 it is s+4. Finite windows interpolate (slopes 3.0 → 6.98 as a goes from 1 to 10⁻² on [0.45, 0.9]). Masking in a 2×2 matrix channel likewise depends on the off-diagonal amplitude. The selection is *among* supplied branch classes (CON). |
| **E-22** CA-1 earned elimination of the classical carrier | **CON** | (S) J(ω) > 0 at ω > 4T/ħ (the relaxational response has no upper edge) vs the zero classical noise floor | S-3 (quantum class); S-4 (ħ); S-11/S-13 (T); S-1 (ring K) | **Two models (W4):** the phonon (z = 1) and the Schrödinger boson (z = 2) on the *same* K both survive. Their low-ω J exponents are −1.01 vs −0.50 (1D theory −1, −½) and their upper edges are 2 vs 4. The cut removes one option and fixes no edge value. The classical-carrier violation threshold ω > 4T reproduces the record (T = 0.5: admissible at ω = 1, violating at ω = 3). |

**Tally:**

| label | rows |
|---|---|
| FIX | none |
| READ (primary) | E-1, E-12, E-13, E-14, E-19 |
| CON (primary) | E-4, E-9, E-21, E-22 |
| SUP | E-11 |
| IND / SUP | E-10 |
| IND | E-2, E-3, E-5, E-6, E-7, E-8, E-15, E-16, E-18, E-20 |
| NF / IND | E-17 |

E-2 and E-15 become READ, and E-5 has two READ entries, under the extended (M) class.

---

## 2. Two-model witnesses (fresh code; numbers from the logs)

### W1 — E-13 against the corrected quotient

Files: `v0_4_sr_w1_sf1_edge_quotient.py` / `.log`.

**Shared predicate.** In every parent, class(dense) ≠ class(dilute).

**Differing primitive.** The one-particle dispersion ε(k), which comes from S-1/S-7. The sector also differs
inside a parent (S-7).

| parent | sector | edges (L), \|v\| (V), r_e (J) | soft set (P), quotient [0, π] | z | class |
|---|---|---|---|---|---|
| A: −2cos k | D, ν = ½ | ±π/2, \|v\| = 2, r_e = 1 | {0, π} | 1.000 | (1, 2) |
| A | D-¼ | ±π/4, \|v\| = √2, r_e = 1 | {0, π/2} | 0.999 | (1, 2) |
| A | E, N = 1 | 0, v = 0, r_e = 2 | {0} | 1.999 | (2, 1) |
| B: −2cos k + ½cos 2k | E, N = 1 | 0, v = 0, **r_e = 4** | {0} | **3.998** | **(4, 1)** |
| B | D-¼ | ±π/4, \|v\| = **0.413** | {0, π/2} | 0.996 | (1, 2) |
| C: −2cos k + 2cos 2k | D, ν ≈ ½ | **4 edges** (±0.427, ±1.993), \|v\| ∈ {2.19, 4.81} | **{0, 0.85, 1.57, 2.29, 2.42}** | 0.998 | (1, 5) |
| C | E, N = 2 | 2 pockets at ±arccos ¼, v = 0, r_e = 2 | **{0, 2k₀ = 2.636}** | 1.98 | **(2, 2)** |

**Effective laws differ:**

- ω⁻ ∝ q² vs q⁴ in the dilute sector;
- one vs two inequivalent edge velocities in the dense sector;
- 1 vs 2 soft points in the dilute sector.

**Inside parent A:** D and D-¼ share the class but differ in edge position and velocity. This confirms the class is
a strict compression of the edge data.

### W2 — E-1, E-14 and E-19 on the pinned conservative chain

Files: `v0_4_sr_w2_pinned_chain.py` / `.log`.

The record's numbers are reproduced at pin = 0.3: r = 0.582109; X_J = 1.16943, 9.16943, 0.33057, 0.73057;
κ = 0.329690 = 1/(2√2.3).

**Part (a): vary the spectral bottom (S-1 pin).**

| pin | r | κ | sign change of X_J(∞) at T_s/T_b |
|---|---|---|---|
| 0.1 | 0.7298 | 0.3450 | 0.734 |
| 0.3 | 0.5821 | 0.3297 | 0.831 |
| 0.6 | 0.4693 | 0.3101 | 0.890 |
| 1.0 | 0.3820 | 0.2887 | 0.927 |

- The predicates are unchanged at every pin: NET-ARROW on the declared pairs, M-1 YES / M-2 NO, and P_memory.
- E-1 factorization: k_pin/(e^{−pin τ}k₀) = 1.0 at every τ. The late decay rate is pin + O(log τ/τ).

**Part (b): vary the edge exponent** at fixed band [0.3, 4.3], with a nearest-neighbour, conservative chain
(Gegenbauer-graded couplings near the end).

| α | end-site measure | native tail exponent | theory |
|---|---|---|---|
| 1 | (1−x²)^{½} | 1.492 | 3/2 |
| 2 | (1−x²)^{3/2} | 2.478 | 5/2 |

"Non-Markovian dissipation only, branch-cut tail" holds in both cases. The record's t^{-3/2} is a read-out of the
supplied chain's square-root edges.

### W3 — E-21 per-branch counting with masking

Files: `v0_4_sr_w3_branch_masking.py` / `.log`.

- **Shared:** branch exponents s = 3 and s+4 = 7.
- **Varied:** the amplitudes (sector-supplied).
- **ω→0 exponent:** 3 for a > 0 and 7 for a = 0.
- **Window slopes:** on [0.1, 0.2] the slope is 3.002 / 3.200 / 6.087 / 7.000 for a = 1 / 10⁻² / 10⁻⁴ / 0. In
  the matrix channel the slope depends on the off-diagonal amplitude c (subject to cone positivity c² < ab).

### W4 — E-9 and E-22 (cone and carrier cut)

Files: `v0_4_sr_w4_carrier_cone.py` / `.log`.

- **Classical relaxational carrier:**
  - T = 0: ν − J/2 = −0.15 at ω = 1 (violates the cone);
  - T = 0.5: +0.15 at ω = 1 (admissible), −0.058 at ω = 3 (violates);
  - general threshold: ω > 4T.
- **The two survivors on the same K:**

  | carrier | J ~ ω^x at low ω | upper edge |
  |---|---|---|
  | phonon | x = −1.01 | 2 |
  | Schrödinger | x = −0.50 | 4 |

- **E-9 alone:** two cone members with spectral-bottom exponents 1 vs 3.

**Not run numerically:**

- E-4: the analytic directed-ring witness in the table suffices.
- E-12: the analytic d-dependence of 2d suffices. The record itself supplies the window witness (6.2142 vs 6.0454).

---

## 3. Headline verdict

> **No row is FIX.** Every row that involves an edge or spectral datum either reads it from a supplied primitive
> (E-1, E-12, E-13, E-14, E-19; plus E-11 as transported input) or constrains its sign, existence or selection
> while leaving its value free (E-4, E-9, E-21, E-22). The other rows carry no edge datum (IND), and E-17's
> question cannot be posed at fixed invariants (NF).

**Why no FIX can exist, by pattern.** This is a sketch of the impossibility, not of a FIX.

1. **Ingredient → property implications (E-1, E-2, E-4, E-5, E-9, E-22).** The banked content is "if supplied
   ingredient X, then property P", or "X is necessary for P". A model without X satisfies the implication
   vacuously and can carry any datum value. A model with X carries a datum that varies with the size of X. So the
   predicate never isolates a value. The most it isolates is a sign or existence bit, which is CON, and that bit
   is itself a function of the supplied antecedent.
2. **Parent-scoped theorems (E-12, E-13, E-14, E-19).** Each is proved for one declared parent, or one family. The
   datum is an explicit function of that parent's supplied K / ε / pin / d / sector. Perturbing the parent inside
   the class the theorem's *mechanism* covers keeps the banked predicate and changes the datum (W1, W2). The
   record forbids promoting these theorems beyond their parent scope, so the predicate has no content that could
   pin the datum independently of the parent.
3. **Class or selection statements (E-21, and E-13's CON aspect).** These pick a class among options generated by
   supplied structure. Amplitudes, state and sector, all supplied, decide which class is observed (W3).
4. **IND, SUP and NF rows** contain no edge value to fix.

**Closest calls (adversarial search; each fails the decisive test):**

- **E-13, dilute sector: v_b = 0 for every smooth parent.** The dilute edge sits at a band extremum, so its
  velocity is zero for any parent. This is the strongest "universal value" in the record. It still fails:
  - it is conditional on the *supplied* sector scaling (fixed N), and the dense sector gives v_b ≠ 0;
  - the datum that decides the law under the corrected quotient, r_e, is not fixed: it is 2 for parent A and 4
    for parent B (W1).

  At most it is a conditional determination of a compressed datum: READ, not FIX.
- **E-4: the ¬CM bit.** It is fixed only once affinity ≠ 0 is supplied, so it is a function of the supplied
  antecedent (READ / CON).
- **E-21: the +4 increment per order-2 cancellation.** The increment is fixed per branch, but the symmetry order
  is supplied, and the *observable* edge exponent depends on supplied amplitudes (W3).
- **E-12: "2d = 6 remains the dimensional contribution".** It is fixed once d = 3 and the acoustic edge are
  supplied (READ). d = 2 gives 4.
- **E-1: rate = pin, exactly.** This is an identity between two forms of the same supplied datum (READ).

Under the rejected literal reading of FIX (§0), E-12, E-13, E-14 and E-19 would be "member-scoped FIX". That
reading would also make every numerical result in the record a FIX, which empties the label. I record it as an
interpretation fork, not as a FIX finding.

---

## 4. E-13 against the corrected V0-3 edge-data quotient

**Corrected quotient (spec §3).** The relevant object is the edge-data set:

- positions and orientations of the relevant occupation edges;
- their local dispersion orders / jets r_e;
- soft momenta from the inequivalent allowed edge differences;
- z and the soft-point structure as separate invariants.

The scalar v_b = 0 / ≠ 0 test is only a compression for the single-pocket cosine-like case.

**Classification of E-13: READ (primary), CON (secondary). Not FIX.**

**READ.** SF-1's class (I-z, I-q) is a function of the edge-data set, and every element of that set is computed
from supplied primitives:

- edge positions k_e from ν (sector, S-7) through exclusion (S-3/S-7; P-02 shows exclusion, not exchange sign, is
  the carrier);
- jets r_e and velocities from ε(k), the one-particle hopping of the supplied ring (S-1);
- soft momenta from the differences k_a − k_b, under the SF-1 BZ quotient (an owner-ruled convention).

The record itself says "the IR class follows once the sector scaling is supplied" and "a supplied sector does not
imply a supplied effective law". The second sentence is about the *law* being derived given the sector. It is not
a claim that the edge data are derived.

**CON.** The banked finding constrains the map from sector to edge region: "the conserved sector selects which
region of the one dispersion becomes the IR reference structure". This is a selection among regions (edges)
supplied by ε. It is selection among supplied options, not value-fixing.

**What SF-1 actually records, in quotient terms.** For parent A, SF-1's (z, n_soft) is a faithful compression of
the edge-data set *up to the class*:

- z = r_e at the active edges (1 at Fermi points, 2 at the band bottom);
- n_soft = |{0, 2k_F}/quotient|.

It does not record:

- edge positions (π/2 vs π/4 in D vs D-¼);
- velocities (2 vs √2);
- the orientation data.

So even inside its own scope it compresses, rather than fixes, (L) and (V).

**What the corrected quotient adds against the single-v_b formulation:**

1. **The jet r_e is independent of v_b = 0.** Parent B has v_b = 0 like parent A, but r_e = 4 and z = 4, so the
   class is (4, 1) and not (2, 1). The scalar v_b test cannot see this.
2. **Multiple pockets.** Parent C at ν ≈ ½ has four occupation edges with two distinct |v| (2.19 and 4.81) and
   five soft momenta. The dilute parent C has two pockets: v = 0 at both, r_e = 2, and soft set {0, 2k₀}, giving
   class (2, 2) and not (2, 1). One v_b cannot represent either case.
3. **C-6.** The two vanishing scales (k_F ~ L^{−1/2}, q ~ L^{−1}) make z depend on the path. The edge jet at the
   moving edge, together with the order of limits, is what the quotient needs. This is consistent with SF-1's own
   MULTISCALE report-only label.

**Verdict on E-13.** Against the corrected quotient, E-13 READs the edge-data set (positions, jets and soft
differences supplied by ε, ν and exclusion) and CONstrains which region of the band is IR. It FIXes none of it:
W1 gives three parents with the same banked predicate and different edge data and law classes. The P-02 single-v_b
theorem was used only as a pointer to the mechanism and was not relied on. Its "law class = function of v_b^∞"
fails exactly where the corrected quotient says it does: parent B (jet) and parent C (pockets).

**Ambiguity, both readings given:**

- If "edge data" is restricted to the class actually banked, (z, I-q), then E-13 is still not FIX. The class is
  determined only after the sector and the parent are supplied, and W1 parent B changes it.
- If "SF-1" is read as the single member ε = −2cos k, as frozen, the class values for that member are
  member-scoped computations. That is the literal-reading fork of §0, giving READ (member-scoped "FIX*").

---

## 5. Upstream trace summary (READ / CON / SUP rows)

| row | primitive(s) the datum traces to |
|---|---|
| E-1 | S-1 (pin, C1-a K); S-2 (relaxational generator) |
| E-4 | S-2 (affinity / orientation) on S-1 (ring) |
| E-9 | S-4 (ħ); S-3 (quantum class) |
| E-11 | S-13 |
| E-12 | S-1 (lattice, d, BC, dispersion) and the declared window |
| E-13 | S-1 (ε); S-3/S-7 (exclusion statistics, conservative realization); S-7 (sector) |
| E-14 | S-1 (K, pin, K₁₁ = ω_s²); S-2 (declared conservative parent); the admitted L-vH limit (a declared deformation, not in Table 3, priced by "conditional on the admitted L-vH deformation") |
| E-19 | S-1; S-11 (preparation, T_s, T_b); S-12 (declared observables / bookkeeping) |
| E-21 | S-1 (locality, d); S-6 (symmetry / coupling class); S-7/S-11 (amplitudes, state) |
| E-22 | S-3; S-4; S-11/S-13 (T); S-1 |

**Trace gap.** Table 3 has no explicit row for "admitted limits / deformations", which matters for E-14. It also
has no row for "declared one-particle dispersion of a non-GRUT parent", which matters for E-13: SF-1's
free-fermion parent is not GRUT's substrate. I map that parent to S-1 and S-7 ("conservative realization") and
flag the mapping as interpretive.

---

## 6. Rubric-adequacy note

1. **FIX is under-specified.** "With no supplied primitive varied" does not say *within what class* the primitives
   are held. On a literal member-level reading, every parent-scoped computation (E-12, E-13, E-14, E-19) is FIX.
   The decisive two-model test fixes this and was used throughout. Recommendation: define FIX as "invariant over
   all models satisfying the banked predicate within its declared mechanism class".
2. **READ and CON are not mutually exclusive.** Necessity implications (E-1, E-4) both *read* a supplied datum
   (the conclusion is a function of it) and *constrain* the admissible data (gap > 0, ¬CM). E-13 and E-14 have
   the same overlap. Dual labels are needed and were used. **Flag: E-1, E-4, E-13 and E-14 depend on whether the
   primary emphasis is the functional dependence or the restriction.**
3. **SUP vs READ.** E-11 transports the input (SUP) but its output is a function of it (READ). E-10 is IND or SUP
   depending on whether "holding (K, N) fixed" counts as involving a datum.
4. **NF vs IND.** For E-17, the adjudicated crossed cell is NF, while the report-only R₂ note is a READ identity.
   The label depends on whether report-only notes count.
5. **Datum classes.** The strict five classes miss the edge jet r_e, which V0-3 needs, and low or inverse moments.
   E-2, E-15 and two E-5 entries change from IND to READ if (M) is admitted. E-12's 6.2142 is a window slope, not
   an edge exponent: only its "2d" part is an edge datum, and the +0.2142 is window curvature.
6. **"Relevant to an effective law".** E-2 (geometry) and E-8 (partition) concern no effective law, which supports
   IND whatever the spectral reading.
7. **Upstream labels.** Table 3 lacks a row for admitted limits (L-vH) and for non-GRUT declared parents. The
   S-labels for E-13 and E-14 are therefore partly interpretive (§5).

**Net assessment.** The labels are well defined enough to decide the headline (no FIX) robustly under every
reading considered, except the rejected literal reading. Row-level primary labels are interpretation-sensitive for
E-1, E-2, E-4, E-5, E-10, E-11, E-13, E-14, E-15 and E-17, in the READ/CON, IND/READ, SUP/READ and NF/IND
directions. None of those forks changes a FIX status.

---

## 7. Firewall incidents

No allowed file quoted an edge-data-audit classification matrix. No reading was stopped for that reason.

Directory listings (`ls`) of `/home/user/VER0`, `playground/`, `playground/SCOUT_0/` and `playground/SCOUT_0/probes/`
showed deny-listed *filenames*. None was opened.

A content grep for the literal string `0.678` (not a forbidden string) printed one-line matches from three other
files (§8). None was opened beyond those lines.

**Interpretation of "cited by Table 1/3".** I read verdicts and maps that a Table-1 ruling directly rules on or
cites (one hop). All are top-level primary records, none is under `verification_0/` and none is deny-listed.

---

## 8. Files read

**Full:**
- `verification_0/specs/V0_4_EDA01_TARGET.md`
- `playground/SCOUT_0/BASELINE_MAP.md`
- `playground/SCOUT_0/probes/P02_RESULT.md` (E-13 only)
- `SF1_OWNER_RULING_01.md`, `SF1_OWNER_RULING_02.md`, `SF1_FORMATION_VERDICT_01.md`
- `S5_OWNER_RULING_01.md`, `S5_OWNER_RULING_03.md`
- `S6_OWNER_RULING_02.md`
- `GR2_L6_OWNER_RULING_01.md`, `GR2_L6_VERDICT_01.md`
- `P3_P4_OWNER_RULING_01.md`
- `SYN1_OWNER_RULING_01.md`
- `L0_1A_OWNER_RULING_01.md`, `L0_1A_VERDICT_01.md`
- `L0_1B_OWNER_RULING_01.md`, `L0_1C_OWNER_RULING_01.md`, `L0_1D_OWNER_RULING_01.md`
- `L0_1_FLOOR_DEPOSIT_01.md`, `L0_1_FLOOR_O7_OWNER_RULING_01.md`
- `L0_ACCESS_BRIDGE_OWNER_RULING_02.md`
- `LEVEL0_FOREST_SYNTHESIS_OWNER_RULING_01.md`
- `P2_S1_OWNER_RULING_01.md`, `D1_P2_OWNER_RULING_01.md`
- `GR2_SYNTHESIS_OWNER_RULING_01.md`
- `KERNEL_TRANSPORT_OWNER_RULING_01.md`, `CLOCK_MISMATCH_OWNER_RULING_01.md`
- `S2_OWNER_RULING_03.md`
- `S3_OWNER_RULING_01.md`
- `L0_LIFT_SELECTION_OWNER_RULING_02.md`
- `L0_1G_OWNER_RULING_01.md`

**Partial:**

| file | portion read |
|---|---|
| `S2_HB_OWNER_RULING_02.md` | lines 1–60 |
| `L0_STRUCTURAL_DIAGNOSTIC_REVERSAL_OWNER_RULING_01.md` | lines 1–70 |
| `CA1_CARRIER_VERDICT_01.md` | lines 1–80 |
| `RELATIONAL_ONTOLOGY_MAP_01.md` | lines 1–60 |
| `GRUT_WORKING_THEORY_SYNTHESIS_01.md` | heading grep; lines 237–406 |
| `SFG0_OWNER_RULING_01.md` | lines 40–89 |
| `SECTOR_SELECTION_VERDICT_01.md` | lines 1–60 |
| `P3_NC_LIFT_VERDICT_01.md` | grep context only: about lines 41–55 and 111–124 |
| `L0_1B_VERDICT_01.md` | lines 1–60 |
| `uploads/GRUT_Consolidated_Theory_PUBLIC_RECORD.md` | one line (1046), via the `0.678` grep |
| `LEVEL0_FOREST_SYNTHESIS_01.md` | one line (401), via the same grep; line 169 matched but was omitted as too long |
| `GRUT_WORKING_THEORY_01.md` | one line (216), via the same grep |

**Metadata only (no content):** `wc -l` line counts of about 30 top-level record files; directory listings as in §7.

**Not read:** every deny-listed file; every other `probes/*_RESULT.md` (allowed, but not needed); everything else
under `verification_0/`; the other working trees. No git commands were used.

## 9. Files created

All paths are under `/home/user/VER0/verification_0/`.

- `V0_4_EDA01_SECOND_READER.md` (this file)
- `code/v0_4_sr_w1_sf1_edge_quotient.py`, `code/v0_4_sr_w1_sf1_edge_quotient.log`
- `code/v0_4_sr_w2_pinned_chain.py`, `code/v0_4_sr_w2_pinned_chain.log`
- `code/v0_4_sr_w3_branch_masking.py`, `code/v0_4_sr_w3_branch_masking.log`
- `code/v0_4_sr_w4_carrier_cone.py`, `code/v0_4_sr_w4_carrier_cone.log`
