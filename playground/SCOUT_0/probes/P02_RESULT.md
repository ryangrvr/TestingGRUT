# SCOUT_0 W3 P-02 RESULT — the SF-1 statistics swap (one controlled pass)

**Charter:** `PROBE_CHARTERS.md` P-02, sharpened by the auditor. Control map:
`P02_SF1_CONTROL_MAP.md` (every SF-1 ingredient extracted from the record; the only intended change is
fermionic → bosonic statistics). Script: `p02_boson_sectors.py` (log `p02_boson_sectors.log`).
NEW HYPOTHESIS CLASS — DOES NOT ALTER OLD TERMINAL. SF-1 and SFG-0 are untouched; nothing migrates
into earned Level-0.

## 0. Verdict

> **Main model: BOSON-COLLAPSE** (charter terms: STATE-NOT-LAW). Free bosons give the class
> `(I-z, I-q) = (2, 1)` in **every** SF-1 sector family — D, E, D-¼, D-¾, E-3, Ē and C-6 — under both
> prescriptions. The collapse is exact and stronger than class level: the normalized object
> `S_{L,N}(q, ω)` is **independent of `N`** at every `L`.
>
> **Minimal carrier: occupancy exclusion, not exchange antisymmetry.** Hard-core bosons (commuting
> operators, `n_j ≤ 1`) reproduce the SF-1 support **exactly** (all odd-N sectors checked at L = 6,
> 10, 12). Statistics in the sense of exchange sign is not load-bearing; exclusion is.
>
> **Sector invariant → law class (theorem, free number-conserving single-band parents):** the IR
> class is fixed by where the sector's ground-state **occupation edge** `k_b^∞` sits on the dispersion.
> `z = 1` iff the edge velocity `ε′(k_b^∞) ≠ 0` (edge inside the band); `z = 2` iff `k_b^∞` is a
> band extremum. Statistics enter **only** through the map `ν ↦ k_b^∞`: exclusion gives
> `k_b^∞ = πν` (it moves with density); free bosons give `k_b^∞ ≡ argmin ε` (pinned at the band
> bottom for every density).

That is the auditor's anticipated carrier, made exact: **how a conserved density forces the occupation
distribution against the spectral edge.**

## 1. Exact model and checks

- **Parent:** `H_B = −Σ_j (b_j†b_{j+1} + h.c.)` on the periodic ring, full bosonic Fock space, conserved
  `N`. No μ, no sector term, one builder. Readout, object, prescriptions and invariants are
  identical to SF-1 (§§4–5).
- **V-BFOCK / A-BOSE** (exact sector blocks, L ∈ {6, 10}, N ∈ {1, 3, 5}, dimension up to 2002):
  the sector ground state is unique (gap 1.000 at L = 6, 0.382 at L = 10), with `E₀ = Nε(0)` exactly,
  i.e. the condensate `(b₀†)^N|0⟩/√N!`. The grouped `ρ_q` support is the single line
  `ω = 4 sin²(q/2)` for **every** q — all True.
- **Closed form (every family):** `ω⁻(q) = 4 sin²(q/2)`, `z_P = z_{P′} = 2`, soft set `{0}` ⇒ class
  `(2, 1)`. The SF-1 fermion classes are `(1, 2)` for D-type and `(2, 1)` for E-type.

**Proof of the collapse (all L, all N).** For free bosons `H_B = Σ_k ε(k) n_k`. In sector `N` the
minimum `Nε_min` is reached only by putting every boson at `argmin ε = {0}` (unique on the ring), so
the ground state is the condensate. Then `ρ_q = Σ_k b†_{k+q} b_k` acts only through `k = 0`:
`ρ_q|N₀⟩ = √N |(N−1)₀, 1_q⟩`, at energy `ε(q) − ε(0)` with normalized weight `N/N = 1`. Neither
quantity depends on `N`. ∎

## 2. Zero-mode firewall

Classification: **not UNCONTROLLED.** The zero mode `k = 0` is the band minimum. With `N` exactly
conserved there is no divergent occupation, no μ pinning ambiguity and no degeneracy: the condensate
is the unique sector ground state at every L (checked), and the IR limit is regular. The zero mode
**is** the collapse mechanism, not an artifact to be removed. No regulator, gap or zero-mode
separation was needed. The listed controls were therefore not run, as unnecessary rather than
skipped: the main comparison is well posed exactly.

## 3. Minimal-carrier test (assumptions removed one at a time)

| Variant | Exclusion? | Antisymmetry? | Result | Reading |
|---|---|---|---|---|
| Free fermions (SF-1 record) | yes | yes | split `(1,2)` vs `(2,1)` | — |
| **Free bosons (main)** | no | no | **collapse** `(2,1)` everywhere | conservation alone is insufficient |
| **Hard-core bosons** (NEW ASSUMPTION, diagnostic) | yes (`n_j ≤ 1`) | **no** | **identical to SF-1** for every odd N (L = 6, 10, 12) | exclusion suffices; exchange sign is not needed |
| **Bosons in the Fermi pattern** (diagnostic) | no | no | an exact eigenstate but not the ground state (`E − E_gs = 372.8` at L = 1026); support has `ω < 0` (`min ω = −1.99, −0.76, −0.20` at `q = π/2, π/8, π/32`) | exclusion is used **twice**: to make the edged occupation the ground state, and to block transitions into occupied levels |

- **HCB scope note.** In 1D, hard-core bosons map to free fermions by Jordan–Wigner with boundary twist
  `(−1)^{N−1}`. That is periodic for odd N (SF-1's convention), and the density readout is
  JW-invariant. This is why the match is exact. Even N gives the twisted ring (L = 12, N = 6: not
  used). In d > 1 hard-core bosons are not free fermions, so the HCB statement is **1D-scoped**.
- **Fence honoured.** HCB is an interaction (infinite on-site repulsion). It is used only to isolate
  the carrier, never to overwrite the free-boson verdict.

## 4. Theorem (sector invariant → law class)

**Class.** Free (quadratic), number-conserving, single-band parent `ε(k)` on a 1D ring, with
smooth `ε`. Reference state: the sector ground state. Readout: density; SF-1 prescription P.

**Statement.** Let `O_N` be the ground-state occupation set and `k_b(L)` its edge momenta, with limit
`k_b^∞`.

- **(i)** The support is `{ε(k+q) − ε(k) : k ∈ O_N, k+q admissible}`, where admissible means
  unoccupied under exclusion, and anything for free bosons.
- **(ii)** **Exclusion:** `O_N` = the N lowest levels, so `k_b^∞ = πν` for `ε = −2cos k`. At a
  regular edge, `ω⁻(q) ≈ |ε′(k_b^∞)| q`, so `z = 1`; at an extremal edge, `z = 2`. `I-q` counts the
  soft momenta `{0, 2k_b^∞}` modulo the quotient.
- **(iii)** **Free bosons:** `O_N = {argmin ε}` for every N, so the class is that of a band extremum
  for every sector: one class, collapse.

**Hence:** the law class is a function of the single sector invariant `k_b^∞` (equivalently the edge
velocity `v_b = ε′(k_b^∞)`). Particle statistics act only through `ν ↦ k_b^∞`.

- **Record reproduced** (script table): D, D-¼, D-¾ have `v_b = 2, √2, √2 ≠ 0` ⇒ `z = 1`; E, E-3, Ē
  have `v_b = 0` ⇒ `z = 2`. All six match the SF-1 verdict.
- **C-6 explained.** `k_b(L) ~ L^{−1/2}` gives `v_b(L) → 0` at a rate that competes with `q ~ L^{−1}`.
  The path dependence is the two-scale approach of `v_b` to 0.
- **Proof:** (i) is the free-particle form factor of `ρ_q`; (ii) and (iii) follow from the ground-state
  occupation of each statistics and a Taylor expansion of `ε` at `k_b^∞`. Grade: exact for the
  declared class (1D, one band, free, sector ground state, SF-1 prescription).
- **Classification:** REDISCOVERED-KNOWN as physics (Fermi points vs band bottom; free Bose
  condensate has quadratic single-particle excitations). KNOWN-BUT-NEW-IN-GRUT: the reduction of the
  SF-1 mechanism to one sector invariant, with exclusion — not exchange statistics — as the minimal
  carrier.

## 5. Universality alarm — answered

Fermions and hard-core bosons share the sector → law-class map. Free bosons do not. The quotient
that decides the class is **`k_b^∞` relative to the band extrema**, i.e. **whether the conserved
density produces a nonzero IR velocity scale at the occupation edge**. The class can be written
without mentioning statistics:

> `class = (z, n_soft)` with `z = 1 ⇔ v_b^∞ ≠ 0`, `z = 2 ⇔ v_b^∞ = 0` (quadratic extremum), and
> `n_soft` = number of inequivalent edge-pair momenta.

**Conjecture (spawned, not run — interaction fence):** a sector law-class split needs a mechanism that
converts conserved density into a nonzero IR velocity. Exclusion is one (`v_F = ε′(πν)`). Repulsive
interaction is another: weakly interacting bosons have Bogoliubov sound `c = √(gn)`, which vanishes
as `n → 0`, so dense and dilute sectors would split again (`z = 1` vs `z = 2`). This is textbook
physics, but it would make the carrier **"density → IR velocity"** rather than any specific
statistics. Spawned as **P-02b**. Not authorized in this pass.

## 6. What this changes

- **SF-1's dependence is pinned down:** conservation is necessary but **not sufficient**. The split
  needs exclusion (or, per the conjecture, another density → velocity mechanism).
- **SFG-0's Tier-C qualifier sharpens (candidate, owner's call):** "conditional on supplied
  conservative **fermionic statistics**" can be narrowed to "conditional on supplied conservative
  **exclusion**" in 1D. The tier and terminal are unchanged.
- **For the target set in ZOOM_OUT_02:** this *is* a mechanism that restricts the effective law —
  `sector invariant v_b^∞ → law class` — but the restricting power sits in a **supplied** ingredient
  (exclusion, or interaction), not in conservation alone. That is consistent with the quotient
  principle: the law class factors through `v_b^∞`, and statistics are fibre data except where they
  move `v_b^∞`.

**Status: P-02 COMPLETE (one controlled pass). Main verdict BOSON-COLLAPSE (exact, all sectors).
Carrier = exclusion (HCB = SF-1 exactly, 1D). Theorem: law class = function of the occupation-edge
velocity `v_b^∞`. P-02b (interacting bosons) spawned, not run. Next: P-01.**
