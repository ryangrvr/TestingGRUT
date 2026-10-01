# SCOUT_0 W1 P-09 RESULT — locality of representation (one reconnaissance pass, budget-capped)

**Charter:** `PROBE_CHARTERS.md` P-09 (frozen). **Budget (auditor, binding):** one decisive pass, no
definition-tuning. Question as sharpened by the auditor:

> Is there any locality notion already fixed strongly enough by the earned substrate and recorded
> lift maps to prune the lift set **without** introducing another representation/algebra/net choice?

Branches that need operator-vs-Poisson, tensor-vs-graded locality, `*`-vs-non-`*` images, a new
quantum local net, a new factorization, or new access/readout structure are marked
**NEW ASSUMPTION / CONDITIONAL** or **BLOCKED** and stopped.
**Script:** `probes/p09_locality_checks.py` (pure Python; substrate = the recorded anchor bath block
`K_b`, `build_K(24,0)` `(1:,1:)`, `N = 23`; small abstract models labelled as mathematics).
**Governing texts:** `L0_LIFT_SELECTION_01.md` §1, §3 (D-6), §5 (R-0, R-1);
`L0_LIFT_SELECTION_OWNER_RULING_01.md` (D-6 = CRITERION, item C); `..._VERIFICATION_01.md` V-1, V-6,
V-7; `..._CORRECTIONS_01.md` LS-8, LS-9; `L0_LIFT_SELECTION_OWNER_RULING_02.md` §3–4;
`EA0_OWNER_RULING_02.md` §3–4; `L0_1B_OWNER_RULING_01.md` (E-2).

## 0. One-line verdict

> **No.** The only locality the earned record fixes is (i) the declared commuting site net
> (commuting classical site coordinates) and (ii) the support graph of `K_b`. Both are
> **descended one-particle data** that every lift reproduces, so as a test they are passed by all
> four lifts and prune nothing. Every locality notion that *could* prune needs a choice the record
> does not make: a bath factorization (Sz.-Nagy), a site-local complex structure (Λ-B), tensor
> vs graded locality (Λ-F), or generator k-locality (already ruled CRITERION-SMUGGLED).
> **P-09 reinforces P-08: selection at the lift interface is controlled by representational
> choices the earned theory does not fix.**

## 1. Frozen definition (stated before evaluation)

**L-E (earned locality).** A lift is *locality-compatible* if (a) the images of the site coordinates
lie in site-labelled subalgebras whose physical (observable) parts commute across disjoint sites,
and (b) its one-particle generator has off-diagonal support inside `G(K_b)`, the support graph of
`K_b`.

- (a) is the earned fact behind D-6 ("commuting classical site coordinates; the declared net",
  `L0_LIFT_SELECTION_01.md` §3). The net itself is **declared, not derived** (`EA0_OWNER_RULING_02.md` §3).
- (b) is the geometry support already inside D-1 (E-2: locality → P_geometry, NECESSITY-CERTIFIED).

**Record fence found on reading (decisive by itself).** R-0 is settled and forbids any discriminator
outside D-1…D-7. A "locality of lift" selector **is** D-6, and the owner ruled D-6 **CRITERION**:
"usable as a compatibility check only". So at earned scope P-09 *cannot* be a selector: it can only
re-run D-6 as a compatibility check. This pass does exactly that, for the record, on all four lifts.

## 2. Earned-scope check (L-E) on all four lifts

| Lift | (a) commuting site net | (b) generator support ⊂ `G(K_b)` | L-E |
|---|---|---|---|
| Λ-H cotangent (`H = pᵀf`) | site Poisson algebras `{x_i, p_i}` commute across sites (canonical bracket) | `∂²H/∂p_i∂x_j = −(K_b)_ij − 12βx_i²δ_ij`: support = `G(K_b)` exactly (C1; all β) | **passes** |
| Λ-H Sz.-Nagy | the source-descended part is the pullback `f∘P_H` (classical function reading) — commutative | system part of the generator is `−K_b` (support `G(K_b)`); the **bath** has no earned net, so (b) is not formulable for the coupling | **passes on the system; bath part NOT-FORMULABLE** |
| Λ-B complex | under the **recorded** doubling (`π_i` conjugate to `x_i`, V-7/LS-9) `σ(x_i, x_j) = σ(x_i, π_j) = 0` for `i ≠ j`: site Weyl algebras commute (C3) | `dΓ(K_b) = Σ (K_b)_ij a_i†a_j`: hopping on `G(K_b)` | **passes** (under its already-priced doubling) |
| Λ-F complex | graded: odd generators anticommute across sites, the **even** (physical, parity-superselected) net commutes (C4; V-6, LS-8) | even hopping on `G(K_b)` | **passes** (graded reading; D-6 ruling) |

**Why this was predictable (the structural reason).** (a) and (b) are functions of the descended
one-particle data `(K_b, site labels)`. Every lift reproduces `e^{−K_b t}` on its one-particle
sector (I-1 as verified, V-1), so every predicate built only from those data takes the **same value
on every lift**. Finite-time locality is shared too: `−K_b` is Metzler and irreducible, so
`e^{−K_b t}` is entrywise positive for all `t > 0` in every lift alike (C2 tail: `[e^{−20K_b}]_{1,23}
= 4.2·10⁻⁷`). An earned locality selector would have to see structure **not** descended from `K_b` —
and the record lists all such structure as supplied (ruling 02 §4).

## 3. The pruning branches — each needs a new choice → stopped

| Stronger notion | Who could fail | What it requires | Status |
|---|---|---|---|
| Strict locality of the dilation's system–bath coupling | Sz.-Nagy | a **bath factorization** (channel basis). The minimal coupling `C = (2K_b)^{1/2}` is **dense** (529/529 entries; `|C_{1,j}|`: 2.1, 0.50, 6.4e-2, … 3.1e-8 at `j = 23`; decay ≈ 0.7–0.9 per site — quasi-local, not strictly local). A strictly 2-local channel basis exists (`K_b = BBᵀ`, pins + springs, 46 channels) and equals the minimal one up to a channel-space isometry `W` (`WᵀW = I` to 4e-15, `W C_min = C_loc` to 8e-15). Sz.-Nagy is unique only up to unitary equivalence, and **strict locality is not a unitary invariant** | **NEW FACTORIZATION → BLOCKED at earned scope** |
| Tensor-product (ungraded) locality of all field generators | Λ-F | tensor-vs-graded choice: `{a₁,a₃} = 0` but `‖[a₁,a₃]‖ = 2` (C4) | **already D-6 CRITERION** (ruling 01 item C) → stop |
| Site-local complex structure | Λ-B | a choice of `J` on the doubled space: a metric-compatible `J'` (rotation of `x₂` into `π₁`) gives `σ(x₁,x₂) = 0.565 ≠ 0`, so distinct-site coordinate images stop commuting (C3). The recorded doubling convention *is* the site-local choice, already priced (V-7) | **CONDITIONAL on the priced doubling**; adds no pruning |
| Locality as an operator-algebra net (vs Poisson) | cotangent, Sz.-Nagy | operator-vs-Poisson reading (P-08's priced algebra-type choice) plus a new quantum net | **NEW ASSUMPTION** → stop |
| k-locality of the lifted generator | (would test Sz.-Nagy minimal only) | "H is k-local" is **CRITERION-SMUGGLED / EXTERNAL-UNVERIFIED** (`EA0_OWNER_RULING_02.md` §4) | **BLOCKED by ruling**. For the record: even if it were adopted, cotangent, Λ-B, Λ-F pass, and Sz.-Nagy passes or fails according to the bath factorization above, so it still prunes nothing without that choice |
| Access/readout locality (which sites are read) | — | multi-site access sets are supplied (S-10, C-7) | **NEW ACCESS** → stop |

## 4. Verdict (frozen criteria)

**P-09: EARNED SCOPE — FORMULABLE, CONSTRAINT-EMPTY** (unlike P-08, the earned notion *is*
formulable — it is D-6's compatibility check — and all four lifts pass it). **Every pruning branch
needs a new representation/factorization choice → CONDITIONAL / BLOCKED, not pursued (budget).**
Survivor set unchanged: {cotangent, Sz.-Nagy, Λ-B complex, Λ-F complex}. Terminal
IRREDUCIBLE/SUPPLIED untouched; no new selector; charter fence ("no post-hoc locality definition to
kill a survivor") honoured. Null hypothesis of the charter (all survivors local in the declared
sense) **confirmed**.

**Classification.** Mostly **REDISCOVERED-RECORD** (D-6, V-6, V-7, LS-8 already contain the Λ-F and
Λ-B content). One item is **KNOWN-BUT-NEW-IN-GRUT**: the Sz.-Nagy coupling's strict locality is a
property of a chosen bath factorization, not of the (unique-up-to-unitary) dilation — elementary as
mathematics (polar decomposition of `C_loc`), not previously recorded. It adds "bath factorization"
to the supplied-structure list next to P-08's image class and algebra type.

## 5. Hostile attack (brief)

- **H1 (was locality defined too weakly so everyone passes?).** No: L-E is the *strongest* notion the
  record earns; every stronger one is in §3 with the choice it needs named. Tuning definitions to
  manufacture a failure is exactly what the budget forbids.
- **H2 (does the dense minimal coupling disqualify Sz.-Nagy?).** No. Density is a basis property; a
  strictly local basis exists and is isometrically equivalent (C2). Choosing the minimal basis to kill
  it would be a post-hoc factorization choice.
- **H3 (nonlinear members, β > 0).** Locality adds nothing there: the quasi-free lifts fail coverage
  (D-2, CRITERION), and the cotangent generator is graph-local for every β (C1).
- **H4 (standard physics).** Graded fermionic locality, doubling, and quasi-locality of matrix
  functions of banded matrices are textbook. REDISCOVERED-KNOWN as mathematics.

## 6. What changes in the scout map

- **Lift neighbourhood closed for reconnaissance.** P-08 + P-09 together show one mechanism: every
  earned criterion at the lift interface is a function of `K`-descended data, which all lifts share;
  any discrimination needs lift-side structure (image class, algebra type, statistics/grading,
  complex structure, bath factorization), all of it supplied. Banked as a **structural observation**
  for the zoom-out review (`ZOOM_OUT_01.md`), not as a theorem. No further lift probes are queued.
- **Supplied-structure list, additions:** bath factorization (P-09), alongside image class and algebra
  type (P-08).
- **Next:** move outward — P-17 (physical autonomous bath / noise origin).

## Script log (`p09_locality_checks.py`)

- C0: `G(K_b)` = 22-edge nearest-neighbour path; closed-form eigenpairs `λ_k = 2.3 − 2cos((2k−1)π/47)`,
  residual 6e-15, `λ_min = 0.3045`.
- C1: cotangent mixed Hessian off-diagonal support = `G(K_b)` (β = 0.7, generic `x`).
- C2: `C_minᵀC_min = 2K_b` (2e-14); 529/529 entries nonzero; local `B` (46 channels, ≤ 2 sites each),
  `BBᵀ = K_b` (2e-16); polar isometry `W`: `WᵀW = I` (4e-15), `W C_min = C_loc` (8e-15);
  `e^{−20K_b}` entrywise positive, min 4.2e-7.
- C3: recorded doubling `σ(x₁,x₂) = σ(x₁,π₂) = 0`; rotated metric-compatible `J'` (`J'² = −1` exact):
  `σ(x₁,x₂) = 0.5646`.
- C4: 3-site Jordan–Wigner: `{a₁,a₃} = 0`, `‖[a₁,a₃]‖ = 2`; `[hop₁₂, a₃] = 0`.

**Status: P-09 COMPLETE (one pass). EARNED SCOPE CONSTRAINT-EMPTY; all pruning branches
CONDITIONAL/BLOCKED on a named new choice; survivor set and terminal unchanged. Reinforces P-08.**
