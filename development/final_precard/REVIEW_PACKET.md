# FINAL PRE-CARD REVIEW PACKET — Issue #9

**READY FOR EXTERNAL CHECK.** DRAFT — pending an actual external checker.
This is same-author Work evidence. It grants no external clearance, canonical
integration/re-freeze, bank acceptance or Card-1 authorization. No candidate was
formulated, scored or priced: **0 cards / 0 IP-12 drafts**.

The reviewed target is the **isolated evaluator proposal** in the privately
supplied, hash-pinned input archive, together with the **mandatory public
`PrevalidatedKit` boundary** around the unchanged CR-5 module. Direct calls to the
original unvalidated module remain unsafe. The reviewer must approve this exact
boundary or an equivalent in-module repair before integration; a clean replay is
not that approval.

| Check | Result |
|---|---|
| P6 four-document scope and document guards | exact clause removal; 7 PASS; no stale-net cases |
| Frozen, subject-free oracle controls | 24 PASS |
| Appendix-B price comparisons | 1,573 agree; 17,024 finite semantic evaluations |
| Partial HB controls | 14 agree; singular targets rejected |
| NR-4 historical/certified regression | 11 agree; pre-fix false negative reproduced |
| Invalid reference-measure entrypoint checks | 32 PASS across 16 invalid lists |
| Original valid CR-5 comparisons | all 183 unchanged |
| Additional valid weighted-measure comparisons | 240 agree, including zero atoms and VOID |
| Exact review inventory | 2,084 unique entries |
| Remaining supported-profile semantic mismatches | **0** |

## 1. P6 closure evidence

Authority is the [owner disposition](https://github.com/ryangrvr/TestingGRUT/issues/2#issuecomment-6092788257).
Relative to `b89a0e607e71450fce07d5cc55daae8224acad18`, the sole book edit removes
`; the 21 vacuum-cluster nodes were retired from the ledger sum in the V2 consolidation`.
The resulting line retains **net +16 / 53 GRUT-scope nodes**. All current sync
markers still show **74 total = 53 GRUT + 21 vacuum-cluster**, with cluster net 0.
README, Agenda and What Survived are byte-identical to that parent. The historical
P6 amendment was also inspected: the remaining changes synchronize numbers or
label historical values; none supplies a new scientific claim under this ruling.

`P6_NARROW.patch` is the exact four-document diff. `P6_EVIDENCE.json` pins all before/
after hashes, the original amendment-diff hash, empty prose enumerator and seven
document guards. `OPEN_PASSES.txt` closes P6 by **owner ruling plus narrow landed
amendment**, preserving the prior failed review as history.

Closure makes two existing live-integration tests pass without code repairs.
Those precise consequences are hard-coded and checked in `owner_transition.py`;
the original observational manifest is unchanged. No outcome table or bank update
is accepted from the transition file.

Actual default: **234 PASS / 7 FAIL / 1 SKIP**, with zero unapproved deltas and zero
infrastructure errors. The unchanged five-file review-only repair remains
**239 PASS / 2 FAIL / 1 SKIP**; its live classifier now has zero undeclared cases
and zero orphan passes. Five of the original seven repair targets still change
FAIL→PASS in that overlay; the other two already pass after P6 closure. The patch
is not integrated here. All **25 bank flags** and the five remaining OPEN passes
remain unchanged. These record checks are distinct from kit review readiness.

## 2. Price-coder reference and derivation

Normative sources are the canonical frozen charter **§4.1–4.2 / Appendix B** at
`97ef325ccb9471f88243caaa7448e4cc8ad500d7`, blob
`8dd7b05ab1ccec191e1ee8c54fcb41fbfd85c052`. G2-14R overrides the relevant
decoupling clauses, not Appendix B pricing.

`price_oracle.py` was written from that text, tested, and hash-frozen in
`ORACLE_FREEZE.json` **before importing either evaluated implementation**. Its
freeze SHA256 is `5d8d6525eb214b0504766a1d4078f7c7fe88c8acaa5076fd5a814152d27829ab`.
The executable replay verifies that freeze before each comparison.

For a positive integer with binary length L, Elias-delta encodes the length with
gamma code and then L−1 payload bits: delta length is
**L + 2 floor(log2 L)**. All computations use exact integer/string lengths, never
floating logs at boundaries. Natural literals cost delta(n+1); rational literals
cost delta(abs(n)+1) + delta(d+1) + 1, with positive denominator; real literals
cost p + delta(abs(exponent)+1) + 1. There are 64 slots at 6 bits each, including
four unavailable reserved slots and distinct set/arithmetic multiplication slots.
Variable-index width is ceil(log2(v+1)); every explicit variable occurrence is
charged. Definition bodies are charged once, uses pay the definition-index width,
and cyclic or duplicate registries reject.

The normalizer first removes implication/equivalence and formula `ite`, pushes
negation, standardizes binding scopes apart, and pulls quantifiers across Boolean
connectives only with the freshness side condition. Lambda terms stay scoped.
Quantifiers hidden inside terms reject. Renaming preserves free variables and
finite truth controls test capture, shadowing, duplication and alpha equivalence.
The standard prenex equivalence is restricted to the declared **nonempty
first-order domain profile**; no empty-domain or arbitrary typed-term equivalence
is certified by these controls.

The profile explicitly distinguishes raw syntax length from normal-form price;
shorter own-form use needs equivalence evidence. ASCII AST spellings are a review
encoding, not added L0 tokens. Appendix B lists the conditional/process tokens
without AST arities: this packet supplies and hashes explicit fixture signatures,
rejects unregistered signatures, and grants no interpretation/type certificate
for arbitrary uses. Formula macros require hygienic unfolding; the proposed
normal-form entrypoint rejects them. Unsupported paths earn no favorable result.
This is not a claimed complete type checker or automatic semantic IP audit.
Existing CR-3 KD-1…KD-5 duties remain in force.

The suite covers reserved tokens, wrong/missing/extra arguments, unsupported
symbols, ignored literal payloads, negative/fractional integer fields, cyclic/
duplicate definitions, Elias boundaries around powers of two through 2^256,
array-length mismatches, empty structures, and **p=6,10,16**. Table/input/selection
hooks enforce statement prices, explicit table entries and certified prior KL
lower bounds. Synthetic audit IDs test arithmetic only, not proof truth. General
ledger bits can be fractional where log2 prices require it; integer counts cannot.

## 3. HB control derivation

For affine entries Y_i=M_i+G_i X with both G_i invertible, the certified map is
A=G_2 G_1^(-1), b=M_2−A M_1. It carries the full driver law, not merely moments.
If covariances are nonsingular, whitening gives residual
O=Sigma_2^(-1/2) A Sigma_1^(1/2), with OO^T=I. Thus the frozen BL quotient distance
is **exactly zero** and the family admits one orbit reference. A singular target
can satisfy a pathwise identity while failing GL(k); it must not get this result.

`hb_oracle.py` uses independently written exact determinant and finite-support
arithmetic. Controls check support coverage, matching dimensions, positive measure,
nonsingular covariance, common representation changes and degenerate maps. The
1D witness gamma1^2=mu3^2/mu2^3 is affine/reflection invariant; unequal values
prove non-equivalence, while equal moments alone never certify equivalence.
An exact perturbation of 10^(-40) remains nonzero, not a numerical zero.

The proposed static-tail control uses **one common** globally invertible cubic
s(x)=x+alpha x^3, alpha>=0, on both protocols. Its derivative is >=1 and its range
is all R. Independently computed shape change disappears after exact removal of
the common calibration. This repairs the original protocol-specific/quadratic
fixture. These are the kit's partial HB-3/HB-4 controls, not clearance of the
unbuilt HB library or of a physical observer/carrier.

## 4. NR-4 regression and conservative semantics

The pre-fix source is pinned to
`f6ba47f79eba8ef594c92a618c8e07983a2a665c`, blob
`2a0361544d0ae95ec212e1a5399a6ba825603a26`. A coupling-insensitive empty variant
excludes its own entire denominator and falsely reports NOT-RELOCATED. The
CR-4 source and certified CR-5 boundary instead return **RELOCATED** on the same
abstract operational fixture. This is a required regression, not a new law test.

The independent oracle computes the full-K corner mask once. Instrumented subject
calls confirm exactly one certificate call per point and the **identical mask
object** in K, empty, standard and deletion runs. Missing ablated images stay in
the denominator without a hit. Zero-mass fibers produce no 0/0 verdict. Excluding
all positive mass leaves fraction None / VOID responsibility. The 183-case
regression also retains unresolved-reference rejection, enclosure boundaries,
representation transports and the distinct NR-17(b)/Q3(c) fallbacks.

All receipts here are explicitly **synthetic arithmetic fixtures**. No global
physical limit theorem, approved physical registry or experimental lock is
established by passing them.

## 5. CR-5 reference-measure boundary

`weight_boundary.py` materializes a copy of every point and converts the complete
weight sequence to exact Fractions. It validates finiteness and nonnegativity,
then strictly positive total mass, **before delegating to any supplied appearance,
corner, certificate or NR-4 arithmetic/callback**. Validation's own conversion
and total-mass sum are the validation step. Invalid lists give deterministic
`ValueError('REFERENCE_WEIGHT_CONTRACT')` with zero callbacks. Zero atoms in a
positive measure are allowed; an all-excluded positive measure stays VOID.

Sixteen invalid lists are checked at both appearance and NR-4 entrypoints,
including the original nine, negative early/late weights, mixtures yielding 2 or
−1, zero/empty totals, NaN/±Inf early/late, and an invalid late value after the
original implementation has visibly executed callbacks. All 32 boundary checks
pass. All 183 original valid results are byte-equivalent as comparison records,
and 240 extra exact weighted comparisons stay within [0,1] or return VOID.

The sealed private kit is **unchanged**. This package isolates the defect; it does
not silently fix or clear direct module calls. The mandatory boundary's source
hash is in the reviewer manifest. Canonical integration/re-freeze and actual
independent approval remain separate actions after this review.

## 6. Valid-input inventory and reproducibility

`evidence/TEST_INVENTORY.json` lists all **2,084 unique** replay entries: 24 frozen
oracle methods, seven document methods, 1,573 pricing comparisons, 14 HB comparisons,
11 NR-4 checks, 32 invalid entrypoint checks, 183 original comparisons and 240
additional measures. Method bodies disclose their exact finite subgrids. Separate
development integrity tests pass **65/65**, including corrupt inputs, frozen-oracle
tampering and a failed replay leaving no READY result. The default/overlay record
inventories and source outcomes are retained separately.

One command from the pinned repository checkout, with **Python 3.12.14** available:

```bash
bash development/final_precard/review.sh /absolute/path/GRUT_FINAL_PRECARD_REVIEW_INPUTS.zip /fresh/output/directory
```

The driver needs only Python's standard library and Git. It verifies the private
archive/member hashes, public review envelope, frozen oracle, document scope and
controls before importing subject code. No latest ref is fetched. A mismatched
input, missing evidence or infrastructure failure emits BLOCKED with nonzero exit;
a READY result is written only after inventory/provenance validation succeeds.
Fresh output directories are mandatory. Subject code is loaded in a disposable
directory, not copied into TestingGRUT.

## 7. Exact semantic differences

All supported-profile comparisons now agree. `SEMANTIC_DIFF.json` and the full
machine rows preserve every before/after result. The original coder disagrees on
619 repeated fixtures, chiefly boundary arithmetic and malformed-input acceptance;
that is an affected-fixture count, not 619 distinct bugs. The HB original
disagrees on seven fixtures.

The fresh comparison also found a real remaining defect in the prior isolated
proposal: alpha-equivalent lambda terms could price **106 vs 114** in normal form.
The repaired proposal standardizes every binding scope after NNF duplication.
The three p=6/10/16 discrepancies are preserved. Two additional underpricing
fixtures expose missing explicit prior-KL and table-entry lower-bound hooks;
those hooks now reject underpriced rows. No scientific threshold was changed.
The full source diff is in the private input archive; only its hash and the
independent test/contract recommendations are published here.

## 8. Pins

`REVIEW_MANIFEST.json` is the complete reviewer envelope: approved parent,
reviewed work-branch source SHA, public source/oracle/test hashes, all original/
proposed/private member hashes, canonical charter/STATE/CHECKS blobs, historical
pre-fix pin, G2-14R ruling hash, Python/Git versions, inventory and result digests.
The envelope is added after its source commit to avoid a self-referential hash;
the later envelope commit changes no reviewed executable bytes.

- Private inputs SHA256: `dc7e4437612c45b023b391350040f4557bec51c1369f3fab85ed784c1da7d856`
- Original CR-5 bundle SHA256: `f0de023add045314620ee21dddba760de5c3b5e350b4117cb03a46bf2b6eed4b`
- Issued G2-14R ruling SHA256: `d3b5f160735a71ce1017e5e01755f49b548aa9070b534d41c6e54283e6b087a0`
- Canonical STATE blob: `14c41d4487d74b8af86a7fc3cc273e23f7512789`
- Canonical CHECKS blob: `950d4e401bab9618d87cc71be2d3a9eed2a41031`
- Deterministic results SHA256: `76187c6885e13111137331263523ea158f481acb387971c7787fc84846e2bf9d`

## 9. External review boundary

The external checker must identify themselves and the exact reviewed source SHA,
input archive digest, mandatory isolation boundary and reproduction digest. They
must assess the source/proofs and unsupported-profile boundaries, not merely
rerun the builder's command. An unresolved discrepancy blocks approval. Work's
separate implementation, a connected-account rerun or a green hosted public
integrity job is not external authorship. Hosted CI has no private subject archive
and checks only public pins/controls plus the distinct engineering record.

## 10. Status

**READY FOR EXTERNAL CHECK.** No external signer exists in this package. Canonical
Stage-3 files, the frozen charter and the 25 bank dispositions are unchanged.
Card 1 remains CLOSED. No foundations lane or experimental-lock prerequisite was
opened. Actual external clearance of the pinned evaluator/CR-5 implementation
must precede the separate owner Card-1 decision.
