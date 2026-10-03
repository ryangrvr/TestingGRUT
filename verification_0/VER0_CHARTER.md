# VER0 — INDEPENDENT VERIFICATION 0: CHARTER

**Branch:** `grut-independent-verification-0`, created from the frozen scientific parent
**`scout-0 @ ab2da47407bb670d94e2a52c87599fa13fd8ab99`** (unmodified).
- Not branched from BRI0 or DA0.
- Nothing here modifies `scout-0`, `grut-backreaction-identifiability-0`, `grut-directed-autonomy-0` or
  `grut-reflexive-autonomy-0-frozen`.
- No PR, no merge.

**Purpose (verification / reproduction, not exploration):**
1. **VER0-A:** discharge SCOUT-0 criterion-2 reproduction debts.
2. **VER0-B (registered, not opened):** independently reconstruct the closed BRI1-X1 theorem from a sealed,
   statement-only target.

## 1. Provenance

The frozen SCOUT-0 handoff §7 ("Packaging/verification work, not exploration") lists the owed reproductions:

| item | owed reproduction |
|---|---|
| 1 | P-17 identifiability theorem and its exact Taylor checks |
| 2 | P-02 edge-velocity theorem, plus the hard-core-boson ↔ SF-1 match |
| 3 | P-15 Born ⟺ martingale ⟺ decomposition-independence chain |
| 4 | K1-H / K1-HS reproduction, with primary texts read in full (Linder 2020; Peirone et al. 2018; Bellini–Sawicki 2014), including the strong-branch series and the shooting family |
| 5 | EDA-01 row classification checked by a second reader |
| 6 | Secondary primary-text checks: the small-α Horndeski ratio; the Pogosian–Silvestri conjecture; the P-18 homogenisation citations |

**SCOUT-0 remains PROVISIONALLY SATURATED — OWED CHECKS ONLY** until the required grades are reached. VER0 keeps its
own ledger and never edits the frozen SCOUT-0 record.

## 2. Independence grades

**VER-I1 — independent derivation / code path.** The result is reconstructed:
- without copying its proof or implementation;
- without importing intermediate symbolic expressions, beyond the target statement, its definitions and its declared
  assumptions;
- with a fresh derivation and, where computational, a freshly structured implementation.

VER-I1 is the **maximum grade the present automated workflow can earn**. It is **not** "independent external review" and
**not** "independent scientific replication by another researcher".

**VER-I2 — external independent review.** A genuinely separate human or independently operating reviewer. **Not
earnable automatically.** Reserved for an actually supplied review.

**Exposure qualifier (mandatory where it applies).** If the agent that orchestrates a reproduction has prior exposure to
the original derivation, the result is graded **VER-I1 (orchestrator-exposed)**. The reproduction itself must then be
done by a **context-isolated sub-agent** that receives only the frozen target spec. Its file accesses are reported and
audited.

## 3. Anti-leak firewall

For every reproduction:
1. **Target spec** (`specs/…_TARGET.md`). It contains only:
   - the result statement;
   - assumptions;
   - the definitions needed to formulate it;
   - sources;
   - acceptance tests.

   The original result files may be read only to extract the spec. The source file's blob hash is recorded. **The spec
   is committed before reproduction begins.**
2. **Sealing.** After the spec is frozen, the original derivation is **sealed** until the reproduction result is
   committed.
3. **Reproduction.** A fresh derivation in `V0_k_…_REPRODUCTION.md`, and fresh code (only if needed) under
   `code/v0_k_*`.
   - No original source code is copied into `verification_0/`.
   - No target tests are copied blindly.
4. **Comparison phase**, only after the reproduction is committed. Unseal the original and classify:
   - exact agreement;
   - equivalent derivation;
   - correction;
   - discrepancy;
   - unresolved.

   The independent derivation is never edited to resemble the original. Corrections are preserved for owner review.

## 4. BRI1 firewall (registered now; not run)

**Closed target branch:** `grut-backreaction-identifiability-0 @ f2e6999c71b8b621868a019299085f64b80dce46`.

**In VER0-A, its BRI1 proof and code are not inspected**: the fact that BRI1 exists is metadata only. VER0-B may open
after VER0-A, or by owner authorisation.

**VER0-B target:** a separate **statement-only** spec. It **contains:**
- the X1 Hamiltonian;
- the preparation;
- the clamp protocol needed;
- the E₂± definition;
- the theorem statement;
- the allowed standard tools.

It **excludes:**
- the BRI1-T1 and BRI1-R1 proofs;
- the Taylor-coefficient derivation;
- the symbolic scripts;
- the PF4Q values;
- the quadrature code;
- the original proof structure.

**Disclosure (recorded now).** The agent operating VER0 is the author of BRI1, which is maximal exposure. **VER0-B must
therefore be carried out by a context-isolated sub-agent (or a human), never by the orchestrator itself, and is graded
VER-I1 (orchestrator-exposed) at most.**

## 5. VER0-A order (frozen; not reordered after results)

| step | target |
|---|---|
| **V0-1** | P-17 |
| **V0-2** | P-15 |
| **V0-3** | P-02 / hard-core-boson ↔ SF-1 |
| **V0-4** | EDA-01 second-reader audit |
| **V0-5** | K1-H / K1-HS |
| **V0-6** | secondary primary-text checks |

## 6. Primary-text policy

For every literature-dependent check, record:
- exact title, authors, publication, year;
- DOI / arXiv identifier;
- **whether the full text or only an abstract / snippet was inspected**;
- the exact claim checked.

**Not acceptable:** snippets, secondary summaries, or GRUT documents quoting the source. If the primary text cannot be
obtained, the grade is **ACCESS-BLOCKED**, never inferred.

## 7. Per-item grades

| grade | meaning |
|---|---|
| V0-k-A | REPRODUCED, VER-I1 (with the exposure qualifier if applicable) |
| V0-k-C | REPRODUCED WITH CORRECTION |
| V0-k-F | REPRODUCTION FAILED |
| V0-k-I | INDETERMINATE (precise obstruction) |
| ACCESS-BLOCKED | for primary-text items |

**Criterion 2 is not declared met** until every listed item reaches its required grade. The I1 / I2 distinction is kept
visible throughout.

## 8. Scope

No new physics, no exploration, no new campaign.
