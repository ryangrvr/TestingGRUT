# EC01_ITER2_EDGE_AUDIT_PROVENANCE — Resolution of the provenance review

> Provenance-only artifact. No scientific content. Closes the open provenance
> question left at Iteration 3. Canonical GRUT-RAI untouched.

## Classification: OMITTED_REQUIRED_ARTIFACT (now archived)

The file `program/EC01_ITER2_EDGE_AUDIT.json` was determined to be a genuine
contemporaneous Iteration-2 artifact that was inadvertently omitted from commit
`03a18e9`.

## Evidence

1. **Archival commit exists.** Commit `dd47a92` ("E&C-01: archive omitted
   Iteration 2 edge-audit artifact (provenance only; content verified against
   03a18e9; no classifications altered)") added the file with exactly that
   provenance-only purpose.
2. **Semantic verification.** Independent re-verification in this run: the
   working-tree copy and the committed copy at HEAD are *semantically equal*
   when parsed as JSON (`semantic_equal=True` under full parsed-object
   comparison, not textual comparison).
3. **No scientific alteration.** The archival commit message and content
   inspection confirm no classification, edge, or source reference was changed
   relative to the Iteration-2 audit conclusions recorded in `03a18e9`.
4. **Traceability.** The file is consistent with the edge identifiers and
   classifications referenced by `EC01_ITER2_HOSTILE_AUDIT.*` (committed at
   `03a18e9`).

## Disposition

- The provenance question is **CLOSED**.
- The archival commit `dd47a92` is judged adequate; no further archival
  commit is required.
- The file does **not** alter any E&C-01 scientific conclusion.

## Residual anomaly discovered during this run (reported, not repaired)

An uncommitted working-tree diff exists across 34 tracked JSON artifacts
(~3,400 insertions), including historical frozen E&C and Reality-Check
artifacts. Characterization performed before any work:

- **30 of 34 files**: semantically identical to HEAD (whitespace/formatting
  only under parsed comparison).
- **5 files** (`REALITY_COMPARATIVE_FOUNDATIONS_01_CASE_02/03/04/05/07.json`):
  both the committed copy AND the working copy fail JSON parsing
  (malformed, hand-authored JSON-like text); they additionally show real
  non-whitespace differences (~145 insertions).
- Nothing is staged (`git diff --cached --name-only` empty).
- Protocol hash verified unchanged:
  `6b012117…27a00a8` (full value recorded in the freeze report).

This anomaly is OUT OF SCOPE for DP-R0-01. It is recorded here and in the
freeze report so it cannot be silently laundered into any future commit.
Recommendation for a later, separately authorized run: forensically identify
the origin of the working-tree rewrite before any further commit; treat the
5 malformed files as requiring a dedicated validity repair pass with its own
audit trail.
