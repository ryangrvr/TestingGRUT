# E&C-BORN-01 — Evidence Assembly Report (Phase 1: Source Collection)

> Source-level evidence package only. NO route verdicts, NO premise adjudication,
> NO circularity classifications, NO scores or rankings were assigned.
> Canonical GRUT-RAI untouched. TestingGRUT laboratory record.

## Scope

Phase 1 of E&C-BORN-01: collect and verify primary derivation sources and
identified critiques for the eight frozen route families (R1–R8), with
statement-level records, target-claim mapping, and descriptive provenance.

Per protocol: critique/adjudication of circularity, independence, and
legitimacy are explicitly deferred to a later authorized phase.

## Source coverage by route

| Route | Family | Primary source | Status | Critique identified | Status |
|---|---|---|---|---|---|
| R1 | Gleason / Busch / noncontextual measure | Gleason 1957 (J. Math. Mech. 6, 885) | PARTIALLY_VERIFIED (theorem content canonical; verbatim pending) | — | — |
| R1 | (cont.) | Busch 2007 (arXiv:0712.1207) | UNVERIFIED, pending | — | — |
| R2 | Envariance / symmetry | Zurek 2003 (PRL 90, 120404) | **ABSTRACT VERIFIED** — verbatim quote recorded | Hemmo–Pitowsky 2007 (arXiv:quant-ph/0703184) | UNVERIFIED, pending |
| R3 | Deutsch–Wallace decision theory | Wallace 2009 (arXiv:0906.2718) | **ABSTRACT VERIFIED** — verbatim quote recorded | Kent 2009 (arXiv:0905.0624) | **METADATA VERIFIED** via arXiv API: "One world versus many: the inadequacy of Everettian accounts of evolution, probability, and scientific confirmation" (OUP *Many Worlds?* 2010) |
| R4 | Typicality / frequency / self-locating uncertainty | Sebens–Carroll (arXiv:1405.7577v3; BJPS 69, 25–74, 2018) | **METADATA VERIFIED** via arXiv API (was GAP; now filled) | Adler 2003 (arXiv:quant-ph/0208169) | UNVERIFIED, pending |
| R5 | Quantum reconstruction / operational axioms | Masanes–Müller 2011 (arXiv:1004.1483) | UNVERIFIED, pending | — | — |
| R6 | Bohmian quantum equilibrium | Dürr–Goldstein–Zanghì 1992 (J. Stat. Phys. 67, 843) | PARTIALLY_VERIFIED (equivariance/typicality content canonical) | — | — |
| R7 | Bohmian nonequilibrium | Valentini 2001 (arXiv:quant-ph/0106098) | UNVERIFIED, pending | — | — |
| R7 | (cont.) | Valentini–Westman (arXiv:quant-ph/0403034; Proc. R. Soc. A 461, 253–272, 2005) | **METADATA VERIFIED** via arXiv API | — | — |
| R7 | (cont. — supporting) | Valentini (Phys. World 22N11, 32–37, 2009, "Beyond the Quantum") | VERIFIED (retrieved earlier in phase) | — | — |
| R8 | (cont. — supporting context) | Percival (arXiv:quant-ph/9508021, "Quantum space-time fluctuations and primary state diffusion") | **METADATA VERIFIED** via arXiv API | — | — |
| R8 | Collapse stochastic law | Pearle 1984 (Found. Phys. 14, 1139) | PARTIALLY_VERIFIED (postulated noise-measure content canonical) | — | — |

DP-R0-01 admitted as **CONTEXT only** per the frozen protocol; its verdict was
not imported into any route record.

## Representative statement-level evidence (verbatim, from fetched text)

### R2 — Zurek (envariance)
Source: S-R2-01, abstract, [arXiv:quant-ph/0211037](https://arxiv.org/abs/quant-ph/0211037):
> "I introduce environment-assisted invariance -- a symmetry related to causality that is exhibited by correlated quantum states -- and describe how it can be used to understand the nature of ignorance and, hence, the origin of probabilities in quantum physics."

Author-claimed target: WEIGHT_REPRESENTATION + OBJECTIVE_CHANCE.

### R3 — Wallace (Deutsch–Wallace)
Source: S-R3-01, abstract, [arXiv:0906.2718](https://arxiv.org/abs/0906.2718):
> "I develop the decision-theoretic approach to quantum probability, originally proposed by David Deutsch, into a mathematically rigorous proof of the Born rule in (Everett-interpreted) quantum mechanics. ... lastly consider a number of proposed ``counter-examples'' to show exactly which premises of the argument they violate."

Author-claimed target: WEIGHT_REPRESENTATION + RATIONAL_CREDENCE; the author
himself flags premise-dependence via the counter-example analysis.

## Assumption extraction (descriptive, pre-adjudication)

No assumption-independence or circularity classification was performed. The
ledger records author-characterized assumptions where sources state them and
marks what remains to be extracted from full texts.

## Access limitations recorded

The execution environment provided abstract-level access for arXiv items.
Statements quoted verbatim in this report are actually observed. All remaining
sources are recorded as UNVERIFIED/PARTIALLY_VERIFIED with their exact claims
pending full-text retrieval. This is explicit in the ledger.

## Known gaps (must be filled before adjudication)

1. ~~R4 requires a genuine primary derivation source~~ **RESOLVED: Sebens–Carroll
   (1405.7577v3, BJPS 69, 25–74) verified 2026-09-12 via arXiv API and registered
   as the R4 primary.**
2. R2/R3 critiques (Hemmo–Pitowsky, Kent) need full-text fetch (Kent metadata
   now VERIFIED).
3. R1/R5/R7 primaries need full-text fetch (R5 Hardy quant-ph/0101012v4 and
   R7 Valentini–Westman quant-ph/0403034 metadata VERIFIED via API).
4. No R1, R5, R6, R7, R8 critique sources collected yet.
5. Two initially wrong arXiv IDs were caught and replaced during the API
   verification pass; the intended Kent critique (0905.0624) is confirmed
   correct; the extra fetch returned Percival (9508021), retained as R8
   supporting context.

## Firewalls respected

- decoherence ≠ outcome selection ≠ Born weights — preserved throughout.
- Cbp boundary recorded as scope information only; no substantive application.
- No CIRCULARITY_STATUS assigned. No route verdicts. No rankings.
- No escape-door map modification. No existing artifact modified.
- Nothing staged, committed, or pushed.

## STOP

Execution of Phase 1 is complete. Adjudication phases require separate
authorization.
