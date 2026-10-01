# SCOUT_0 CANDIDATE MATRIX — triage of the 24 routes

Scoring: 0–5 per dimension. **Assumption cost** and **Import risk**: higher = worse (so a *good* score on
those two is LOW). Everything else: higher = better. Scored from the BASELINE_MAP + RECON_MAP only; no
physics performed yet.

| # | ID | Leverage | InfoGain | Falsifiability | Distinctiveness | AssumpCost↓ | Tractability(maths) | Tractability(comp) | Empirical link | DependencyReach | ImportRisk↓ | Total* | Note |
|---|----|---------|---------|----------------|-----------------|-------------|--------------------|--------------------|----------------|-----------------|-------------|--------|------|
| 1 | R-01 | 3 | 4 | 4 | 3 | 1 | 3 | 4 | 0 | 4 | 4 | — | Cheap and clean; tests whether floor edges are coupling-family-specific |
| 2 | R-02 | 4 | 4 | 4 | 4 | 2 | 3 | 3 | 0 | 5 | 4 | — | Directly attacks the SFG-0 "does not migrate" wall with a different statistics |
| 3 | R-03 | 3 | 3 | 4 | 2 | 1 | 4 | 3 | 0 | 3 | 5 | — | High tractability; result either way is small but sharp |
| 4 | R-04 | 5 | 5 | 4 | 4 | 3 | 2 | 2 | 0 | 5 | 3 | — | Attacks the biggest open seam (noise origin) with a genuinely new class |
| 5 | R-05 | 4 | 5 | 5 | 4 | 2 | 4 | 3 | 0 | 5 | 4 | — | Clean analytic; different admitted-limit family than S5-1 |
| 6 | R-06 | 4 | 5 | 5 | 4 | 1 | 4 | 3 | 0 | 5 | 5 | — | Locates the CM/non-CM boundary; directly serves the unselected S-7 |
| 7 | R-07 | 3 | 3 | 4 | 2 | 2 | 5 | 4 | 0 | 3 | 4 | — | Cheap but low novelty (already a preserved option, partially known) |
| 8 | R-08 | 4 | 4 | 5 | 3 | 1 | 4 | 4 | 0 | 5 | 4 | — | The owner's own preserved option; audit-grade |
| 9 | R-09 | 4 | 4 | 5 | 4 | 2 | 3 | 3 | 0 | 5 | 4 | — | Genuinely unpriced selector; highest structural upside in C |
| 10 | R-10 | 4 | 4 | 5 | 4 | 2 | 3 | 3 | 0 | 4 | 4 | — | Same family as R-09; different physical principle |
| 11 | R-11 | 3 | 3 | 4 | 4 | 3 | 3 | 2 | 0 | 4 | 4 | — | Formally heavier; payoff contingent |
| 12 | R-12 | 3 | 4 | 4 | 3 | 2 | 3 | 3 | 0 | 4 | 4 | — | Interesting but risky to formulate |
| 13 | R-13 | 5 | 5 | 3 | 5 | 5 | 2 | 2 | 0 | 5 | 2 | — | Highest stakes, highest assumption cost, highest import risk |
| 14 | R-14 | 3 | 3 | 3 | 3 | 2 | 4 | 3 | 0 | 3 | 4 | — | Conceptual groundwork; possibly representational-only |
| 15 | R-15 | 5 | 5 | 4 | 5 | 3 | 3 | 3 | 0 | 5 | 3 | — | The frontier terminal says CLOSED at tested class; a *structural* change is legal |
| 16 | R-16 | 3 | 4 | 4 | 3 | 2 | 4 | 3 | 0 | 3 | 4 | — | Connects X-3; moderate |
| 17 | R-17 | 4 | 5 | 5 | 3 | 3 | 3 | 3 | 0 | 5 | 3 | — | Textbook bath; high tractability; risk of just reproducing standard QBM |
| 18 | R-18 | 3 | 4 | 4 | 3 | 2 | 3 | 3 | 0 | 4 | 4 | — | Clean; but the martingale property may fail for a reason worth knowing |
| 19 | R-19 | 4 | 4 | 4 | 3 | 2 | 3 | 2 | 0 | 5 | 3 | — | Weinberg's theorem is *known*; the GRUT-specific work is the classification |
| 20 | R-20 | 3 | 4 | 4 | 3 | 2 | 4 | 3 | 0 | 4 | 4 | — | Directly targets the SETTLED-NEGATIVE bridge's c₀-freeness |
| 21 | R-21 | 3 | 4 | 5 | 3 | 2 | 4 | 3 | 0 | 3 | 4 | — | Tests whether the S6 arrow is quantum-structural |
| 22 | R-22 | 4 | 4 | 4 | 4 | 3 | 3 | 3 | 0 | 4 | 3 | — | Tests whether criticality is a native formation variable |
| 23 | R-23 | 4 | 4 | 4 | 4 | 2 | 3 | 3 | 5 | 3 | 3 | — | Only route with an empirical link; fights INVISIBLE-BY-SUPPRESSION |
| 24 | R-24 | 3 | 4 | 4 | 3 | 3 | 3 | 3 | 4 | 3 | 3 | — | Empirical; kernel transport is earned so the input cost is low |

\* Total shown as "—" deliberately: **selection is by Pareto frontier, not by sum** (master instruction §6).

## Pareto frontier

Dominating dimensions vary. Reading the matrix, four clusters emerge on the frontier:

**Cluster F (the access/lift seam — R-08, R-09, R-10):** high falsifiability, high dependency reach,
moderate assumption cost. These attack the seam the canonical Q4 named, using constraints that were
*never priced*. Cheap to run. **All three go.**

**Cluster G (the generator/CM boundary — R-05, R-06, R-07):** R-06 has the best combination of
information gain, falsifiability, and low assumption cost in the whole matrix (it is the only route
that directly serves the unselected S-7 *and* produces a boundary theorem rather than a single
yes/no). **R-06 and R-05 go. R-07 is dominated** (everything R-07 gives, R-05/06 give with less
novelty risk).

**Cluster N (the noise parent — R-04, R-17, R-18):** R-04 has the highest leverage but high
assumption cost and low tractability; R-17 is the tractable version; R-18 is the clean-limit version.
**R-17 and R-18 go** (the tractable/clean pair), **R-04 is parked** — it cannot be executed honestly
without first knowing what R-17/18 say.

**Cluster E (empirical — R-23, R-24):** R-23 has the only genuine empirical link in the matrix and
directly attacks an INVISIBLE-BY-SUPPRESSION entry from the *rate-scaling* angle. **R-23 goes.**
R-24 is dominated by R-23 for this campaign (kernel-transport computability is real but the
finite-T interferometry observable is more model-dependent).

**Cluster S (substrate — R-01, R-02, R-03):** R-02 is the most interesting (attacks the SFG-0 wall
with different statistics); R-01 is a solid cheap probe; R-03 is dominated by R-01. **R-01 and R-02
go.**

**Foundational cluster Q (quantum structure — R-13, R-15, R-19):** R-13 has the highest leverage
*and* the highest assumption cost; it cannot be honestly executed in a single campaign. **R-13 is
parked with a note.** R-15 is the strongest structural-change-legal route into the closed
Experiment-P frontier. **R-15 goes.** R-19 is a classification exercise with known literature;
deferred.

## Selected probes (9; balanced per the master instruction)

**Foundational (3):** R-06 (CM boundary of the memory kernel), R-08 (faithful-representation audit),
R-09 (locality-of-representation constraint).
**Constructive / model-building (3):** R-01 (floor-edge portability across coupling families), R-02
(boson-statistics sector-formation), R-17 (caldeira-leggett bath as C-B candidate).
**Empirical / predictive (2):** R-23 (decoherence-rate scaling observable), R-15 (nonlinear
branching → Born-weights question, analytic formulation + a small numeric check).
**Clean-limit complementary (1):** R-18 (ergodic fast-variable averaging; runs as the null to R-17).
**Hostile / null route (folded into the hostile phase, §9 of the instruction):** R-03 serves as the
null discipline — the coarsening-stability test applied to any survivor that claims structural
status.

Parked, with reasons recorded so the next scout does not re-derive the decision: R-03 (dominated),
R-04 (needs R-17/18 first), R-07 (dominated), R-11, R-12 (contingent), R-13 (assumption cost too
high for an honest single pass; the scout notes it as the single most important *future* route),
R-14 (possibly representational-only), R-16 (moderate), R-19, R-20, R-21, R-22 (deferred to the
second-wave list), R-24 (dominated by R-23).

## Execution order (dependency-driven)

R-06 → R-08 → R-09 (foundational, fast, sharpen the language)
  → R-17 → R-18 (bath pair; R-18 runs as R-17's null)
  → R-01 → R-02 (substrate portability)
  → R-23 (empirical, uses whatever survives from R-06/R-17)
  → R-15 (analytic formulation; numeric check only if time allows)
  → hostile pass on the survivors.
