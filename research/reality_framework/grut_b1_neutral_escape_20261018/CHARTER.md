# GRUT-B1: neutral-composite escape, frozen bounded cycle

Date: 2026-10-18 UTC. Parent: A1 `22916ef985231ddc4b0f5170abc64687cc1d1fe7`.
This is a separately named candidate, not an A1 variant or a rescue of its
failed extensive, every-time gate. Canonical records and all old gates stay fixed.

## Question and hypotheses

Can initially excited fuel discharge into a mobile neutral composite and
transmit an internal qubit before substantial return to fuel, with conversion
and its adjoint present throughout the experiment?

Two mechanisms are assessed, not two numerical searches:

1. A gap that merely detunes fuel from matter suppresses production as well
   as loss. Give the exact two-level obstruction; do not launch a variant.
2. The primary B1 construction couples fuel resonantly to an extended
   low-energy composite band. Spatial escape, rather than a removed reverse
   operator, can permit a finite operating interval. Test this one candidate.

Pre-execution analysis motivates the candidate and bounds below; this is not
a blind discovery claim. No primary matrix trajectories have run at this freeze.

## Fixed primary specification

Use P1's hard-core charges, positive capacity-one electric links, mass 1 per
constituent, density coefficient 1/4, and hoppings j+=1/64, j-=1/96. Replace
adjacent-charge conversion by *on-event neutral* conversion at a supplied
fuel site. Add an explicitly supplied degree counterterm. Include a two-state
internal label on the plus carrier and excited fuel. See `GRUT_CANDIDATE.md`.

Let A=13/46080, J=1/3840, Omega0=2-2A. Primary class:

    |Omega-Omega0| <= J/1000, |lambda-J| <= J/1000.

All other microscopic matter parameters and counterterms are fixed. This is
a narrow, size-independent interval, not broad perturbation universality.
The neutral fuel is initially excited; matter and fields initially are vacuum.
No clock, switch, cooling, reset, dissipative law, or favorable matter seed.

First analyze one resource on a fixed recorded path of at least 64 events.
Then prove the conditional full-geometry extension from an accepted record,
retaining its actual receiver state. The latter uses N>=2^28 and m>=N/2,
not a new claim that this record was obtained autonomously by B1.

## Frozen event and gates

Use the *same* Hamiltonian and units throughout; put theta=J t.

* For every theta in [2,4], matter probability must exceed 1/2.
* At theta=4, the joint event of preserved retained prefix, a co-located
  neutral composite at recorded distance >=4 from the source, and successful
  decoding of its internal qubit must have probability >=3/5, per prepared
  fuel trial conditional on the original geometry record. Failure is an
  erasure, not renormalized away. Prove arbitrary-qubit fidelity on success.
* Reverse conversion and every allowed charged separation remain present.
* Derive a time-domain error, not only fixed-energy Schur eigenvalues.

The primary theorem is one-resource. Two neighboring fuel cells on a common
short path are a required secondary exact diagnostic, not an extensive
theorem or another primary pass criterion. No parameters change after it.

## Analysis and controls

Use an explicit approximate intertwiner into the complete one-pair Hilbert
space, with a residual integrated by exact unitary Duhamel. Bound end effects
against the half-line walk. Check scalar load-bearing constants with rational
arithmetic; small matrix checks use ordinary double precision.

Freeze these diagnostics: one-resource paths 2, 8, 16; two fuel sites at the
ends of one 8-event path with exactly two resources; the five fixed theta
values 0,1,2,3,4. On the 16-event one-resource matrix also inspect 41 evenly
spaced points in theta in [2,4]. Two-color factorization, conservation,
Gauss/capacity, actual full-state intertwining residuals, and energy are checked.
No disorder samples or favorable-time search. These finite paths do not
replace the analytic >=64-event or full-geometry theorem.

Numerical ceiling: 60 seconds total across primary execution and one replay,
1 GiB process RSS, 100 MiB new files. Single CPU worker, no new packages,
allocation, persistent job, or paid resource. If a ceiling is reached, preserve
the stop and return unresolved; no fallback or revised budget.

## Decision and stopping

Stop after the analytic challenge and these controls. Select: controlled
conditional one-composite regime; failed mechanism; invalid/unproved bound;
or computationally unresolved. Regardless of a pass, do not claim autonomous
geometry, extensive production, solved collisions, selected dimension,
gravitational dynamics, independent validation, or evidence about nature.

Builder/challenger/integrator/recorder are functions of this one assistant,
sharing context and implementation. Their agreement is not external review.
Archive the five requested integrated documents, exact source pins, code,
raw numerical output, hashes and execution limits. Publish only on this
isolated research branch. End the cycle rather than opening a repair series.
