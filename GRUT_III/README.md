# GRUT III — the responsive vacuum, principle-first

**Status:** PROPOSAL on the experimental branch (`claude/friendly-hopper-7zkwf5`, branched from the SCOUT_0
freeze `ab2da47`). Nothing here alters the canonical record; every claim below is tiered and every
computation is runnable. Date: 2026-10-11.

## The statement in one paragraph

GRUT I/II asked a *mechanism* ("the vacuum is a responsive medium with finite memory") to select
physics, and its own audits proved — honestly and exhaustively — that a mechanism without a principle can
only classify: at every interface the observable factors through a quotient object whose value the theory
supplies (SCOUT_0 handoff §4). GRUT III replaces the mechanism by a **principle of identification**:
*the bath of any observer is the algebra of their causal complement.* The system/bath split — the record's
deepest supplied item — becomes the causal structure itself. Three theorems then do the work that was
previously assumed: the medium's state is KMS at the horizon temperature (Bisognano–Wichmann, Gibbons–Hawking),
its noise and memory kernels are fixed in closed form (derived here, §D2), and its equilibrium equation of
state is Einstein's equation with Λ as an integration constant (Jacobson), bought with the two inputs the
record was already paying for. In the matter sector the dissipative kernel is derived from the optical theorem:
Im Π(s) = c_eff(s) G s², with the "spectral edge data" (edge locations 2mᵢ, exponents 5/2, 3/2, 1/2, weights
1 : 3 : 12) fixed by the particle spectrum. The responsive-vacuum program, done right, is GR + quantum field
theory in curved spacetime with horizon thermodynamics as the organizing principle: the kernel is nowhere
supplied, GR is hosted, the record's frontier question ("what bath Hilbert space was integrated out") is
dissolved rather than solved, and the supplied ledger falls from +16 to +5. The price is candour: the
single-pole finite-memory vacuum of GRUT I is gone (the local response is instantaneous; the memory lives in
the noise and in the threshold structure), and the sharp, parameter-free predictions at accessible scales are
those of GR + QFT — stated here as falsifiers, not as novelties. The genuinely open items are the ones every
physicist knows: the area-law coefficient (G), the value of Λ, the outcome rule, the arrow, the UV completion.
One bold conjecture (a₀ ∼ cH₀) is carried, fenced, with its discriminating test.

## Files

| file | content |
|---|---|
| `00_LETTER.md` | Why the program stalled and what the principle is (the argument, in prose) |
| `01_PRINCIPLES.md` | Postulates P0–P5, data, tiers; what each replaces in the old ledger |
| `02_DERIVATIONS.md` | D1–D8 with equations; every equation points at a verification script |
| `03_PREDICTIONS_AND_FALSIFIERS.md` | What GRUT III predicts with zero free parameters, what kills it, what it retracts |
| `04_LEDGER_MIGRATION.md` | Old register node → new status; the +16 → +5 accounting |
| `05_ROAD_TO_COMPLETION.md` | The open primitives, posed sharply, with the next computation for each |
| `06_CONJECTURE_A0.md` | The one fenced conjecture (a₀ ∼ cH₀) with its discriminating test |
| `RESULTS.md` | Output of the verification suite (generated) |
| `calc/` | `v1`–`v6` + `run_all.py` (Python 3.13, numpy, sympy, mpmath) |

## Run

```
pip install sympy mpmath        # numpy is present; stdlib otherwise
python3 GRUT_III/calc/run_all.py   # writes GRUT_III/RESULTS.md; exit 0 iff every check passes
```

## Discipline

Tiers used: **THEOREM** (standard mathematics/physics, proof or citation given, verified here where
feasible), **DERIVED** (follows from the postulates by a computation in `calc/`), **POSTULATE** (assumed;
counted in the ledger), **DATUM** (a measured constant every theory must take; not counted, per the
charter's Λ ruling), **CONJECTURE** (not derived; carried with its falsifier). Every "new in GRUT"
statement is a *placement* of known physics unless explicitly marked otherwise; two closed forms in §D2
and the species-resolved kernel in §D5 are, to this author's knowledge, not on the GRUT record.
