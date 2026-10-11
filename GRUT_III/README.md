# GRUT III — the responsive vacuum, principle-first (revision 2.4)

**Status:** PROPOSAL on the experimental branch (`claude/friendly-hopper-7zkwf5`, branched from the SCOUT_0
freeze `ab2da47`). Revision 2 (2026-10-11) incorporates three external reviews and an owner's note:
`RESPONSE_TO_REVIEWS.md` is the point-by-point record, `CHANGELOG.md` the diff. Nothing here alters the canonical
record; every claim is tiered and every computation is runnable.

## The statement in one paragraph

GRUT I/II asked a *mechanism* ("the vacuum is a responsive medium with finite memory") to select physics, and its own
audits proved that a mechanism without a principle can only classify. GRUT III replaces the mechanism by an
identification stated algebraically: *the medium an observer sees is the triple (local algebra, restricted state,
modular flow); the bath is the commutant.* Imported theorems then supply what the record had assumed — horizon
thermality (Bisognano–Wichmann, Kay–Wald, Gibbons–Hawking), Einstein's equation as an equation of state with Λ an
integration constant (Jacobson; G inserted as a datum), the horizon's transport coefficients — and computations supply
the kernels: the de Sitter stress-tensor noise spectrum in closed form (response local, memory in the noise), the
matter-sector absorptive kernel Im Π(s) = c_eff(s) G s² from explicit amplitudes, the gravitational sector as a
stationary thermal process plus a Brownian zero mode with D = H³/4π². No kernel is supplied anywhere, and nothing is
earned by that alone: every kernel is standard physics. The reviews established the structural fact that organises
revision 2: because GRUT III builds the KMS state in, the fluctuation–dissipation identity N = ½coth(βω/2)ρ forces any
"extra noise" to be extra dissipation, so **in its own regime GRUT III is standard by construction**. A distinctive
GRUT statement must therefore be a *residual* Δρ (or a two-time ΔW, or a cumulant Δκₙ) for a named observable and
state, computed where KMS fails; the conventional O(ε) baseline under slow roll is computed here (`v11`), and the only
open sector is the interacting gravitational infrared in (quasi-)de Sitter. GRUT III has five framework postulates, a
list of data, a short list of open problems posed sharply (G's coefficient — Sakharov's route fails with Standard-Model
content; Λ; the outcome rule; the arrow), and no distinctive prediction. It is a correct place to stand, not a theory of
reality.

## Files

| file | content |
|---|---|
| `00_LETTER.md` | Why the program stalled and what the principle is (the argument, in prose) |
| `01_PRINCIPLES.md` | Postulates P0–P5, data, tiers; what each replaces in the old ledger |
| `02_DERIVATIONS.md` | D1–D8 with equations; every equation points at a verification script |
| `03_PREDICTIONS_AND_FALSIFIERS.md` | What GRUT III predicts with zero free parameters, what kills it, what it retracts |
| `04_LEDGER_MIGRATION.md` | Old register node → new status; the +16 → +5 accounting |
| `05_ROAD_TO_COMPLETION.md` | The open primitives, posed sharply, with the next computation for each |
| `06_CONJECTURE_A0.md` | Appendix: the a₀ ∼ cH₀ coincidence (Milgrom's conjecture, not GRUT's) |
| `07_ACCESS_RESOLVED_RESPONSE.md` | The Δρ residual framework: FDT identity, trichotomy, memory dictionary, what a GRUT claim must be |
| `08_PRIORITY_AND_SOURCES.md` | Literature priority checks (what PNRV 2009, Cho–Hu 2014, MGM 2019, DESI DR2, HLZ 1999 … actually contain) |
| `09_TWO_TEMPERATURES.md` | Detector vs apparent-horizon temperature at O(ε): which one P1 calls a temperature; massless sector done (`v14`) |
| `10_SECULAR_IR_SURVEY.md` | Literature check for O5 item 3: do the secular-IR resummations disagree? (scalars no; gravitons: a QG observables dispute); the `v13` result |
| `11_VERIFICATION_ROUND.md` | Sixth-round independent re-derivation and adversarial audit of `v8`, `v11`, `v12`, `v13` (workflow report) |
| `RESPONSE_TO_REVIEWS.md` | Point-by-point response to the reviews (R1–R20) |
| `CHANGELOG.md` | Revision history |
| `RESULTS.md` | Output of the verification suite (generated) |
| `calc/` | `v1`–`v14` + `run_all.py` + `requirements.txt` (Python 3.13, numpy, sympy, mpmath) |

## Run

```
pip install -r GRUT_III/calc/requirements.txt   # numpy, sympy, mpmath
python3 GRUT_III/calc/run_all.py   # writes GRUT_III/RESULTS.md; exit 0 iff every check passes
```

## Discipline

Tiers used: **THEOREM** (standard mathematics/physics, proof or citation given, verified here where
feasible), **DERIVED** (follows from the postulates by a computation in `calc/`), **POSTULATE** (assumed;
counted in the ledger), **DATUM** (a measured constant every theory must take; not counted, per the
charter's Λ ruling), **CONJECTURE** (not derived; carried with its falsifier). Every result here is a *placement* of standard physics on the GRUT record; imported theorems earn GRUT nothing and
are labelled IMPORTED. The closed forms of D2, D9, D10 and the species kernel of D5 are corollaries of standard
results (priority: `08`), verified here by computation.
