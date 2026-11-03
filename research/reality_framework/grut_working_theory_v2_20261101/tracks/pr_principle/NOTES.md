# Can a GRUT principle fix the holder? Result

Claude, for the owner (issue #13). The candidates and the rule were committed first, in
[PREREGISTRATION.md](PREREGISTRATION.md) (commit c697fae). Every choice remains a model choice. Nothing here validates
or excludes anything.

This is my own analysis on the existing solvers and scans; no separate skeptic was run, to save budget. New code:
- `code/edge.py`, which writes `results/edge.json`;
- `code/phase_soft.py`, which writes `results/phase_soft.json`.

Both reuse fx_holder's stationary solver read-only, through the e2e_chain skeptic's import path.

## Answer

**No candidate principle determines c0, δ or s_* without being tuned to the branch.** The only one that lands in the
windows is "minimise S/Q", and that is the branch criterion itself.

**The reheating closure is different: it is not a free parameter at all.** Fermi statistics decides it. It needs a
computation, a fermion field evolved on the lattice background, not a principle.

## Candidate by candidate

| | candidate | outcome | counts? |
|---|---|---|---|
| A | closure as defined | Constrains none of the four. c0 = \|g_h\|s_*/(12.30 H_i) is a coupling times a state datum; the closure pass showed closure fixes neither couplings nor the state [D, SYNTHESIS.md] | no value |
| B | maximal CP torque (held X–Σ relative phase = −90°) | Gives δ = −11.0 / −11.9 / −12.8 / −13.5 / −14.3° at c0 = 1.40 / 1.45 / 1.50 / 1.55 / 1.60 [C]. Misses the M2-Q window (−4…+4°) and the M2-X window (−2…+6°) by 7–10°. At δ = 0 the relative phase is −77° to −76°. Being conjugation-even, it fixes \|phase\| only, never the sign (preparation theorem) | miss |
| C | typicality of s_*: Starobinsky–Yokoyama equilibrium, a Hartle–Hawking-type ground state, the Boltzmann measure behind recurrence | The tree-level V_min(\|u(s)\|) dominates s's potential, with m²_eff ≈ −6×10⁻⁴ H² at M2-X. No equilibrium sits inside the window. With the required λ_s the attractor is at c ≈ −0.1 (holder nearly off), and P(window) ≲ exp(−6×10⁷) [C, e2e_chain skeptic sk_natural.json]. Typicality **anti-predicts** the window | anti-predicts |
| D | maximal yield (ε or \|N_X\|) | Both rise monotonically over the whole scanned range: ε = 0.08 → 0.79 from c0 = 1 to 1.60 (GM), 0.09 → 0.70 to 1.50 (P0) [C, e2e_chain scans]. There is no interior maximum. The held state stays stable to at least c0 = 2.6 (E), so the principle says "as large as possible" and gives no value | no value |
| E | edge of stability (a held mode turns tachyonic) | No edge for c0 = 1–2.6 at δ = −4, 0, +4°. The lightest non-gauge held mass² stays ≥ 8.75 H² and then grows (24.7 H² at 1.6, 93 H² at 2.6); the Y amplitude grows smoothly [C] | no value |
| F | extremal S/Q | Maximum S/Q points to a weak holder (c0 ≈ 1, S/Q = 15–61): free X, outside every window. Minimum S/Q gives the floor at c0 ≈ 1.50, δ ≈ +2.5°, inside the windows, but only because the windows are defined by S/Q | branch-defining by construction |
| G | reheating closure | Not a choice. The produced singlets have k ≈ m_χ/2 ≈ 125 H_i (wavelength ≈ 0.05 H_i⁻¹). The χ fragments at the gradient peak are ~2π/m_χ ≈ 0.025–0.035 H_i⁻¹ [E, from the lattice spectrum scale]. The fermions are therefore not localised inside one fragment, so the per-site closure P0 is not literally right, and the truth lies between P0 and GM [A]. It is decided by evolving a fermion field on the lattice's χ(x, t) | computation, not principle |
| H | observer selection | Both branches have CDM-like P(k) (to 10⁻⁵), the same Higgs safety and an independent η_B. Nothing observers need differs, so there is no selection [A] | no value |
| I | UV fixed point for g_h | Even a predicted g_h leaves s_* free; s_* is an initial condition with no equilibrium (C). So c0 stays undetermined | no value |

**Post hoc, not pre-registered, and not counted.** The softest held mode has its minimum (m² = 8.75 H²) at
c0 = 1.08 / 1.13 / 1.17 for δ = −4 / 0 / +4° [C]. That sits next to the c0 ≈ 1.10 at which M1r's own single-field
history turns Q-matter at e = 10⁻⁶. It is outside every M2 window (1.40–1.60), and M1r's single-field history is not
realisable (inflaton pass). It is recorded as a coincidence to keep in view, not a result.

**Chance-hit baseline.** Take a log-uniform prior on c0 over two decades. The windows then cover about 3% of it, and
about 5% of the phase range (half-circle). With five value-giving candidates (B, D, E, F, plus the post-hoc one), the
expected number of chance hits is about 0.15–0.25. Observed non-tuned hits: 0.

## Why, structurally

- **c0 is a coupling times an initial condition.** A principle that fixes it must also fix the pre-inflationary state
  of a light field. The only principles about states that we have are typicality and equilibrium measures, and they
  point the wrong way by an enormous factor (C).
- **The phase is CP-odd.** Any conjugation-even principle can at most fix |phase| at a CP-even extremum, and the
  natural one (maximal torque) misses by about 10°. A principle that fixed the sign would itself have to break CP. It
  would then just be another supplied CP phase.
- **The windows are defined by S/Q.** Any principle stated in terms of entropy or S/Q is branch-defining, not
  predictive.

## What would change this answer

- **A computed reheating closure.** It is physics, not a principle. It moves the windows (M2-Q is robust on P0 and
  marginal on GM), so it should be computed before any principle is judged.
- **A stated GRUT law about the initial state of the closed whole** that is not a typicality measure, gives s_* (and
  g_h) directly, and is written down before it is evaluated against these windows. None is on the table.
