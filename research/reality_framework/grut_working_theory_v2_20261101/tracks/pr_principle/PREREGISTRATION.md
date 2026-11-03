# Can a GRUT principle fix the holder? Pre-registration

Written and committed **before** any of the candidates below was evaluated against the M2 windows. Claude, for the
owner (issue #13).

**The question.** Can a principle determine any of the following without being tuned to the desired dark-matter
branch?
- the holder strength c0;
- the holder phase δ;
- the spectator's initial amplitude s_*;
- the reheating closure: Pauli blocking local (P0) or box-mean (GM).

## The targets (from the end-to-end pass)

| target | c0 | δ | closure |
|---|---|---|---|
| M2-Q | 1.40–1.50 | −4…+4° | P0 |
| GM pocket | 1.485–1.52 | +2…+3.5° | GM with fragment spread |
| M2-X | 1.45–1.60 | −2…+6° | GM |
| M1r's own history: Q-matter at e = 10⁻⁶ | ≳ 1.10 | — | — |

The holder strength is defined by μ_eff(slice) = μ + c0(10i − μ), with μ = 7.159456 H_i.

## The rule, fixed in advance

A principle **determines** a parameter only if all four hold:
1. It is stated without reference to S/Q, the branch, Ω_DM or the windows.
2. It can be evaluated on the model as it stands.
3. It returns a value or a narrow range. A broad distribution with the window in its tail does not count.
4. It is not one of many principles tried until one fits. Every candidate below is reported, hits and misses alike,
   and the expected number of chance hits is quoted.

A principle that is **branch-defining by construction** (one that extremises S/Q itself) is reported, but it does not
count as a non-tuned determination.

## The candidates (all will be reported)

- **A. Closure as defined** (unitarity, no exterior, purity of the whole). Does it constrain c0, δ, s_* or the
  closure?
- **B. Conjugation-even principles and δ.** The preparation theorem says such principles cannot fix the sign. Test
  whether a CP-even extremum, maximal CP torque (|sin| = 1 of the held relative phase), fixes |δ|.
- **C. Typicality of the state for s_*.** Starobinsky–Yokoyama equilibrium, a Hartle–Hawking ground state, or the
  Boltzmann measure behind Poincaré recurrence. Compute P(window).
- **D. Extremal yield.** Maximise the charge per energy ε, or |N_X|, over (c0, δ) at fixed couplings and history.
  Where is the maximum?
- **E. Edge of stability.** The held hidden state sits at the threshold where a held mode turns tachyonic (Y
  condensation). Where is that threshold in c0?
- **F. Thermodynamic extremum** of S/Q, maximum or minimum. Branch-defining by construction; reported only.
- **G. The reheating closure.** Not a free parameter: Fermi statistics fixes it. Ask what decides between local and
  box-mean blocking, and whether a computation, not a principle, settles it.
- **H. Observer selection.** Does either branch change anything observers need?
- **I. A UV fixed point for g_h** (asymptotic-safety type). Can it fix c0 = |g_h| s_*/(12.30 H_i) without s_*?

Data used: `tracks/e2e_chain/results/scan_*.json` (the (c0, δ) scans on the final histories), the prepared-state
solver, and `tracks/e2e_chain/skeptic/results/sk_natural.json`.
