# V0-2 TARGET SPEC — P-15 (Born ⟺ martingale ⟺ decomposition-independence chain)

**Definitions, assumptions, target statements and acceptance tests ONLY. No proof.**

**Source:** frozen `scout-0 @ ab2da47`, `playground/SCOUT_0/probes/P15_RESULT.md` (blob `0e0bb08920de…`).

**Lines inspected during extraction:** a heading grep, then lines 49 – 52, 56 – 70, 76 – 90 and 110 – 113 (generator,
theorem statement, controls A / B table and control F).
- Lines 56 – 63 include the original's one-paragraph proof sketch. It was **seen by the orchestrator but is NOT
  reproduced here.**
- No P-15 script or log was opened.

**After this spec is committed, the following are SEALED until the reproduction is committed:**
- `P15_RESULT.md`;
- `p15_branch_hitting.py`;
- `p15_leak_check.py`;
- their logs;
- all other P-15 material.

## 1. Setting (definitions)

**Branch coordinate.** p_t ∈ [0, 1] obeys the Itô SDE dp_t = b(p_t)dt + σ(p_t)dW_t.

**Generator.** 𝓛 = b ∂_p + ½σ² ∂_p², with σ > 0 on (0, 1) and σ(0) = σ(1) = 0, i.e. zeros only at the branch endpoints.

**Hitting / fixation weight.** h(p) := P_p(outcome = 1), where "outcome 1" means p_t → 1 (absorbed or fixed at 1),
either in finite time or asymptotically.

**SUPPLIED (not derived).**
- The identification p = |α|² is SUPPLIED.
- The stochastic law, including the probability carried by dW, is SUPPLIED.
- The chain does **not** derive probability itself.

## 2. Target claims (to be established or rejected independently)

**T1 — backward equation.** h solves b h′ + ½σ² h″ = 0 on (0, 1), with h(0) = 0 and h(1) = 1, under appropriate
absorption / asymptotic-fixation hypotheses. **State exactly which boundary assumptions are needed.**

**T2 — Born ⇒ zero drift.** h(p) = p for all p ∈ (0, 1) implies b = 0, via the backward equation. **Expose the
regularity used.** If the conclusion is only almost everywhere, or only on accessible states, say so.

**T3 — zero drift ⇒ Born.** If b = 0, then p_t is a bounded martingale and h(p) = p.
- Handle both finite-time absorption **and** asymptotic fixation (endpoint approached only as t → ∞).
- Do not assume a finite hitting time.
- Expose the condition that forces the martingale limit into {0, 1}.

**T4 — the meaning of "Born ⟺ martingale".** Distinguish:
- b = 0 as a generator statement;
- p_t being a martingale from every interior initial condition;
- exceptional or inaccessible sets.

Give the strongest statement actually justified.

**T5 — decomposition independence.** For a finite decomposition {(w_i, p_i)} with Σw_i = 1 and p̄ = Σw_i p_i, the
outcome-1 frequency is Σw_i h(p_i).
- **Claim:** this depends only on p̄, for every finite decomposition, **iff h is affine on [0, 1]**.
- Use the **minimum** regularity of h needed. The original states "midpoint-affine plus bounded".
- With h(0) = 0 and h(1) = 1, this gives h(p) = p. Then combine with T2 / T3.

**T6 — complete chain.** Determine whether the correct theorem is

  Born hitting law ⟺ zero drift ⟺ bounded martingale branch coordinate ⟺ decomposition-independent outcome statistics

under **one** common, explicitly stated hypothesis class. If not, state a set of conditional implications. **This
assumption audit is primary.**

**The original's stated hypothesis class (for audit):** σ² bounded away from 0 on compact subsets of (0, 1), and b
locally bounded. The original's corollary claims: "Born is the unique hitting law under which outcome statistics are a
function of the density matrix."

## 3. Exact / analytic controls (acceptance tests)

**C1 — driftless Wright–Fisher type.** σ²(p) = 2Dp(1 − p), b = 0. Then:
- derive h(p) = p;
- derive or verify E_p[τ] = −[p ln p + (1 − p) ln(1 − p)]/D, when the endpoints are exit boundaries (finite-time
  absorption).
- Also: for the variant σ(p) = 4√κ·p(1 − p) (natural boundaries, asymptotic-only absorption), the original claims h = p
  still holds.

**C2 — nonlinear-drift hostile control.** b(p) = λp(1 − p)(2p − 1), σ²(p) = 2Dp(1 − p), with frozen λ = 2, D = 1.
- Derive the scale-function expression for h.
- Evaluate h at p ∈ {0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9} to at least 5 decimals.
- **Original table to check:**

| p | 0.1 | 0.25 | 0.4 | 0.5 | 0.6 | 0.75 | 0.9 |
|---|---|---|---|---|---|---|---|
| h(p) | 0.07793 | 0.21955 | 0.38390 | 0.50000 | 0.61610 | 0.78045 | 0.92207 |

  Key check: h ≠ p away from the symmetry points.

**C3 — decomposition counterexample (under C2).** Compare pure p = 0.4 with a 50/50 mixture of p = 0.2 and p = 0.6
(same p̄).
- **Original:** pure h = 0.38390; mixture 0.39271. So the frequencies differ.
- Under C1, both are 0.4.

## 4. Hostile assumption audit (required)

1. σ has an interior zero.
2. A boundary is inaccessible.
3. Absorption / fixation is not almost sure.
4. b ≠ 0 only on a polar / inaccessible set.
5. Is boundedness alone enough for the decomposition argument?
6. Does h = p characterise zero drift pointwise under the stated regularity?
7. Is p = |α|² derived anywhere? (It must remain SUPPLIED.)
8. Does the chain derive probability itself? (It must not: the stochastic law is supplied.)

**Any necessary scope repair is grade-bearing.**

## 5. Out of V0-2 scope (deferred; not criterion-2 obligations)

- Monte Carlo counts;
- control D Euler leakage;
- dephasing / QSD simulation percentages;
- the full Lindblad unravelling comparison;
- the multibranch extension.
