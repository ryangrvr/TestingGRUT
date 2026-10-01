# S2-1b RESULT — can locality be selected without supplying graph distance, local dimension or k?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-1b. Pre-registered in REPAIR 01, before the run.

**Files:** `probes/S2-1b/s2_1b_criteria.py`, with log `probes/S2-1b/s2_1b_criteria.log`.

**Wall attacked:** D.

**Method.**
- One abstract Hermitian H on ℂ^N is scored in several candidate TPSs. A TPS is a unitary frame plus a grouping of
  slots into factors.
- Every score is computed from the exact-factor-support decomposition of H, which is LU-invariant. So each score is a
  TPS invariant, not a basis artefact.
- Criteria:
  1. minimal k;
  2. sparsest interaction graph (edge fraction);
  3. Lieb–Robinson (mean normalized ‖[P_a(t₀), Q_b]‖², t₀ = 0.3);
  4. MDL (dimension of the support-closed operator space);
  5. stability;
  6. locality-max (‖H‖² fraction on ≤ 2-factor supports);
  7. predictive autonomy (‖H‖² fraction on single-factor supports).

**Labels:**
- **LOCALITY NOT SELECTED — CRITERION-PRICED (D)**;
- the dynamical / operational branch (Lieb–Robinson, autonomy) is **ACCESS-PRICED (A)**;
- **STABILITY MEASURE-PRICED**;
- **EMERGENT BUT NONUNIQUE**;
- KNOWN RESULT IMPORT: Cotler–Penington–Ranard (local uniqueness given k and d); Zanardi (observable-induced TPS) —
  SECONDARY.

## 0. Verdict

**0. Without a supplied number of factors or local dimension, minimal k is ill-posed.** For a random H on ℂ¹⁶:

| factorization | 1 factor | 4⊗4 | 2⊗8 | 2⊗2⊗4 | 2⁴ |
|---|---|---|---|---|---|
| k | 1 | 2 | 2 | 3 | 4 |

The trivial factorization always wins. Locality-max has the same defect: any 2-factor TPS scores 1.0. Both need
"number of factors" or "local dimension d" supplied. "Finest = prime factorization" is itself an inserted rule.

**Case A — the local dimension conflicts.** N = 16, with H = h₁₂ + h₃₄ + 0.05·X₂X₃ (generic strong pair terms and a
weak bridge).

| TPS | k | edges | MDL | LR | loc2 | autonomy |
|---|---|---|---|---|---|---|
| qubits 2⁴ | 2 | 3/6 | **39** | 0.0815 | 1.0 | 0.506 |
| ququarts 4⊗4 | 2 | 1/1 | 255 | **0.0002** | 1.0 | **0.9997** |
| 2⊗8 | 2 | 1/1 | 255 | 0.1579 | 1.0 | 0.809 |

The criteria split:
- **qubits** win on sparse graph and MDL;
- **ququarts** win on Lieb–Robinson and autonomy;
- k and locality-max tie.

The split is a **scale** split. Structural criteria pick the fine factorization. Dynamical / operational criteria
pick the factorization whose factors are weakly coupled, a choice set by the coupling g. Which is "the" subsystem
structure depends on a supplied scale.

**Case B — same d (qubits), same number of factors (5); a Clifford-rotated frame.** H is a ZZZ 3-body chain plus
X fields. Of 1500 random Clifford frames, one strict disagreement pattern was found:

| TPS | k | edges | MDL | LR | loc2 | autonomy |
|---|---|---|---|---|---|---|
| computational | 3 | **7** | 159 | **0.074** | 0.700 | **0.700** |
| Clifford-rotated | 3 | 8 | **141** | 0.104 | **0.852** | 0.227 |

- MDL and locality-max pick the rotated frame.
- Graph sparsity, Lieb–Robinson and autonomy pick the computational frame.
- k ties.

**Even with d and the number of factors fixed, different criteria pick different TPSs.**

**Positive control.** For the transverse-field Ising chain (2-local nearest-neighbour), **no** strict disagreement
appeared in 1500 random Clifford frames. The computational frame is dominated on no criterion. **The criteria agree
when H already has a dominant, strongly local frame.** This matches CPR local uniqueness, which also needs k and d
supplied.

**Stability.** Perturb H → H + εV, with ‖V‖ = ‖H‖. The answer depends on the measure chosen for V:

| V drawn from | effect |
|---|---|
| GUE (TPS-neutral) | exact criteria collapse: qubit k → 4 already at ε = 10⁻³, while the ququart k stays ≤ 2 trivially. Exact-k stability therefore favours the *fewest factors* |
| 2-local in the qubit TPS | qubit k stays 2; MDL 66 |
| local in the ququart TPS | qubit k stays 2; MDL 39 |

A perturbation measure that is local in a chosen TPS is circular. Under every measure tried, weighted criteria survive
but still disagree: loc2 vs autonomy. **STABILITY IS MEASURE-PRICED.**

## 1. What selects the objective function?

Nothing in H. Every criterion is a functional F(H, TPS). Choosing F is an extra input. The criteria fall into two
clusters:

- **Structural** (k, graph, MDL, loc2). These require supplying d and the number of factors, or a "finest" rule.
  They are **D-items inserted** as a description-length or locality convention.
- **Dynamical / operational** (Lieb–Robinson, autonomy). They ask which factors can be **isolated, predicted and
  controlled separately** over a time scale t₀ or a coupling scale g. This is an **A-item**: what an agent can treat
  as a separate system. It also needs a supplied scale.

The clusters agree only when H has a dominant local frame, which is itself the D-property to be explained.

## 2. Structure inserted vs forced

| Item | Status |
|---|---|
| the TPS, given F, d and the number of factors (large n, dominant frame) | **forced** (CPR; numerics agree) |
| the objective F | **inserted**: D (structural) or A (operational) |
| local dimension / number of factors | **inserted** (minimal k is degenerate without them) |
| the scale (t₀, g) at which autonomy is judged | **inserted** |
| the robustness measure | **inserted** (a measure) |

## 3. Accounting against T2-1′

- **D is not eliminated. It is sharpened.** The D-wall contains an inserted **objective functional plus a local
  dimension plus a scale**.
- Part of what looked like D ("locality") is **A in disguise**: autonomy and light-cone criteria pick out subsystems
  an agent can separately prepare and predict.
- The **D/A boundary is itself criterion-dependent**. This is the first concrete sign that D and A may not be
  independent walls. The joint-forcing question (does primitive structure force D and A together?) is exactly where
  the criteria clusters must be reconciled.

**Status: S2-1b COMPLETE — locality not selected without a supplied objective, local dimension and scale. Criteria
disagree across factorizations (Case A: qubits vs ququarts) and within a fixed factorization (Case B: Clifford
frame). They agree only for H with a dominant local frame. Stability is measure-priced. D sharpened, partly → A.**
