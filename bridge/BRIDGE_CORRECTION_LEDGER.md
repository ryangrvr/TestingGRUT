# BRIDGE CORRECTION LEDGER (Bridge-1)

**Ordering.** Corrections are additive and owner-directed. Historical script logs (`bridge/b1/*.log`) are kept exactly as
emitted. Where a log line's wording is broader than the corrected scope, the correction below takes precedence.

## BRIDGE REPAIR 01 — B1 scope wording (owner, after `add9f79`; none of these changes the verdict)

| ID | location | original wording | corrected reading | effect on the verdict |
|---|---|---|---|---|
| **BR1-01** | `B1_SIGMA_LOCAL_NET_RESULT.md` (chain row, memory-kernel row); `BRIDGE_ZOOM_OUT_01.md` §1; ledger BI-05; log line "Lanczos from ANY unit start vector" | "every q gives a connected tridiagonal chain frame" | **Every cyclic start vector q gives a full connected Jacobi / tridiagonal representation.** For a simple spectrum the cyclic vectors form an open, dense, full-measure set. Exceptional non-cyclic q (orthogonal to some eigenvector) terminate Lanczos early | none: the continuum / non-uniqueness result is unchanged |
| **BR1-02** | `B1_SIGMA_LOCAL_NET_RESULT.md` (sparsity row); `BRIDGE_ZOOM_OUT_01.md` §1 | "optimum of every factorization criterion" | **"optimum of the tested entry-sparsity / autonomy / response-factorization criteria."** No claim is made about arbitrary possible criteria | none |
| **BR1-03** | `B1_SIGMA_LOCAL_NET_RESULT.md` (C1-a row); `BRIDGE_ZOOM_OUT_01.md` §1; ledger BI-02 | "known theorem: path DLS" | the fact that the unit-weight path is determined by its Laplacian spectrum is a **KNOWN-RESULT IMPORT**. Bridge-1 gives no independent proof | none: still **CONDITIONAL COMPRESSION / CRITERION-PRICED** |

**Status:** B0 and B1 are **accepted provisionally** by the owner, with BR1-01 … BR1-03 applied.

## BRIDGE REPAIR 02 — B3 bookkeeping scope (owner, after `046c2d7`; neither changes a B3 verdict)

| ID | location | original wording | corrected reading | effect on the verdict |
|---|---|---|---|---|
| **BR2-01** | `B3_ACCESS_RESULT.md` §B3-13; `B3_ACCESS_COMPONENT_LEDGER.md` (sharpened residual); `BRIDGE_ZOOM_OUT_02.md` §§2, 7; ledger tally line | `A_res → A_seed ⊕ A_partition ⊕ A_resolution ⊕ A_time` (and "A_seed(+readout identification)") | **`A_res → A_interface ⊕ A_partition ⊕ A_resolution ⊕ A_time`, with `A_interface := (A_seed, A_readout)`.** B3 separately shows that both the seed and the readout are supplied, and does **not** derive the readout from the seed. In the classical observability representation the readout vector plays the seed's role. That identification is **not** established across the canonical classes (P-5 / P-6 separate the coupling seed from the downstream probe / readout). The closure stays downstream of the **seed component only** | none |
| **BR2-02** | `B3_ACCESS_COMPONENT_LEDGER.md`; `BRIDGE_ZOOM_OUT_02.md` §§2, 7 (wherever `A_closure = f(D, A_seed)` appeared) | `A_closure = f(D, A_seed)` | **`A_closure = f(D, A_seed; R_closure)`.** The P-5 vs P-6 hidden-sector control (dimension 16 vs 4) shows the rule matters | none: still **CONDITIONAL COMPRESSION** |

**Status:** B3 is **accepted provisionally** by the owner, with BR2-01 and BR2-02 applied. TRUE COMPRESSION remains 0.

## BRIDGE REPAIR 03 — B2 scope (owner, after `4d313fa`; no verdict changes, no S6 result downgraded)

| ID | location | original wording | corrected reading | effect on the verdict |
|---|---|---|---|---|
| **BR3-01** | `B2_HCORR_RESULT.md` §B2-7 | "Every Q = F(K), P = K·F(K) with F > 0 is stationary, because AΣ + ΣAᵀ = 0 **if and only if** P = KQ with [Q, K] = 0." | **"For every positive function F(K), the zero-q–p-cross covariance Q = F(K), P = K·F(K) is stationary. Thus the harmonic chain admits an infinite non-Gibbs stationary family."** That suffices for STATIONARITY DOES NOT SELECT GIBBS OR A UNIQUE TEMPERATURE. No characterization of *all* stationary covariances is claimed (cross blocks and spectral degeneracies are not treated) | none |
| **BR3-02** | `B2_HCORR_RESULT.md` §B2-6 (incl. B2-CUC-1); `B2_BOUNDARY_COMPONENT_LEDGER.md` HB-05; `BRIDGE_ZOOM_OUT_03.md` §6 | S6 "forward" described as a single temperature-ordering convention | **"S6 contains two declared notions of forward: a temperature-ordering sign for J and descent toward a supplied reference for D. Neither is selected by the time-symmetric dynamics itself."** J: f_J = ±J on L1 / L2. σ = −Ḋ: reference-relative to S_ref = T_b·diag(r, 1). At T_s = T_b the L1 / L2 ordering disappears for J, while the entropy-reference construction keeps its supplied reference | none; B2-CUC-1 now carries the split |

### B2-P1 review note (recorded separately)

> An independent logic review finds Proposition B2-P1 sound at its stated scope: self-adjoint trace-class covariance
> perturbations of the S6 Gibbs state, positive total covariance, fixed local observables, and the S6 purely absolutely
> continuous half-line operator. The proof establishes local covariance convergence, not global trace-norm convergence.

**Grade:** **BRIDGE-THEOREM REVIEWED AT STATED SCOPE**, not CANONICAL GRUT THEOREM.

**Status:** B2 is **accepted provisionally** by the owner, with BR3-01 and BR3-02 applied.

## BRIDGE REPAIR 04 — accounting atoms are not automatically physical primitives (owner, after `2b874a3`)

| ID | location | original wording | corrected reading | effect on the verdict |
|---|---|---|---|---|
| **BR4-01** | `B4_NINE_LAYER_COMPRESSION_MATRIX.md` §B4-11; `BRIDGE_ZOOM_OUT_04.md` §5 and Status; `BRIDGE_INFORMATION_LEDGER.md` BI-17 and summary | "atomic supplied-primitive count ↑" (six previously implicit items) | **"explicit supplied / declarative accounting items increase"** / "the dependency graph resolves into more atomic information items". The six are of different kinds: A_resolution and A_time are access / protocol declarations; R_closure is a mathematical convention / rule choice; the absolute coordinate value of H_epoch is gauge under global time translation, though the existence of a special-form boundary event is physical; Gibbs-vs-GGE is a supplied state-class / postulate choice; H_cross is genuine physical boundary data. **PHYSICAL COMPRESSION = 0; DEPENDENCY COMPRESSION > 0; EXPLICIT ACCOUNTING-ITEM COUNT ↑; PHYSICAL PRIMITIVE COUNT CHANGE = NOT ESTABLISHED** | none; no B4 scientific verdict changes |

**Status:** B4 is **accepted provisionally** by the owner, with BR4-01 applied.

## B5 — no repairs required

B5 introduced no correction to earlier Bridge files. The success bar was preregistered from the owner's B5-0 text before
any candidate was evaluated. Literature for P6 (heat-flow reversal by initial correlations) was originally graded
**KNOWN-RESULT IMPORT, cited from memory and not re-fetched**; it is upgraded by BR5-05 below.

## BRIDGE REPAIR 05 — final handoff / provenance cleanup (owner freeze ruling, after `2cba90c`; no scientific verdict changes)

| ID | location | original wording | corrected reading |
|---|---|---|---|
| **BR5-01** | `BRIDGE_1_HANDOFF.md` (B0 row); `BRIDGE_INFORMATION_LEDGER.md` BI-13; `GRUT_SCOUT_CROSSWALK_01.md` (additive note) | "the 10 supplied layers mapped"; BI-13 "layers 6–10" | **"the nine canonical supplied layers mapped."** Canonical GRUT has exactly nine supplied layers. Layer 9 combines the gravitational branch and the cosmological transport inputs. No tenth canonical layer |
| **BR5-02** | `BRIDGE_1_HANDOFF.md` header; `B2_HCORR_RESULT.md` scope | "B2-P1 was logic-reviewed by the owner" / "independent logic review (owner)" | **"B2-P1 received an independent logic review in the Bridge audit and was found sound at its stated scope; it has not received independent human peer review."** The grade stays **BRIDGE-THEOREM REVIEWED AT STATED SCOPE**, not CANONICAL GRUT THEOREM (the B2-P1 review note above is read with this attribution) |
| **BR5-03** | `BRIDGE_1_HANDOFF.md`; `BRIDGE_ZOOM_OUT_04.md` §7; `B4_NINE_LAYER_COMPRESSION_MATRIX.md` B4-8 | "Supplied in both programs: Lorentz / causal-cone structure" | **"GRUT explicitly supplies its universal causal-cone / Lorentz structure. SCOUT-2 found Lorentz structure NOT DERIVED / not selected within its premise envelope. The correspondence is convergence on a residual boundary, not a shared derivation."** No compression implied |
| **BR5-04** | `BRIDGE_1_HANDOFF.md`; `BRIDGE_ZOOM_OUT_05.md` §9 | "A change could come only from a quantum lift, a non-Gaussian physical bath, or gravity" | **"A change would require a materially enlarged or altered premise envelope, for example a physical quantum lift, a non-Gaussian physical bath, gravity, an infinite-dimensional / field-theoretic primitive structure, or a genuinely new selector principle."** Not exhaustive |
| **BR5-05** | `B5_PREDICTION_PAYOFF_RESULT.md` P6 | `KNOWN-RESULT IMPORT — cited from memory / not re-fetched` | **`KNOWN-RESULT IMPORT — PRIMARY-SOURCE VERIFIED`** (owner-checked against the publication pages): Partovi, PRE 77, 021110 (2008); Jennings & Rudolph, PRE 81, 061130 (2010); Micadei et al., Nat. Commun. 10, 2456 (2019). They support **only** the qualitative background that initial correlations can reverse energy-flow direction, **not** Bridge-1's classical-Gaussian quantitative formulas |

**Status:** B5 is accepted. **BRIDGE-1 FROZEN** by owner ruling with BR5-01 … BR5-05 applied.
