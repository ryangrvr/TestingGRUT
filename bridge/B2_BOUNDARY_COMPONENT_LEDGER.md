# B2 BOUNDARY COMPONENT LEDGER

**Fixed throughout:** D = K_∞ (S6 parent); Σ = site 1 | bath; A_interface, A_partition, A_resolution and A_time (retained
2×2 block plus J, same times). No access variation.

| ID | component | SCOUT counterpart | canonical GRUT value | parent inputs | derived (and grade) | residual choice | classification | evidence |
|---|---|---|---|---|---|---|---|---|
| HB-01 | **H_sys** (system marginal) | H_corr\|Σ marginal part | T_s/2.3, T_s (bare Gibbs) | supplied | locally forgotten (LS-1; B2-P1) | T_s, bare vs dressed form | **SUPPLIED** (memory compressed downstream, retained in X_J) | B2-0, B2-7, B2-8 |
| HB-02 | **H_bath** (bath marginal) | environment state / freshness (D4) | T_b·K_BB⁻¹, T_b·I | Gibbs postulate + T_b | the local asymptotic state S_ref = T_b·diag(r, 1) (LS-1, theorem) | T_b; Gibbs vs GGE | **PARTIAL RELOCATION INTO ENVIRONMENT** (conditional on the Gibbs postulate + T_b); **THERMAL CLASS CONDITIONALLY DEFINED; TEMPERATURE SUPPLIED** | B2-7, B2-8 |
| HB-03 | **H_cross** (S–B correlations) | **H_corr\|Σ** proper (D3 / D6) | 0 (product) | supplied | X_J(∞) depends on Q12(0) exactly (bridge algebra) | the cross block | **SUPPLIED; INDEPENDENT; LOAD-BEARING** (the sign of X_J is set at equal T) | B2-1, B2-2 |
| HB-04 | **H_epoch** | special moment (Janus) | t = 0, where the cross blocks vanish | supplied | — | the epoch / form | **SUPPLIED**: selected by the preparation form; a gauge relative to the history otherwise | B2-5, B2-12 |
| HB-05 | time orientation (forward convention) | orientation NOT selected | L1 / L2 sign rule, σ = −Ḋ | convention | — | the sign convention | **NO BRIDGE / CONVENTION**: Σ(−t) = RΣ(t)R; the reversed preparation anti-relaxes | B2-4, B2-5, B2-6 |
| HB-06 | relaxation mechanism | Markov / fresh-bath mechanism (D4 / D5) | a.c. spectrum of K_∞ | D + infinite bath | local return to equilibrium for **any trace-class** Δ (B2-P1, bridge sketch) | — | **CONDITIONAL COMPRESSION** of memory downstream ("D COMPRESSES MEMORY OF H_corr"); not a derivation of H_corr | B2-3, B2-9 |
| HB-07 | stochastic-noise stationary state | — | L0-1e, Q = 2·diag(T_i) | T_i | stationary covariance (correlated) | T_i | **RELOCATION** of the reference state into the noise layer; initial H_cross erased, not selected | B2-10 |
| HB-08 | "low entropy" label | (SCOUT: not H_env, not low entropy) | — | — | — | — | **NOT THE PRICE** (the higher-entropy marginal-matched state relaxes; equal-entropy reversed states anti-relax) | B2-11 |
| HB-09 | "temperature difference" label | — | S6 L1 / L2 members | — | — | — | **NOT NECESSARY** (the equal-T transient ½Tr²; the correlation-only transient Tr²) | B2-1 |

## Count

| class | count |
|---|---|
| TRUE COMPRESSION | **0** |
| CONDITIONAL COMPRESSION | 1 (HB-06: memory compressed downstream) |
| RELOCATION | 2 (HB-02 partial into the environment; HB-07 into the noise layer) |
| SUPPLIED | 3 (HB-01, HB-03, HB-04) |
| NO BRIDGE / CONVENTION | 1 (HB-05) |
| NOT-THE-PRICE findings | 2 (HB-08, HB-09) |

**Correspondence (GRUT S6 ↔ SCOUT reviewed H_corr\|Σ):** the S6 product preparation is **one explicit member** of the
H_corr|Σ class. It plays the same role, so this is **RENAMING**. The two programs reach it independently. GRUT supplies an
exact closed-form witness that the correlation boundary, not temperature or entropy, is load-bearing.
