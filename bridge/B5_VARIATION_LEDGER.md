# B5 VARIATION LEDGER (supplied-input vector and per-candidate variation)

## B5-1 Residual-input vector Θ (still supplied after Bridge-1)

| block | supplied inputs |
|---|---|
| **Dynamics / substrate** | K and its coupling strengths (κ = K₁₁, g, bath pin c, dispersion); the allowed net within its class (L0-1b non-unique, B1); the generator class; optional nonlinear drift (on-site or not) |
| **Σ** | the local net / decomposition; the retained / bath split (not fixed by K, B1-3) |
| **H** | H_marginals (T_s, T_b, bare vs dressed form); H_cross (any PSD-admissible S–B block); the special-form boundary event (its coordinate value is gauge); the bath-state class (Gibbs vs GGE); temperatures |
| **A** | A_seed, A_readout, A_partition, A_resolution, A_time; R_closure wherever a closure appears |
| **Fenced** | lift / ħ / outcome / gravity / cosmology. **No candidate consumes them.** |

## Candidate rows: candidate → relevant free inputs → variation test → surviving relation

| ID | candidate | relevant free inputs | variation test (evidence) | surviving relation | terminal |
|---|---|---|---|---|---|
| **P1** | frame(net) = frame(drift) = frame(noise) | drift frame, noise frame, net frame (each a supplied declaration) | A: common frame → overlap 1.000000 / 1.000000. B: drift written on-site in an independently rotated frame → overlap with the net frame **0.349**; the model stays admissible (convex V, relaxes to 7e-4). C: an ordinary nonlinear Langevin network with site-local V₄ and noise is the same equation (`b5_payoff.log` P1) | only "if all layers are declared in one frame, they agree" | **CONSISTENCY-ONLY / STANDARD-STRUCTURE** |
| **P2** | X(∞) = E₁(0) − E₁_G + g·Q12_G − g·Q12(0) | κ, g, c, T_s, T_b, α, bath class | 5 chains × 4 preparations: numeric = closed form to ≤ 1.5e-4. The value moves with every input. The class boundary κ = 1.0 (bound state at λ = 0.231): X(T) wanders over [1.35, 2.00], so there is no limit. Comparator: a random 30-node network with a 3-site system satisfies ∫flux = ΔE_B to 1e-6 | an identity = **energy conservation** (+ return to equilibrium inside the a.c. class); its value is preparation-dependent | **STANDARD-STRUCTURE** (identity) / **INPUT-DEPENDENT** (value) |
| **P3** | equal-T offset X(∞) = ½T_b r² | T_b, κ, g, c, α | closed form X(∞) = −⟨E_int⟩_G (½ − α), confirmed in 5 chains × α ∈ {0, 1} (ratio ±0.4998 … ±0.5000). The ½ is the static equipartition identity E_S(prod) − E_S(Gibbs) = ½⟨E_int⟩_G, exact in 3 random 25-node networks with a 4-site system. "r²" → (K⁻¹)₁₂ takes values 0.339 / 0.241 / 0.167 / 0.319 / 0.224 | ½ survives only at α = 0 and T_s = T_b, and is a generic quadratic-network identity; r is model-specific; the sign reverses for α > ½ | **INPUT-DEPENDENT + STANDARD-STRUCTURE** |
| **P4** | A_closure = f(D, A_seed; R_closure) | seed, R_closure, readout | B3-1: the closure dimension changes with the rule (16 vs 4) and with the seed (24 / 12 / 2 classical). No measurement maps a closure dimension to an operational quantity without a further supplied readout | none | **NO-OBSERVABLE** (bookkeeping only) |
| **P5** | access-relative geometry | seed, partition, resolution | B3-4 + GS1 (canonical): geometry is selected under full access, underdetermined under single-site access, and the dynamic hierarchy beats static data at the boundary; the topology horizon is order-relative (B3-7) | no geometry invariant survives all access choices | **ACCESS-PRICED → INPUT-DEPENDENT** |
| **P6** | arrow sign at equal marginal temperature | H_cross over the full PSD-admissible set | extremal sweep at the S6 marginals: (1,1): X ∈ [−0.334, +0.673], **both signs**; (2,1): [+0.458, +1.881]; (10,1): [+7.58, +10.76]; (½,1): [−0.686, +0.025] (sign flip available) | at equal T, no sign survives. Where the sign is fixed, the boundary is a Cauchy–Schwarz / PSD bound \|Q12(0)\| ≤ √(Q11·Q22) with supplied marginals | **BOUNDARY-STATE-DEPENDENT**; the bound is a generic positivity bound (STANDARD) |
| **B5-9** | dimensionless ratios X/(g·Q12_G), X/(T_s − T_b) | α, T_s, T_b, κ, g | X/(g·Q12_G) = ½ only at (α = 0, T_s = T_b); at α = 0.7 it is −0.200; at (2, 1) it is 3.45 / 4.65 / 10.48. X/(T_s − T_b) = 1.169 / 1.120 / 1.050 (K-dependent) | no constant survives | **INPUT-DEPENDENT** |

## B5-8 Quotient test

| O | {O(θ) : θ ∈ Θ} | collapses to |
|---|---|---|
| frame agreement (P1) | any relative frame (overlap 0.35 … 1) | **a broad family** (unless a shared frame is declared) |
| X(∞) (P2) | ℝ-interval set by the marginals and the correlations | **a broad family**; the identity is energy conservation |
| equal-T offset (P3, P6) | [−0.334, +0.673]·T at the S6 marginals; model-specific scale | **a bounded region**, but its boundary is a PSD / Cauchy–Schwarz bound (excluded by the success bar) |
| closure dimension (P4) | {2, 4, 8, 12, 16, 24, 32, …} | **no observable** |
| reconstructed geometry (P5) | from unique to a non-isometric family | **a broad family** |
