# Track sp_vector: the longitudinal hidden vector at the 373d37e leading point

25 October 2026 · Claude (independent review track, issue #13) · packet `373d37e`
(`origin/codex/grut-inflation-scale-scan-20261025`, `research/reality_framework/grut_inflation_scan_20261025/`),
read only through `git archive 373d37e … | tar -x` into the scratchpad. Our own code is in `code/`, outputs in
`results/`. No Codex code was run, and nothing was committed. Our earlier gp_matching model is reused for the
inflationary minimum, as the brief permits.

**This is one supplied tree-level working model. Nothing here is evidence that X, Y, Σ, φ, the vector or the
inflationary inputs exist. It is not a validation, an exclusion or a GRUT prediction.** All agents involved are
same-session AI, so this is not outside replication.

Labels: [D] derived, [C] computed, [Q] quoted, [E] estimate, [A] argued, [R] recalled.

## 1. Question and inputs [Q]

- Couplings: M_Y = 1e12, h = 1e11, μ = 5e11 GeV, κ = 1e-3, e = 1e-7, λ_Y = λ_Σ = 1, λ_X = 1e-3; charges q = 1, 3, 6.
- Leading point: H_i = 1.1211167304108865e13 GeV, c_H = 1, aligned inflationary minimum, Σ = −|Σ| e^{iπ/2}, zero
  velocities, instant radiation release (Codex CLAUDE_HANDOFF.md). Net comoving charge N_X = −238.649 (H_i³ units).
- Allowance: 0.0049296 pairs per net charge (1% of calibrated dark matter).
- Codex's vector statements: VECTOR_INITIALIZATION.json gives the Higgs norm F_g, m_A/H and an assumed Bunch–Davies
  (BD) Goldstone normalization, √P_G/H = 1/(2π), √P_α = H/(2πF_g), P_AL/H² = (k/H)²/[4π²(m_A/H)²]. MATHEMATICS.md
  says the evolving coupled vector/scalar modes must not be replaced by a constant-vacuum-mass abundance, and
  no final occupation is claimed.
- Our prior estimate (gp_leftovers, gp_critic) used Graham–Mardon–Rajendran (GMR) with the vacuum mass: 1053–1061
  pairs per net charge at 4.7e13 GeV, normalized to Y_Q.

## Result in brief

The vector budget at this point has two different pieces, and they behave very differently.

- **Piece A: the inflationary longitudinal mode (GMR type) with the evolving Higgs norm, computed [C].**
  - The superhorizon gauge angle is frozen, so the vacuum-mass GMR value is reduced by (F_vac/F_inf)² = 2.36e-6.
  - Result: 1.6e-6 to 7.1e-6 X X̄ pairs per actual net charge, or 1.0e-4 to 4.7e-4 normalized to Y_Q.
  - The allowance is 0.00493, so **this source alone does not fail the point**.
  - Our earlier constant-mass estimate (about 1060 pairs per net at 4.7e13 GeV) is withdrawn for this source.
- **Piece B: the release-induced Goldstone angle, estimated [C inputs, E].**
  - The compensating Σ charge winds the Higgs phase by about 2300 rad, and the winding depends on the initial
    conditions at 1.3–1.7e6 rad per rad.
  - Patch-to-patch differences become longitudinal vectors.
  - With natural inflationary noise, in either the source phase or the |Σ| modulus, the Hubble-scale part alone gives
    0.03–1.0 pairs per net charge (6–210 times the allowance). The k_* part is unresolved at the supplied phase.
  - **On these estimates the vector sector fails the point** unless both phase and |Σ| are prepared to ≲ 1e-3 per
    patch and ≲ 1e-6 at k_*. The model has no such preparation.
- **e ≳ 0.0287 is not a free remedy.** It removes only piece A, and it opens strong parametric vector production off
  the eccentric X orbit (Floquet exponent ≈ 0.15 per unit a).

## 2. Method (own code)

- `code/vmodel.py` [C]: inflationary minimum and 6×6 Hessian (gp_matching model), Higgs norm F_g² = 2Σq_i²|z_i|²,
  the phase-sector eigenmodes (Goldstone and two physical phases), the BD gauge angle α = H/(2πF_g), the light-phase
  spread, and the GMR scales a_* (H = m_vac) and k_* = a_* m_vac.
- `code/proca_modes.py` [C]: the longitudinal Proca mode with a time-dependent mass,
  S_k = ½∫dτ[z²|B′|² − a²m²|B|²], z² = a²m²/(k² + a²m²) (GMR action [R]), evolved as the pair (B, Π = z²B′) from the
  adiabatic vacuum at k/aH = 1000 during inflation, through release (a = 1), to 30 a_* in radiation. The late number
  is n_k = [|Π|²/2z² + a²m²|B|²/2]/ω − ½. 49 modes, k = 1e-7…1e-1 H_i. Mass profiles: P0 constant m_vac (the GMR
  control), P1 m² = m_vac² + (m_inf² − m_vac²)/a² (X part of the Higgs norm redshifting), P4 the da-averaged F_g²(a)
  measured in our own release, P2 an instant step to m_vac (stress test only), P5 constant m_inf.
- `code/vrel.c` [C]: our homogeneous X/Y/Σ release (conformal fields, Codex's dictionary, equations written
  independently; Dormand–Prince 5(4)). It adds the unwrapped phases θ_X, θ_Σ, the late Σ phase
  θ_Σ(∞) = θ_Σ(a_end) − ⟨N_Σ⟩/(ĥ²a_end), and the gauge-invariant Goldstone angle
  θ_G = Σ_i q_i 2|χ_i|²θ_i / Σ_j q_j² 2|χ_j|².
- `code/phase_transfer.py`, `code/phase_scan_smoothed.py` [C]: linear response of θ_Σ(∞), N_Σ and ε to the Σ phase and
  to |Σ|; Gaussian ensembles of the Σ phase (σ = 1e-6 … 0.489, N = 128 each); the noise-averaged large-scale response.
- `code/floquet_vector.py`, `code/floquet_band.py` [C]: Floquet exponents of transverse vector modes in the late X
  orbit, as the cost check for a larger e.
- `code/vector_ledger.py` [C, E]: assembles the ledger. `code/f2profile.py` extracts F_g²(a).

## 3. Results

### 3a. Inflationary background and Codex's initialization [C]

| Quantity | Codex VECTOR_INITIALIZATION.json (nearest rows 8.72e12, 1.33e13 GeV) | This track (at 1.1211e13) | Verdict |
|---|---|---|---|
| F_g/H_i | 34.869, 34.811 (log-interpolated 34.834) | 34.8322 (F_g = 3.905e14 GeV) | CONFIRMED |
| m_A/H_i during inflation | 3.4869e-6, 3.4811e-6 | 3.4832e-6 (m_A = 3.905e7 GeV) | CONFIRMED |
| BD gauge angle √P_α | 0.004564, 0.004572 | 0.004569 rad per ln k | CONFIRMED |
| √P_G/H, P_AL normalization | 1/(2π); (k/H)²/[4π²(m/H)²] | same Goldstone-equivalence normalization (initial state of proca_modes.py) | CONFIRMED |

- Vacuum mass m_A = 6eh = 6.0e4 GeV, so F_inf/F_vac = 650.8 and (m_vac/m_inf)² = 2.3607e-6.
- Goldstone weights during inflation (q_i 2|z_i|²/F²): X 0.936, Y 0.0177, Σ 0.0018. The X field is the Higgs.
- Phase masses: 0, 0.3716 H², 7.62 H². The light physical mode (m² = 0.3716 H², ν = 1.3706) is almost purely the Σ
  phase. Its BD amplitude at horizon exit is 0.260 rad per ln k; at k_* it has decayed to 0.0758 rad per ln k by the
  end of inflation. Its stochastic rms is 0.522 rad (sp_yledger: 0.52 [Q]), and 0.489 rad comes from modes above k_*.
- a_* = 1.367e4 (H = m_vac), k_* = 7.32e-5 H_i, which exited 9.52 e-folds before the end of inflation.
- The vector is heavy during inflation (m_A > H_i) only for e > 1/(F_g/H_i) = 0.0287 (gp_critic: 0.029 at
  4.7e13 GeV [Q]).

### 3b. Piece A: the inflationary (BD) longitudinal mode with the evolving mass [C]

Physics [D, A]. A constant gauge rotation δθ_i = q_i α is an exact k = 0 solution for any background, so the
superhorizon gauge angle is frozen: the gauge-invariant longitudinal amplitude B_L = kα/e (comoving) does not change
while the Higgs norm falls from F_inf to F_vac. The late energy m_vac²B_L²/2 is therefore smaller than the
constant-vacuum-mass GMR value by (F_vac/F_inf)².

| Profile | N_A (H_i³) | n_k/n_k(P0) at k ≤ 1e-4 | Pairs per actual net charge | Pairs per net, Y_Q-normalized |
|---|---:|---:|---:|---:|
| P0 GMR control (constant 6e4 GeV) | 204.3 | 1 | 0.856 | 56.5 |
| P1 X redshift | 1.68e-3 (3.77e-4 from k < 1/3000) | 2.357e-6 | 7.1e-6 (1.6e-6) | 4.7e-4 (1.0e-4) |
| P4 measured F_g²(a) | 1.31e-3 (3.79e-4 from k < 1/3000) | 2.358e-6 | 5.5e-6 (1.6e-6) | 3.6e-4 (1.05e-4) |
| P2 instant step (unphysical) | 1.89 (3.82e-4 from k < 1/3000) | 2.361e-6 | 7.9e-3 | 0.52 |
| P5 constant m_inf (orientation) | 7.88 | — | 0.033 | 2.2 |

- **Control.** P0 gives Ω_A/Ω_DM = 3.74e7 with instant reheating, against the GMR formula's 3.97e7 (ratio 0.94).
  The spectrum peaks at k ≈ 1.4 k_*.
- **The frozen-angle scaling is exact for superhorizon modes.** For k ≤ 1e-4 the ratio is 2.357e-6, against
  (m_vac/m_inf)² = 2.3607e-6.
- Modes that enter the horizon while X still dominates the Higgs norm (k ≳ 1/1000) carry more in P1 and P4. In the
  single-Higgs reduction they count as vectors. In the three-field theory they are mostly X-phase quanta, because the
  Σ and X phase waves separate once they are subhorizon. Counting them all is an upper bound.
- P2 is not physical. The conformal X amplitude does not vanish at release, so the X part of the Higgs norm cannot be
  removed in one step. A step forces the k²-suppressed Higgs charge of each mode into the light Higgs and creates UV
  quanta. The superhorizon part (3.8e-4) is unchanged.
- Convergence: n_k changes by ≤ 2.4e-3 between 20 a_* and 30 a_*; the normalized Wronskian error is ≤ 1.4e-10.
- **Result: 1.6e-6 to 7.1e-6 pairs per actual net charge (3e-4 to 1.4e-3 of the allowance); 1.0e-4 to 4.7e-4
  Y_Q-normalized (2–9% of the allowance).** Our earlier constant-mass estimate overstated this piece by
  (F_inf/F_vac)² ≈ 4.2e5. Rescaled to 4.7e13 GeV (F_g/H = 34.71), the earlier 1053 becomes 1.4e-4.

### 3c. Decay channels [D]

- Γ(A → X X̄) = e²m_Aβ³/(48π) = 3.98e-12 GeV (β = 1 to 1e-11), a lifetime of 1.65e-13 s.
- H = Γ at T_vis = 1683 GeV; the decay-weighted epoch is T_vis ≈ 1343 GeV (sp_yledger 1.3 TeV [Q]).
- Other channels are closed. Y (M_Y = 1e12) and the Σ radial mode (1.4e11 GeV) are too heavy, and the Σ Goldstone
  is eaten. A → φφ is forbidden by Bose symmetry for J = 1. Kinetic mixing is forbidden by dark C in the 373d37e EFT
  audit [Q].
- **Each vector therefore gives one X X̄ pair of 3e4 GeV quanta.** At decay it carries 15.2 E_X per vector
  (sp_yledger [Q], same mass and width).
- Vectors made thermally from an X plasma (X X̄ → A → X X̄) shuttle pairs and add none. Only vectors whose energy
  comes from outside the X sector add pairs.

### 3d. Piece B: release-induced Goldstone angle (coupled scalar/vector transfer) [C inputs, E conversion]

This is the "evolving coupled vector/scalar modes" that Codex flags. It is a separate source from piece A, and
it is not in any homogeneous ledger.

**Mechanism [D, C].**
- In every homogeneous patch the total gauge charge is zero, N_X + 3N_Y + 6N_Σ = 0. The X charge made at the point
  is compensated mainly by the Σ condensate: N_Σ = 38.75 (a = 600–1000) and 39.67 (2000–3000), tending to
  −N_X/6 as N_Y → 0.
- A charged Higgs condensate rotates, θ_Σ′ = −N_Σ/(2|χ_Σ|²); late on this is N_Σ(H_i/h)²/a² ≈ 4.9e5/a². The
  Σ phase winds −1824 rad by a = 1000 and −2311 to −2324 rad in total (a_end = 1000 and 3000 estimators; rtol
  1e-10 against 1e-12 changes it by 6e-4 rad).
- Gauge-invariant statement [D]. With zero total charge the charge-weighted phase velocity vanishes, so the
  Goldstone angle changes only as the Goldstone direction rotates: θ_G′ = Σ_i w_i′θ_i. The direction rotates from X
  (w_X = 0.936 during inflation) to Σ. The term w_X r_X falls as 1/a, so θ_G → θ_Σ(∞)/6.
  - Numerical check on two runs 1e-7 rad apart: the window-averaged δθ_G at a = 2000–3000 is 0.0216, against
    δθ_Σ(a)/6 = 0.020–0.023 and δθ_Σ(∞)/6 = 0.028.
- While m_A < H, which holds until a_* = 1.37e4, the gauge field cannot respond on superhorizon scales.
  Patch-to-patch differences of θ_Σ(∞) therefore become longitudinal amplitude, B_L(k) = (k/e)[α + δθ_Σ(∞)/6].
  The energy F_vac²k²δθ_G²/2 is independent of e.
- **Conversion rule [C].** Using the P0 transfer, N_A = ∫dlnk [dN_A^{P0}/dlnk] P_θG(k)/2.974², because GMR with the
  vacuum mass is a Goldstone-angle power of (H_i/2πF_vac)² = 2.974² rad² per ln k. The 1% allowance (at
  |N_X| = 238.6) therefore needs **δθ_G(k_*) < 0.226 rad per ln k**. The BD angle α = 0.0046 is 49 times below this.
  The same rule reproduces P0's subhorizon spectrum, dN/dlnk = 0.0127/k.

**The transfer is extreme [C].**

| Response at the supplied point | Value | Linear range |
|---|---:|---|
| dθ_Σ(∞)/dψ (source phase) | 1.31e6 rad/rad (a_end 1000); 1.70e6 (a_end 3000) | \|Δψ\| ≲ 1e-6 rad |
| dθ_Σ(∞)/d ln\|Σ\| | −3.65e6 | ≲ 1e-7 |
| dε/dψ | 257 (Codex ≈ 260, sp_yledger 257 [Q]) | ≲ 1e-6 |
| dθ_Σ(∞)/dc along the light Σ-radial mode (m² = 1.11 H²) | −5.96e6 per canonical H_i | — |
| same, X-radial mode (m² = 1.99 H²) | 3.45e5 | — |
| same, Y-radial mode (m² = 21.6 H²) | 2.43e6 | — |

| Σ-phase ensemble (N = 128) | ⟨ε⟩ | std ε | ⟨N_X⟩ | std N_Σ | std θ_Σ(∞) [rad] | min\|Σ\|/vacuum |
|---|---:|---:|---:|---:|---:|---:|
| σ = 1e-6 | −0.2406 | 2.5e-4 | −238.6 | 0.039 | 1.3 | 0.118 |
| 1e-4 | −0.2411 | 0.0049 | −238.0 | 0.86 | 73 | 0.118 |
| 1e-3 | −0.2441 | 0.0106 | −238.9 | 1.8 | 158 | 0.107 |
| 0.01 | −0.2495 | 0.0245 | −242.7 | 4.4 | 194 | 0.047 |
| 0.03 (Codex's test) | −0.2498 | 0.0428 | −242.8 | 7.8 | 314 | 1.5e-3 |
| 0.1 | −0.189 | 0.098 | −181.0 | 15.7 | 1071 | 1.4e-4 |
| 0.26 | −0.082 | 0.131 | −76.6 | 20.5 | 1483 | 1.2e-5 |
| 0.489 (natural, above k_*) | −0.058 | 0.122 | −56.0 | 19.3 | 1379 | 2.1e-5 |
| **radial noise only, phase fixed** (H/2π on both light radial modes) | **−0.048** | 0.109 | −45.9 | 17.3 | 1372 | 2.3e-5 |

- Codex's σ = 0.03 mean is reproduced (−0.2498 ± 0.0038; Codex −0.2491 ± 0.0033). Our σ = 0.489 net
  (−0.058 ± 0.011) matches sp_yledger's −0.053 at 0.52 rad [Q].
- **New: Hubble-scale radial fluctuations alone decorrelate the release.** The light inflationary radial mode is
  almost pure |Σ| (m² = 1.11 H², u_Σ = 0.9998). Its BD amplitude is a 26% fluctuation of |Σ| per Hubble patch.
  With the phase held exactly at π/2 it reduces the net charge fraction to −0.048 ± 0.010, the same as the natural
  phase spread. A phase-preparation interaction alone does not keep ε ≈ 0.24; |Σ| must be prepared as well
  (sp_yledger's 2e-7 radial sensitivity [Q] is the linear version of this).
- For σ ≥ 0.03, 7–40% of members bring |Σ| below 1% of its vacuum value, and the minimum is 1e-5. Σ-phase slips, and
  hence string formation at release, are possible. They are not evaluated.

**B1: the large-scale piece (k ~ k_*) [C inputs, E conversion].**
- If the map were linear, the allowance would need δψ2(k_*) < 1.0e-6 rad and δ ln|Σ|(k_*) < 3.7e-7 per ln k.
- Natural inflationary amplitudes at k_* in linear extrapolation:
  - the ψ2 phase (0.0758 rad) is 7e4 times over;
  - the Σ-radial mode (2.6e-3 canonical) gives δθ_G = 2.6e3 rad, 1.1e4 times over;
  - the X-radial mode gives 0.76 rad, 3.4 times over;
  - the Y-radial mode gives 0.023 rad, 0.1 times.
- With natural Hubble-scale noise the map is saturated, so the operative quantity is the noise-averaged slope
  d⟨θ_Σ(∞)⟩/dψ0. At the supplied ψ0 = π/2 it is unresolved:
  - −1467 ± 1213 and 12 ± 471 rad/rad (common random numbers, ±0.05 and ±0.15);
  - 65 ± 132 rad/rad over ±0.6 (640 releases); the ⟨ε⟩ slope there is 0.014 ± 0.011.
  - ⟨ε⟩ is nearly flat there: −0.052 and −0.054 at ±0.15.
  - The criterion at the natural net charge (|N_X| = 56) is |slope| ≲ 9 rad/rad.
- **B1 is NOT RESOLVED at ψ0 = π/2.** With the scan slope:
  - central: δθ_G(k_*) = 0.82 rad, 0.24 pairs per net charge (49 times the allowance);
  - 2σ upper: 6.1 pairs per net charge (1.2e3 times);
  - the slope is consistent with zero.
  - At a generic ψ0 the structural size ⟨θ_Σ(∞)⟩ ≈ 9.6e3⟨ε⟩ implies slopes of a few hundred rad/rad, which would
    fail by 1e2–1e3.

**B2: the Hubble-scale piece [E, order of magnitude].**
- Model:
  - patch-to-patch white noise in N_Σ, set at a_s = 50–100 with the ensemble variance and correlation coefficient
    c_w = 0.1–0.8 (P_N = c_w Var(k a_s)³);
  - each mode accumulates δθ_Σ = δN_Σ[θ_Σ(a_k) − θ_Σ(a_s)]/N_Σ until horizon entry at a_k, from the nominal history;
  - it then counts (ĥ²/2) a_k P_θΣ quanta per ln k (the P0-consistent rule above).

| Patch noise | Pairs per net charge | Times the allowance |
|---|---:|---:|
| natural phase (0.489 rad), net \|N_X\| = 56 | 0.030–1.0 | 6–210 |
| radial noise only, phase fixed, net 46 | 0.030–1.0 | 6–207 |
| Codex's σ = 0.03, net 243 | 0.0011–0.039 | 0.23–7.9 |
| σ = 1e-3, net 239 | 6.2e-5–2.1e-3 | 0.013–0.44 |

- Each vector carries 15.2 E_X at decay. B2's upper natural value (~58 vectors per H_i³) would add about 900 E_X to
  sp_yledger's 510–542 E_X injection, so it also worsens the equilibrium relic.

### 3e. What would remove the vector source, and at what cost

**Piece A needs no remedy.** At e = 1e-7 it is ≤ 1.4e-3 of the allowance per actual net charge, and ≤ 0.09 with
the Y_Q normalization.

**Larger e: heavy vector during inflation (e ≥ 0.0287) [D, C].**
- This removes piece A altogether: m_inf ≥ H_i, and m_vac = 6eh ≥ 1.72e10 GeV.
- [A] A vector heavy through the release (m_A ≫ H and ≫ the X orbit frequency) would also screen piece B, because
  patch-to-patch Higgs-phase differences are then gauged away adiabatically instead of freezing into B_L.
- **Cost: parametric vector production off the X orbit** (`floquet_vector.py`, `floquet_band.py`) [C].
  - The late X condensate is an eccentric conformal orbit: r_min/r_max = 0.15, radial frequency ω = 2.24 per unit a.
  - The vector mass a²m_A² = 2e²|χ_X|² is therefore modulated by a factor 44.
  - Transverse Floquet exponents (per unit a; growth exp(2μΔa)):

    | e | m_A/ω | μ_max | 2μ × 400 |
    |---:|---:|---:|---:|
    | 1e-7 (stated) | 1.8e-6 | 1.6e-12 | 1e-9 |
    | 1e-3 | 0.018 | 1.85e-4 | 0.15 |
    | 3e-3 | 0.055 | 1.66e-3 | 1.3 |
    | 0.01 | 0.18 | 0.0165 | 13 |
    | 0.0287 | 0.53 | 0.152 | 122 |
    | 0.05 | 0.92 | 0.224 | 179 |
    | 0.1 | 1.8 | 0.057 | 45 |
    | 0.3 | 5.5 | 0.171 | 136 |
    | 1 | 18 | 0.066 | 53 |

  - Narrow first band: μ ≈ 185e², which matches the analytic narrow-resonance estimate to better than 1%. Staying
    resonance-free over the release needs **e ≲ 2.6e-3**, where the vector is light and piece B is unscreened.
  - Every tested e from 0.0287 to 1 is strongly resonant. Vectors would be produced explosively until backreaction
    and would then decay to X X̄. This conclusion is the transverse, X-only Floquet result: the growing Σ mass term
    (which detunes) and the longitudinal modes are not included [E].
- Other consequences [A, not computed]:
  - for e > 0.118 (m_A > ΔM_Y, m_σ/2), sp_yledger's σ → AA and Y_H → Y_L A channels close, so the σ quanta
    (N_σ = 4.9) need another decay;
  - the Z6 string type changes (β = m_σ²/m_A² = 2/(36e²));
  - the screened gauge correction to the object binding, Q²/(36h²R³), does not depend on e.
- **There is no e-window that removes both pieces without a new resonance problem.**

**Phase and radial preparation.** Piece B would need a preparation that fixes the source phase ψ2 and the Σ
modulus, which is a light inflationary mode in this model. It would have to hold:
- at k_*: to ≲ 1e-6 rad and ≲ 4e-7 in ln|Σ| per ln k in the linear regime. A preparation mass of m ≳ 1.5H for both
  would suppress the k_* amplitude by about e^{−1.5 N_*} ≈ 6e-7 relative to horizon exit [E];
- at the release Hubble scale: to a patch noise of ≲ 1e-3–1e-2 (B2 table).

The model contains no such interaction, and Codex states this openly.

## 4. Verdict table

| # | Claim checked | Codex or prior value | This track | Verdict |
|---|---|---|---|---|
| 1 | VECTOR_INITIALIZATION Higgs norm, m_A/H and BD normalization at the point (Codex) | F_g/H ≈ 34.83 (interpolated), m_A/H ≈ 3.483e-6, √P_α ≈ 0.00457 | 34.8322, 3.4832e-6, 0.004569 [C] | CONFIRMED |
| 2 | "Do not replace the evolving coupled modes by a constant-vacuum-mass relic" (Codex) | statement | Constant-mass GMR (0.856 pairs per net; 56.5 Y_Q-normalized) overstates the superhorizon BD result by (F_inf/F_vac)² = 4.2e5 [C] | CONFIRMED |
| 3 | BD longitudinal abundance with the evolving mass (not claimed by Codex) | — | N_A = 3.8e-4 to 1.7e-3 H_i³: 1.6e-6 to 7.1e-6 pairs per net (1.0e-4 to 4.7e-4 Y_Q-normalized), ≤ 1.4e-3 (≤ 0.09) of the allowance [C] | new; this piece does not fail the point |
| 4 | Our gp_critic/gp_leftovers: the longitudinal vector gives about 1060 pairs per net and is sufficient on its own at 4.7e13 GeV | 1053–1061 | For the inflationary BD source, 1.4e-4 at 4.7e13 GeV (frozen-angle rescaling) [C, D] | CORRECTED (our error: vacuum-mass GMR) |
| 5 | sp_yledger: the frozen angle gives about 2e-6 pairs; the ψ2 transfer was not computed | 2e-6 | 1.6e-6 (superhorizon modes) to 7.1e-6 (all modes) [C]; transfer now estimated (row 7) | CONFIRMED |
| 6 | "A homogeneous zero-Gauss solution does not supply a longitudinal occupation" (Codex) | statement | True for one patch. Patch-to-patch differences fix B_L = (k/e)[α + δθ_Σ(∞)/6] superhorizon [D] | CONFIRMED, extended |
| 7 | Release-induced Goldstone angle: Σ Higgs-charge winding of 2.3e3 rad (not in any ledger) | — | dθ_Σ(∞)/dψ = 1.3–1.7e6 rad/rad. B2 gives 0.03–1.0 pairs per net with natural phase or radial noise (6–210 times the allowance); B1 at k_* is unresolved (central 49 times) [C, E] | new source; fails the point on estimates unless phase and \|Σ\| are prepared to ≲ 1e-3 per patch and ≲ 1e-6 at k_* |
| 8 | Vector decay (assumed A → X X̄ in gp_leftovers) | Γ = 3.98e-12 GeV | The only open channel; T_vis 1.3–1.7 TeV; one pair of 3e4 GeV quanta per vector [D] | CONFIRMED |
| 9 | Remedy e ≳ 0.029 (gp_critic suggestion) | 0.029 | The threshold 0.0287 is confirmed, but it removes only piece A, which is already negligible. The X-orbit Floquet exponent is 0.15 per unit a at e = 0.0287 (0.06–0.22 up to e = 1); staying resonance-free needs e ≲ 2.6e-3 [C] | CORRECTED (not a free remedy) |
| 10 | Parametric vector production at the stated e = 1e-7 | not discussed | μ = 1.6e-12, so none [C] | CONFIRMED safe (new) |
| 11 | Codex's supplied exact radial minimum | supplied | The \|Σ\| modulus is a light inflationary mode (m² = 1.11 H²). Its Hubble-scale noise alone, with the phase fixed, gives a net ε of −0.048 [C] | new: a preparation must also fix \|Σ\| |

## 5. Limits

- **Piece A** uses a single-Higgs-direction Proca reduction.
  - It is exact for the frozen superhorizon gauge angle (k ≤ 1e-4 reproduces (m_vac/m_inf)² to 0.2%).
  - Modes entering the horizon while X dominates the Higgs norm are counted as vectors, which is an upper bound.
  - The inflationary background is stationary (as Codex's initialization assumes), followed by instant release and
    radiation.
- **Piece B** rests on:
  - homogeneous k = 0 runs read as separate universes;
  - a conversion through the P0 transfer, valid for frozen superhorizon angles;
  - for B2, a white-noise model whose coefficients (a_s, c_w) are bracketed, not derived.
  - No finite-k coupled X/Y/Σ/A lattice was run. Fragmentation of the X condensate near a ≈ 600 (sp_yledger [Q]),
    the phase slips (min|Σ| ~ 1e-5 of vacuum) and string formation can change B. They are not included.
- The θ_Σ(∞) estimator (vacuum-|Σ| tail) has a 30% systematic in derivatives (1.31e6 at a_end = 1000 against
  1.70e6 at a_end = 3000). The ensembles use a_end = 1000 and rtol = 1e-10.
- The noise-averaged slope at ψ0 = π/2 is statistically unresolved: 128–256 members per point, chaotic members.
- **Floquet** covers transverse modes, the X part of the Higgs norm only, and the late periodic orbit. It does not
  include early-release resonance, longitudinal modes, the Σ detuning or backreaction.
- Not included anywhere: vector-production nonlinearities, the vector-mediated X–X̄ interactions, a thermal-plasma
  modification of the decay, or strings.
- Codex's code was not run. The comparisons use its published JSON and handoff numbers.

## 6. Credit

- Graham, Mardon & Rajendran, PRD 93 (2016) 103520 (longitudinal vector from inflation; action and GMR formula) [R].
- Goldstone equivalence: Cornwall–Levin–Tiktopoulos; Lee–Quigg–Thacker [R].
- Separate-universe (δN-type) reasoning: Starobinsky; Sasaki–Stewart; Wands–Malik–Lyth–Liddle [R].
- Stochastic spectator equilibrium: Starobinsky–Yokoyama [R].
- Floquet/parametric resonance in λφ⁴ backgrounds: Greene, Kofman, Linde & Starobinsky, PRD 56 (1997) 6175 [R];
  gauge-field preheating literature (for example Lozanov & Amin) [R].
- Affleck–Dine; Q-ball formation: Kusenko–Shaposhnikov, Enqvist–McDonald, Kasuya–Kawasaki [R].
- Codex 373d37e (point, VECTOR_INITIALIZATION.json) [Q]; sibling tracks sp_yledger and sp_point [Q where cited].

## 7. Reproduction

```
cd research/reality_framework/grut_formation_review_20261025/tracks/sp_vector/code
S=<scratch dir>; gcc -O3 -o $S/vrel vrel.c -lm
python3 vmodel.py ../results/background.json                                   # 4 s (uses ../../gp_matching/code)
$S/vrel 1.5707963267948966 1 1000 600 1e-11 1.1211167304108865e13 ../results/hist_nominal.txt   # 1.4 s, nominal + history
python3 f2profile.py ../results/hist_nominal.txt ../results/F2_profile_nominal.json
python3 proca_modes.py ../results/background.json ../results/proca_modes.json ../results/F2_profile_nominal.json  # 2 min
VREL=$S/vrel python3 phase_transfer.py ../results/phase_transfer.json          # ~25 min on 4 cores (2,070 releases)
VREL=$S/vrel python3 phase_scan_smoothed.py ../results/phase_scan_smoothed.json # ~15 min (640 releases)
VREL=$S/vrel python3 radial_transfer.py ../results/radial_transfer.json        # ~5 min (140 releases; vrel.c with xfac yfac args)
python3 floquet_vector.py ../results/floquet_vector.json && python3 floquet_band.py ../results/floquet_vector.json  # ~10 min
python3 vector_ledger.py ../results ../results/vector_ledger.json ../results/hist_nominal.txt   # 1 s
```

All numbers in this file come from our own code [C], except where marked [Q] (Codex 373d37e, sp_yledger, gp_critic).
No reproduction of Codex's code was attempted.
