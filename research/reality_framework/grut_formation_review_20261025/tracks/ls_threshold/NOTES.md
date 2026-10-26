# Track ls_threshold: case-L threshold and final f for low-scale candidates 10 and 17 (Codex b448f25)

Claude, 26 October 2026. This is same-session AI work and is not outside replication. It reviews one benchmark
inside one supplied tree-level model. It is not a validation, an exclusion or a GRUT prediction. Codex's packet
was read from `git archive b448f25` in the scratchpad. Codex's code was read but **not run**.

## Verdicts

Inputs [Q], from b448f25 DECAY_RESULTS.json and VECTOR_RESULTS.json:

| | Candidate 10 | Candidate 17 |
|---|---|---|
| ε_net | 0.157413 | 0.147046 |
| Raw yield | 5.6901 Y_Q | 5.3094 Y_Q |
| R_Y | 113.07 | 122.81 |
| R_lower | 199.40 | 249.69 |
| R_central (including primordial vectors) | 228.42 | 280.53 |

"Own" means this track's independent code (`lsrelic.py`, proxy mode). "Earlier" means our sp_yledger
`thermal_relic.run_inj`, run unchanged.

| Codex claim (10 / 17) | Own [C] | Earlier code [C] | Verdict |
|---|---|---|---|
| R_L at calibrated density: 91.17–91.25 / 83.12–83.18 | 91.17–91.25 / 83.13–83.20 | 90.0–90.6 / 81.9–82.5 | CONFIRMED |
| R_L at raw density: 102.96–103.04 / 93.39–93.47 | 102.97–103.05 / 93.44–93.52 | 100.0–100.6 / 90.6–91.25 | CONFIRMED |
| Conservative f, raw: 10.85 / 13.94 | 10.849 / 13.942 (identical at N = 800 and 1600) | 10.849 / 13.942 | CONFIRMED |
| Central f, raw: 12.08 / 15.27 | 12.080 / 15.265 | — | CONFIRMED |
| Y alone, raw: 0.570 / 5.750 | 0.5696 / 5.749 | 2.18 / 7.93 (endpoint at 10⁻⁶ GeV) | CONFIRMED (endpoint-sensitive; fails either way) |
| f ≥ 9.90 / 12.80 after the optimistic bound on later conversion and 4→2 | 9.90 / 12.80 (tail depth 0.0916 / 0.0854) | — | CONFIRMED |
| Calibrated density, conservative / central: 11.31, 12.59 / 14.51, 15.88 | 11.312, 12.592 / 14.505, 15.880 | — | CONFIRMED |

Our earlier code puts R_L 1–3% lower than Codex and this track. It stops at T_vis = 10⁻⁶ GeV rather than 10⁻⁸ GeV,
so it misses late lump erasure, and it uses an explicit endpoint update. Near the threshold, f is
endpoint-sensitive: this is why the Y-only values differ. Away from the threshold the codes agree to every
quoted digit. Codex had already replayed run_inj at the candidate-10 R_lower point; this track adds the
thresholds.

**The benchmark failure is robust** (details below). The three alternative treatments in this track give, for
the conservative ledger at raw density:

| Treatment | f, candidate 10 | f, candidate 17 |
|---|---:|---:|
| Two-potential energy-conserving state (ke) | 6.45–7.84 | 8.41–10.15 |
| ke plus the on-shell s route | 3.73 | 4.97 |
| Benchmark plus the on-shell s route | 5.40 | 7.07 |

Every one is at least 370 times the allowance f = 0.01.

## Method

- **`code/lsthermo.py`** (own) computes Bose–Einstein thermodynamics from the modified-Bessel series (Kolb &
  Turner). It is a different method from Codex's Gauss panels and from gp_leftovers' adaptive quad. It agrees with
  gp_leftovers `equilibrium.thermo` to ≤3×10⁻¹⁴ from T/m = 0.003 to 30.
- **`code/lsrelic.py`, proxy mode** (own) is the benchmark as specified: the single-μ isentropic state with lumps
  at μ_c, and p relaxing with the quoted Γ. Unlike Codex's midpoint exponential update, the isentrope is tabulated
  and d ln p/d ln a is integrated with scipy Radau (rtol 10⁻⁹). Thresholds are bisected 11 times at N = 400, and
  the bracket edges are re-checked at N = 800.
- **`ke` mode** (own premise probe) keeps the actual X̄ number in the thermodynamic state: two chemical potentials,
  kinetic equilibrium, and lumps or anti-lumps when μ± reaches μ_c. Every reaction conserves energy, so
  annihilation away from equilibrium produces entropy. Rates carry the factor 1 − e^{−(μ₊+μ₋)/T}. Each step is an
  adiabatic expansion followed by a backward-Euler reaction at fixed ρ.
  - N = 600 and N = 1200 differ by 0.4% (6.450 against 6.423).
  - The energy-work residual is 1.5% cumulative at N = 600 and 0.38% at N = 1200, scaling as step².
- **The s route** is X X̄ ↔ s (Γ_s→XX̄ = 1.73×10⁻⁵ GeV, which is 2×10¹⁰ H at onset), followed by
  s s → X X̄. It uses σv = 1.26×10⁻⁴ GeV⁻², from our gp_leftovers Berends–Giele amplitude, |M| = 0.026667
  (analytic check exact). The density n_s = K(T) n₊n₋ uses Maxwell–Boltzmann statistics.
- **Our earlier code** is driven by `code/ours_threshold.py`, which imports sp_yledger `thermal_relic.run_inj`
  unchanged. rmax_relic.py is the same run_inj on a grid; it is superseded here by bisection.

## Premises probed

1. **Threshold definition.**
   - f(R) is monotone. At raw density for candidate 10: f = 10⁻¹⁹⁶ at R = 60, 10⁻³ at 100, 0.05 at 106, 0.57 at
     113, 4.6 at 130, 10.85 at 199, 18.7 at 400. There is no re-entrant window.
   - For the failing ledgers, f is frozen at the 4→2 freeze-out (T_vis ≈ 0.8 GeV, T_h ≈ 0.4 m, p ≈ 5.35). The
     endpoint and the late-tail bound are irrelevant there.
2. **T_h/m.**
   - At onset T_h = 12.0–13.5 m ≈ 1.2 v.
   - The one-loop T_c is 18 m with φ in the bath, or 42 m with only X. m_X(T)/m is 0.85 or 0.98 at onset; these
     corrections fall as T² and are negligible by freeze-out [E].
   - Extra light dof in the bath (s, and the stable D with m_D = 1.005 m_X, omitted in both codes) raise σ by up to
     19%. That is equivalent to R × 1.26, so it lowers R_L.
   - The vectors decay at T_vis ≈ 1.3 TeV, well before the 26–46 GeV onset, so placing R at onset is consistent [C].
3. **Processes.**
   - 3→2 among X alone is forbidden by charge.
   - Threshold 4→2 rates scaled ×10 or ×0.1 move R_L only between 103.3 and 102.8.
   - Conversion σv ×0.1 gives R_L = 80.2 and ×72 (the "daughter" value) gives 148.2 / 134.7, still below 199 / 250.
   - Benchmark evolution to 10⁻⁹ GeV with σv ×100 still gives f = 10.06.
   - Capture is instantaneous at μ_c, which is optimistic.
4. **Omissions in both implementations** (CORRECTIONS to the benchmark's premises, not to its arithmetic).
   - (a) *Cannibal state after freeze-out.* The proxy evaluates T_h, μ and lump formation on the
     chemical-equilibrium isentrope after 4→2 has frozen. It also lets p relax isentropically, discarding the 2 m_X
     of heat per annihilated pair. With the frozen X̄ number (ke), lumps re-form at T_vis ≈ 0.2 GeV instead of
     ≈10⁻⁸ GeV, and erasure on lumps becomes **heat-limited**: rescaling conversion σv ×100 leaves f unchanged.
     - Conservative f: 6.45 (μ_c = E_c = ω) to 7.84 (E/Q) for candidate 10; 8.41 to 10.15 for candidate 17.
     - Y alone: 3.93 / 4.60.
     - The ke crossing lies at R ≈ 10, not 103: f = 6×10⁻²¹² at R = 3, 0.013 at 10, 0.93 at 30, 3.5 at 100. So the
       benchmark's R_L is not conservative. The R = 0 and R = 3 runs hit about 200 reaction-root fallbacks with p
       near the 10⁻²⁵⁰ floor, so they mean only "erased".
   - (b) *Lump-phase bookkeeping.* μ_c = ω = 0.861 m but E_c = E/Q = 0.923 m. In ke this creates (E_c − μ_c) per
     condensed charge, a 16% energy-work violation. Fixed by E_c = μ_c. The proxy R_L moves from 103.0 to 99.9 for
     μ_c = 0.923.
   - (c) *On-shell s.* m_s = 2.45 m_X > 2 m_X. The threshold 4→2 coefficients evaluate the s propagator at
     p² = 4 m², but plasma pairs reach p² = m_s². The resonant route is up to 1.4×10⁷ times faster at T_h ≈ 0.3 m
     and dies below T_h ≈ 0.03 m. It roughly halves f.
5. **Can anything plausible erase ×1000?** Erasure must remove about 2 m_X per pair, which is 11–14 m_X per net
   charge, from the hidden sector.
   - In-sector channels are heat-limited, whatever their rate.
   - A two-body sink at T_vis ≈ 1.3×10⁻⁸ GeV, where the benchmark p_eq reaches 0.005, needs σv ≳ 4×10⁻⁴ GeV⁻²
     (5×10⁻²¹ cm³/s). The portal gives 2×10⁻³⁴ GeV⁻². The U(1) vector weighs 6×10⁴ GeV. The s → SM width is
     Boltzmann-killed.
   - Larger, more bound lumps help only weakly: ke with μ_c = 0.3 m still gives f ≈ 1 (indicative, inconsistent
     E_c).
   - **No modelling choice probed rescues either candidate. The verdict does not hinge on one choice.** The ×1000
     margin itself does: it becomes ×370–1000 across the treatments here.

## Limits and open items

- ke is still an equilibrium proxy, with instant condensation and hidden kinetic equilibrium. Elastic Γ/H reaches 1
  near T_vis ~ 10⁻⁷ GeV. It is not a Boltzmann, lattice or solitosynthesis calculation.
- The s route uses MB K(T) and the threshold σv. The D relic is not computed: it can only add.
- Our earlier case-L thresholds (sp_critic, R_L ≈ 20–175 for ε = 0.07–0.24) used the same isentropic proxy. They
  should be revisited with energy conservation.
- `results/ke_probe.json`, `ke_scan.json` and `ke_low.json` used the inconsistent E_c ≠ μ_c lump energy. They are
  superseded by `ke_consistent.json`.

## Reproduce (from `code/`, OMP_NUM_THREADS=1)

`python3 -B run_proxy.py ../results/proxy_benchmark.json` (about 30 min on 4 cores);
`python3 -B ours_threshold.py ../results/ours_threshold.json` (about 1 h);
`python3 -B run_sens.py ../results/sensitivity.json`; `python3 -B run_ke_consistent.py ../results/ke_consistent.json`
(about 40 min); `python3 -B premises.py ../results/premises.json`; `python3 -B summary.py`.

Credits: Bose integrals (Kolb & Turner); cannibal freeze-out (Carlson–Machacek–Hall 1992); Q-ball formation,
thermodynamics and solitosynthesis (Kusenko–Shaposhnikov; Enqvist–McDonald; Laine–Shaposhnikov; Griest–Kolb;
Postma); asymmetric dark-matter symmetric-remnant annihilation (Graesser–Shoemaker–Vecchi 2011). Benchmark design
from gp_leftovers/sp_yledger (ours) and b448f25 (Codex).
