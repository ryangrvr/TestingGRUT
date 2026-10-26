# Review of Codex 34f23f8: kinetic continuation of our kinetic-threshold build

27 October 2026. Claude, issue #13. Codex continued our unfinished solver (`tracks/kin_solve`, pinned at e71730c) in
its own packet `research/reality_framework/grut_kinetic_continuation_20261026/`. This is a focused review: analytic
checks of the load-bearing formulas, a check of the bugs Codex fixed in our partial code, and a reproduction of four
key moment runs with Codex's code (reproduction, not independent code). No independent solver was written. Same-session
AI work; one supplied model; nothing observational.

## Verdicts

| Item | Verdict |
|---|---|
| Fixed denominator: f = (2m_X Y_Q / 4.4×10⁻¹⁰ GeV)·P/n_Q = 2.0286·P/n_Q | CONFIRMED (matches kin_source's pairs/net → f factor 2.0286) |
| Corrections to our partial code: scalar storage in the inventory ledger, omitted P-changing channels (E_XXbs_XXb and the B/E2/G2/J conjugates), log-combined rates at cold T, matched initial-state evaluation | ACCEPTED as bugs in our unfinished code |
| Visible-temperature interpolation silently froze H below 10 keV | CONFIRMED: our `ks_common.py` table starts at T = 10⁻⁵ GeV and `np.interp` clamps outside it |
| Thermal-state inequality e ≤ 1 + (m_s/m_X)p + 3t(1 + 2p) | CONFIRMED by derivation: rest energy ≤ m_X n_Q + m_s P (worst case puts the whole inventory in s, since m_s = √6 m_X > 2m_X); kinetic ≤ 3T per particle per Bose-series term for μ ≤ m; gas number ≤ n_Q + 2P. It needs E/Q ≤ m_X for objects, so it does not cover objects between the stationary floor (Q = 51,964, E/Q = 1.0146 m_X) and the binding crossing (Q = 65,620) |
| Object thresholds and fragment arithmetic (Q_min 51,964; binding 65,620; fastest fragment ≈ 2.5×10⁴) | CONFIRMED (our kin_thermo values) |
| Hard-source control: energy leaves the > 10T bins by T ≈ 2T_c in a supplied soft MB gas with contact scattering only | CONFIRMED as a control; agreed it does not prove thermalisation of the coherent release |
| Instant-object stopping values (0.946 / 1.195) are not final abundances; the 2.03×10⁻¹⁸ entries are the chosen floor | AGREED |
| Gas-only moment controls (reproduced, see below) | see table |

## Reproduction of the gas-only controls (Codex's code, our run)

Driver `tracks/kin_review/code/repro.py`, output `tracks/kin_review/results/repro.json` (Radau rtol 2×10⁻⁹, inputs from
our own kin_* tracks rather than Codex's extraction; about 30 s per run).

| Case | Total R | Codex f_stop | Our run | Charge / energy-work residual |
|---|---:|---:|---:|---|
| e = 10⁻⁷ band, gas-only | 41.707 | 1.49213 | 1.4921315 | −4.1×10⁻⁸ / −9.5×10⁻¹⁰ |
| e = 10⁻⁶ band, gas-only | 4.9665 | 0.00458554 | 0.0045855 | 4.0×10⁻⁸ / 2.0×10⁻⁹ |
| diagnostic R = 5, gas-only | 5 | 0.00476 | 0.0047584 | 1.2×10⁻⁸ / 6.1×10⁻¹⁰ |
| diagnostic R = 7, gas-only | 7 | 0.02884 | 0.0288444 | 4.1×10⁻⁸ / 1.7×10⁻⁹ |

All four reproduce to the quoted digits. Log-interpolating between R = 5 and 7 puts the gas-only crossing at
R ≈ 5.8. That number is a property of the recipe, not a kinetic threshold of the theory.

## What the controls say, read conservatively

- The gas-only recipe crosses the 1% allowance between total R = 5 (f = 0.00476) and R = 7 (f = 0.0288). That sits
  with our **energy-conserving proxy** (R_L ≈ 7–8, about 15–20 with the s route) and far from the isentropic proxy
  (≈ 82). It is the first kinetic (moments) evidence on the question, and it points the same way as both build-stage
  corrections (the old s route over-erased by up to 170×; the reservoir holds 1.5× more pairs).
- Under that recipe, e = 10⁻⁷ fails by about 150–180× (f = 1.49 band, 1.81 broad) and e = 10⁻⁶ passes
  (f = 0.0046). The vector gate therefore leans toward requiring e ≳ 10⁻⁶ (R_vector ∝ 1/e), consistent with the
  hw_critic statement "e ≳ 3×10⁻⁶ under the energy-conserving threshold".
- This is conditional on: a forced chemical and kinetic equilibrium start (the release spectrum is not matched); no
  object formation (gas-only); MB fugacities with zero-μ Bose factors rather than the full nonzero-μ Bose collision
  integral. Objects can only change the answer through formation kinetics, which nothing here derives.

## Still open (unchanged)

The physical f_sym and kinetic threshold remain undetermined. The two physical inputs that decide them are the
coherent release's actual distribution (matched, not forced thermal) and fragment coarsening into objects above the
binding crossing (merger kinetics, finite-Q energy). The 2cbc7f1 linear-mode packet is still awaiting review.
