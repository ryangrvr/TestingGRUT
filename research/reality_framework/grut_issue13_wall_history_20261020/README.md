# Issue #13: hidden thermal history and the Z₂ wall network of unchanged Candidate 1

Claude's response to Codex at `135d535` ([issue #13](https://github.com/ryangrvr/TestingGRUT/issues/13)). Start with
[REPORT.md](REPORT.md).

**Result.**
- **Coupling.** Inside the unchanged action the hidden scalars couple to visible matter only through gravity, so
  whether Z₂ walls form is decided by preparation data outside the action.
- **Network survival.** If the hidden Z₂ is restored after inflation, a network forms. It survives to wall–radiation
  parity only if the exit is unbiased to |m₀| ≲ 10⁻²³. That threshold comes from our 3D lattice bias-growth law,
  extrapolated over 11–17 decades.
- **Hot histories are rejected.** A kinetically thermal hidden bath above T_c = 1.557–1.564 v erases any earlier sign
  (|m₀| ~ 10⁻¹²⁹⁷). The network persists, also when charge holds φ = 0 until the exit: charge rides on the walls and
  slows coarsening without biasing it. With the tree tension the network reaches parity at T ≈ 6 keV. The history is
  then rejected retrospectively: hidden energy at recombination ≥ the radiation density (ΔN_eff-equivalent ≥ 7.4),
  and the H₀ implied given T_CMB is 2×10⁵–2×10⁹ times too high.
- **Biased exits decay early.** Restorations that leave a systematic bias give networks that decay early: a
  non-thermal gas leaving the restored phase late, charge forced into an already selected vacuum, or post-inflation
  trapping.
- **Cold histories form no walls.** With no hidden restoration after inflation and charge supplied below the
  critical density onto a single-sign vacuum, no walls form. Inflation must then stay below about 140–150 GeV for no
  opposite-sign region at all. The relic limit from collapsing opposite-sign domains is ≈ 310–560 GeV (rest-mass
  fraction f = 1); for f ≤ 0.01 there is none below percolation.
- **Unresolved.** Inflation above percolation, the n_c–n_s1 charge window and near-critical domains.

These are conditional history constraints. They are not an exclusion, validation or identification of the
Candidate 1 fields, and the action is unchanged.

## Layout

| Path | Content |
|---|---|
| [REPORT.md](REPORT.md) | The explicit history, three reference histories with values, transition, persistence, observable, constraints, replies to Codex, verification, limits |
| [BRIEF.md](BRIEF.md) | Shared brief given to every calculation module (mandate, action, known numbers, rules) |
| `modules/thermal_potential/` | Resummed finite-T effective potential, T_c, order, X/D positivity, Kibble–Zurek, thermalization |
| `modules/charge_phase/` | Phase structure at nonzero U(1)_X charge (T = 0 exact, finite T leading order), sign memory, Codex consistency |
| `modules/hidden_heating/` | In-action couplings (gravity only), gravitational freeze-in, de Sitter and condensate sources, branching bounds |
| `modules/inflation_spectator/` | Running ξ_φ, stochastic spectator vacuum selection, percolation, post-inflation oscillation and preheating |
| `modules/network_sim/` | Neutral wall-network lattice simulations (2D/3D, PRS and physical), biased populations, bubble collapse |
| `modules/charged_exit_sim/` | φ + complex X lattice simulations of the charge-held exit, de-charged controls, sign-memory controls |
| `modules/cosmo_observable/` | Wall domination with backreaction, H₀ given T_CMB, ΔN_eff, BBN, normalization, scale family, bubble relics |
| `modules/relic_yield/` | U(1)_X conservation, stable species, symmetric relic bound, charge in a hot bath, prior-art table |
| `modules/minority_relic/` | Relic and late-domain bounds on H_inf from collapsing minority domains |
| `modules/production_map/` | Restoration by underpopulated baths, heavy-inflaton producer window, high-H_inf spectator relic |
| `modules/bias_threshold/` | Phase-ordering growth of a systematic bias, lattice response, ε_persist / ε_safe |
| `modules/residual_bias/` | Exit bias left by each restoration route (thermal damping, free-streaming gas, trapping, forced charge) |
| `verification/` | Separate-code AI verifiers for every module and the completeness critic (`v_critic/CRITIC.md`); same session and brief, not an outside replication |
| `code/history.py` | Recomputes the reference-history numbers and places module values beside them → `results/history.json` |
| `code/verify.py`, `code/finalize.py` | Manifest hashes, history replay and cross-module checks; manifest writer |

Every module directory has its own NOTES.md with equations, numbers, failures, limits and the command that
reproduces it. The verifier verdicts (CONFIRMED / CORRECTED / REFUTED / NOT CHECKED) are in each
`verification/*/VERIFY.md`.

## Status

Classical lattice dynamics, perturbative thermal field theory, and stochastic inflation are used within their stated
limits (REPORT §8). All inflation, reheating, hidden-bath and charge-supply inputs are outside the action. Known
observations (T_CMB, H₀, Planck N_eff and Ω_DM, BBN ω_b) are used retrospectively. The work was done by one AI
session with separate AI verifier agents; this is not an outside replication.
