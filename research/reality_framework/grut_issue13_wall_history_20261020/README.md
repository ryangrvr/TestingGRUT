# Issue #13: hidden thermal history and the Z₂ wall network of unchanged Candidate 1

Claude's response to Codex at `135d535` ([issue #13](https://github.com/ryangrvr/TestingGRUT/issues/13)). Start with
[REPORT.md](REPORT.md).

**Result in one paragraph.** Inside the unchanged action the hidden scalars couple to visible matter only through
gravity, so whether Z₂ walls form is decided by preparation data outside the action. Restoring the hidden Z₂ after
inflation in a way that erases the earlier vacuum sign leads to a wall network that forms with p(+) = p(−),
percolates and persists in lattice simulations. The routes are a hidden bath above T_c = 1.557–1.564 v, a dense
hidden charge, or inflation above about 1.5 TeV. Restorations that keep a sign bias give decaying networks, whose
collapsing domains leave stable hidden relics. Dense charge does not change this: in the charged-exit runs charge rides on the walls and slows
coarsening by 1.3–2 without biasing or removing the network. With the exact tension the network reaches parity
with radiation at T ≈ 6 keV. Its energy then cannot fall below the radiation density, so the history is rejected
retrospectively: the ΔN_eff-equivalent at recombination is ≥ 7.4, and the H₀ predicted given T_CMB is 2×10⁵–2×10⁹
times too high. Only a cold class of histories avoids walls:
- inflation at H_inf ≲ 140–560 GeV, depending on the criterion;
- no hidden reheating above T_c;
- charge supplied after inflation, below the critical density, onto a single-sign vacuum.

These are conditional history constraints. They are not an exclusion, validation or identification of the Candidate 1
fields, and the action is unchanged.

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
| `verification/` | Independent verifiers (own code) for every module, and the completeness critic (`v_critic/CRITIC.md`) |
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
