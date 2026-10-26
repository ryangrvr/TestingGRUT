# Attempts and corrections

Only MOMENT_CONTROLS.json, FOLLOWUP_CONTROLS.json and RELAXATION_CONTROLS.json are the current reported numerical controls. Old trial trajectories are not included as additional results.

- The inherited repeated-inversion integrator stalled with stiff trial states and Radau Jacobian warnings. It was stopped. The active-set susceptibility ODE replaces it for production. Its unsuccessful run supplies no abundance.
- Early cold trials encountered underflow/fugacity products and saturation derivatives. Logarithmic rates and an explicit finite gas-cap regulator resolved those numerical failures.
- Early T=0.03m and T=0.005m trials stopped before the later capture/closure evolution. Their values were not treated as final.
- The first hard contact-cascade implementation lost interpolation weight placed in the lower bin adjacent to the hard/soft cutoff. Its energy drift was about 4%. Depositing that part into the soft reservoir restores number, charge and energy accounting. No physical input was tuned. Only corrected runs are reported.
- Reducing the Bose series cap exposed inconsistent initialization near saturation: TH.net_charge truncated while species used quadrature. The resulting initial charge offset was about 0.02–0.07%. A new initial matching routine uses the same species densities as the evolution. Claude's original build remains unmodified.
- The original background silently clamped visible temperature at 10 keV. The grid now extends to 1e-12 GeV and rejects out-of-range evaluations. The long instant-object controls are recomputed; no old late value is used in the report.
- The initial conservation audit failed because a loose-tolerance mediator-only gas run exceeded the 1e-6 energy-work check slightly. Production was rerun at rtol=2e-9, giving maximum residual 2.07e-9 and charge drift 6.62e-8. NUMERICAL_BASELINE.json retains summaries of the superseded loose attempt after the initialization/background repairs.
- The earlier tentative tail-bound function was not validated and is not used. No numerical inventory floor, temperature stop or loss of kinetic closure is promoted to a final physical abundance.
- A final input audit replaced the reused e=1e-7 mediator pair count in the e=1e-6 hard control with that coupling's own count; all four hard controls were rerun. Moment production already used the appropriate mediator energy and was unchanged.
- One verification invocation read the hard-control output while that rerun was still writing its first case and stopped. Verification was rerun after all four controls completed; a complete-case assertion now reports this dependency explicitly.
- No Claude workflow was invoked, no AI subagents were used, and no independent critic or verifier is claimed. Local CPU worker processes evaluate independent numerical controls only.
