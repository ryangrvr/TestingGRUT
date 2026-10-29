# M1 energy-resolved kinetic continuation

28 October 2026. Input commit bf8ab3b96c07691410145c3b8c21ca18b009d66d.
Ryan asked this builder to continue the stalled transport calculation. The
shared conversation ends at the committed scattering packet; the additional
transport progress is supplied as prose. No numerical transport checkpoint was
recovered. This packet is a new implementation, not a replay of that unseen run.

Keep the M1 action, couplings, reference profile, measured pulse source, yield,
and existing 1% conditional cold-mass tolerance unchanged. No AI subagents or
external reviewer are used. Do not equate same-author numerical checks with an
independent-author review.

## Bounded model, fixed before running

Use an isotropic, spatially homogeneous, dilute gas per supplied reference core.
Start with the measured 10% radial pulse spectrum. Reproduce its archived
history before recovering positive, negative and neutral outgoing energy spectra.
The initial core charge is Q0-Nplus+Nminus. Initially emitted neutral quanta
decay once with the correct boosted two-body daughter distribution.

Hold the scattering profile and chemical potential at Q0=200000. Evolve core
charge through the charge ledger and core energy through its tangent with slope
Omega. This is a fixed-profile kinetic diagnostic; finite-Q backreaction is not
silently certified. Report charge excursion and compare the tangent cold-mass
diagnostic with re-solved stationary masses at the resulting charge.

Include both directions Xbar+core <-> X+core and neutral production from either
charged incident species, followed by prompt h decay. Eliminate real h only in
the explicitly dilute, fast-decay limit. Include the leading narrow-width
X Xbar -> h -> X Xbar redistribution once, with no separate inverse decay added.
Compare it with a core-only control and with ordinary XX and Xbar Xbar elastic
collisions computed from the same quartic and neutral exchange vertices.
Nonresonant opposite-charge interference, higher multiplicity, medium/Bose
effects, finite core recoil, formation and cosmological expansion are outside
this bounded diagnostic. Quantify the fast-decay and expansion assumptions at
the inherited illustrative 20 MeV clock; do not identify that clock with a
derived formation epoch.

Use particle counts on an energy grid, conservative linear deposition, exact
two-body energy sums, reciprocal partial-wave kernels and the full heat ledger.
No momentum clipping, hidden heat sink, artificial bath, coupling retuning or
renaming a rate as final abundance. Any unsupported high-energy transitions
must be counted as omitted rather than projected back onto the grid.

## Checks and decision boundary

Before interpreting population changes, verify inherited source reproduction,
scattering normalization and reciprocity, per-event charge/energy identities,
resonant rate normalization against the decay width, elastic conservation,
positivity and ODE tolerance. Compare energy-grid refinement and a larger energy
range. Test a Maxwell-Boltzmann distribution as an independent equilibrium
diagnostic for gas-only collisions; conservation alone does not establish the
correct collision operator.

Save distributions and ledgers at optical times tau=n_core*t/m_X^2 =
0, 0.1, 1, 10, 100 and 1000. A longer exploratory continuation must be labelled.
Report antiparticle count, positive-particle count, core charge, gas energy,
kinetic energy, fraction above reverse threshold and stationary cold-mass
diagnostic. A reduced antiparticle count is not a claim of irreversible cooling
or cosmological cleanup. Preserve numerical failures and resolution limitations.
