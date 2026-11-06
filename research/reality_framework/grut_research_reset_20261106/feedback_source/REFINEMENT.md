# Registered refinement after the first force miss

Keep the [first protocol](PROTOCOL.md) and its failed result. At the
fine grid the second-order split shifts the finite static linear
susceptibility by 0.0123084. Correcting only the sampled linear kernel
reduces the difference from complete linear replacement to about 1.5%
of the corrected peak. The coarse momentum-spacing shift is still about
42% of that peak. It occurs late in the drive. Resolve both errors before
authorizing the live pilot.

Use a symmetric fourth-order composition of three second-order
mass/free steps. Measure the force at the full-step endpoints.
Linearize the identical composition to make the complete-linear
replacement diagnostic. Keep the invariant matching convention and
the two quadratic scalar terms from the first protocol. Use the same
closed pulse and in-vacuum. Fold the Bloch integral onto positive
quasimomentum by spatial reflection and spin rotation. Verify the
fold against the full quadrature on a small grid.

| Case | Bands to each side | Bloch nodes before folding | Radial nodes per panel | Radial cutoff | Timestep | Boundary |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Base | 7 | 16 | 8 | 8 | 0.04 | Open |
| Time | 7 | 16 | 8 | 8 | 0.02 | Open |
| Spacing | 7 | 20 | 10 | 8 | 0.04 | Open |
| Longitudinal | 9 | 16 | 8 | 8 | 0.04 | Open |
| Transverse | 7 | 16 | 8 | 6 | 0.04 | Open |
| Fine | 9 | 20 | 10 | 8 | 0.04 | Open |
| Tapered | 9 | 20 | 10 | 8 | 0.04 | Tapered |

Keep the same 3% force tolerance on each measured shift, base-to-fine,
fine-to-tapered and fine matching-method difference. The radial
comparison here holds the finer cutoff in the base. Keep 0.2%
net-work closure. On a small independent pulse, require the full vs
folded result to agree to \(10^{-9}\) relative and the fourth-order
timestep error to fall by at least a factor of eight on two halvings.

If this gate passes, the original live holdout and its acceptance
criteria remain. Use the corresponding fourth-order feedback
integrator and the new base, time, fine and tapered grids at the
specified holdout wavenumber. Use the folded quadrature and preserve
the signed substeps in the drive-work integral. Sample the independent
Hamiltonian ledger at intervals of 0.08 and at the endpoint. If this
refinement misses, preserve the miss and decide which physical
track has better observational access.
