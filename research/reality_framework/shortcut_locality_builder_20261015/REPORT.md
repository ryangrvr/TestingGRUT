# GRUT builder: shortcut-locality control and CP-12 acceptance criterion

Parent: `69098c127e33ad545044e524fe57cf80339cca08`. This is a new **independent diagnostic**; it does not change the frozen CP-12 Hamiltonian or its original neutral ray.

## Executed control

`shortcut_control.py` evolves a single particle on a ten-site cycle, then on the same cycle with chords (0,5) and (2,7). It uses exact dense finite-dimensional matrix exponentials and a 1001-point trapezoidal evaluation of the Duhamel integrand. The integral evaluation is numerical, not an outward-rounded bound. No CP-12 fields or dynamical graph states are included.

| time | ring-only receiver-5 probability | full shortcut receiver-5 probability | state distance | numerical Duhamel integral |
|---:|---:|---:|---:|---:|
| 0.10 | 2.76853283454e-14 | 0.00983448990829 | 0.0995430599709 | 0.0996674169020 |
| 0.25 | 2.59443646971e-10 | 0.0562662269806 | 0.242974199589 | 0.244865163965 |
| 0.50 | 2.49515695201e-7 | 0.162306681581 | 0.447002484278 | 0.460716411725 |
| 1.00 | 0.000198225548359 | 0.136672521924 | 0.672349678283 | 0.749971435790 |

All four independently computed state-distance comparisons satisfy the numerically integrated Duhamel inequality. This is a control example showing that retaining a core does not guarantee transport according to core distance.

## Exact inequality for the full candidate

For a fixed self-adjoint decomposition `H=H_core+V_shortcut` with bounded V, a common normalized initial state psi, and time t>=0:

```
||exp(-iHt)psi-exp(-iH_core t)psi||
<= integral_0^t ||V_shortcut exp(-iH_core s)psi|| ds.
```

For any norm-one projector O, the corresponding probability difference is at most twice this state-norm error. The bound applies to the same state and same Hilbert space; it does not justify comparing separately prepared or separately measured geometries.

For CP-12, a graph-dependent 'core' split must itself be specified as a Hermitian operator on the full gauge-invariant Hilbert space, including what happens when the core changes. One cannot drop off-core hopping *after conditioning on a measured core* and claim the unmeasured state is approximated. Gauge-invariant hopping along non-core active edges and their exact state-dependent action must be retained in V; the electric outgoing-budget residual must be added separately.

## CP-12 next computation (not executed here)

1. On each retained N5/N6 physical orbit, compute the corrected opposite-sign charged circulation core, including charged holes with unequal opposite-sign fields. Do not assume E_ie=-E_je.
2. Construct the gauge-covariant hopping operator on active edges **not** in that core. Explicitly verify the core/off-core split preserves Hermiticity, Gauss law and permutation symmetry. Do not drop graph rewiring or silently project out non-core states.
3. Propagate the same unfiltered neutral ray under the full frozen Hamiltonian and the controlled comparison operator. Bound both the state-difference Duhamel residual and all outgoing electric-space residuals. The N6 time-2 gate is not yet controlled.
4. For a chosen charged event and diagnostic graph record, report the unnormalized joint response and initial-event probability. Select fixed source/receiver labels independently of the final graph. Compare core-distance bins only in a joint recorded protocol and report the record energy/entropy change; no autonomous record is derived.
5. Measure whether full response tracks core distance or full active-graph distance, with random/irregular and degree-matched controls. N5 K5 and N6 Q12 hosts have saturated graph diameters; meaningful distance scaling needs larger or more constrained hosts.
6. Accept locality only when a nontrivial response-distance hierarchy and a small state-dependent shortcut error coexist over a time window, with a positive formation probability and full-law error control. Stable core labels alone are not sufficient.

**Scientific verdict:** one-particle control passes; it is not evidence of CP-12 shortcut suppression or failure. The next step is an actual gauge-invariant CP-12 operator/action computation.
