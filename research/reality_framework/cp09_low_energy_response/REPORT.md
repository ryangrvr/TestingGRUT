# CP09 low-energy propagation diagnostic

2026-10-11. Frozen CP09 (1dadbd6) unchanged; this is a response/spectrum test, not CP10 or a Stage3 card. No parameter change, selected geometry or fitted kernel.

## What was computed

Use the exact same Hamiltonian with J=1, Omega=1/4, gamma=pi/4. Compute its lowest 12 energies for all labelled orders at N=2,3,4, and its stationary ground state. Declare source B=X_0 and detector A=Z_1 independently of the geometry distribution. Calculate R_AB(t)=i<[A(t),B]> and detector connected symmetric noise from the full matrix-exponential action, NOT a truncated 12-mode response. Compare R with an independent impulse exp(i epsilon B), epsilon=1e-5. Noise and response use the same state and H.

| N | joint dimension | E0 | gap E1-E0 | cross response t=.4 |
|---|---:|---:|---:|---:|
| 2 | 12 | .09475528 | .31829854 | .01063039 |
| 3 | 152 | .31097923 | .18571221 | .00863717 |
| 4 | 3504 | .60375765 | .09650386 | .00690217 |

Eigenvector residuals <=2.5e-11. Impulse and commutator responses agree within 1e-8. Dense full spectrum at N=2, sparse lowest modes at N=3,4; all configurations included. The small sample shows decreasing gaps, but it does NOT establish gap closure, an exponent, stable dispersion or any large-N scaling law. No finite-size fit is credited.

The same ground geometry distributions give mean heights 1.66416,2.25483,2.66993. N=4 height-two weight .41926 and height-three weight .47797. These are not continuum classes, and cannot establish manifold selection.

## An exact operational limitation of the symmetric state

For a nondegenerate ground state, event-permutation covariance of the supplied H makes its ground-state density matrix invariant under simultaneous permutation of geometry and matter. Therefore the same specified probes X_i,Z_j have identical R_ji(t) for every distinct ordered pair i,j. The numerical lowest spectra are consistent with a nondegenerate ground state at these sizes; nondegeneracy at arbitrary N is not proved. More generally this equality holds for any permutation-invariant prepared density matrix, without requiring nondegeneracy.

This is not a universal no-go for emergent locality: positions themselves can be relational, and localized observables may be defined using quantum geometry. It is a precise warning that site labels in the symmetric state cannot provide an operational separation-dependent propagation profile. Choosing a detector's geometry-dependent support after inspecting response would be circular. One must specify physical reference records or a relational localization protocol first, then calculate its conditional response and noise.

The global labelled probes used here are mathematically defined tensor-factor operations. Their realizability as localized interventions in a putative emergent world is NOT established. A nonzero cross response consequently is not proof of superluminality or of a physical wave.

## Physical verdict

CP09 has healthy finite low-energy excitations and genuine stationary matter cross-response. Positivity has not produced a physical wave equation. There is no reconstructed momentum, propagation speed, relational clock, Lorentz cone or continuum spectral limit. Decreasing finite gaps could be useful, but calling them a massless particle would be unsupported.

The next real obligation for THIS candidate is independently defined relational localization and clocks, accompanied by a controlled larger-N approximation whose errors can be bounded. The same approximation must retain the geometry distribution, not replace it with a selected graph. This diagnostic changes no kinetic principle and launches no hidden replacement law.

## Reproduce

python research/reality_framework/cp09_low_energy_response/check.py, NumPy/SciPy. Expected CP09_LOW_ENERGY_RESPONSE_PASS. The script loads only the frozen CP09 builder definitions before its top-level tests; it does not alter that file. Results JSON records every evaluated time, noise value and lowest energy. These are author numerical and symmetry deductions, not independent scientific review or empirical evidence.
