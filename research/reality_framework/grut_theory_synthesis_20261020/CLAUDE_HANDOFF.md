# Claude construction brief: test the derived theory without refitting

20 October 2026 · same Candidate 1 action · additive work on the isolated research branch

Read [THEORY.md](THEORY.md) for assumptions, proofs and limits. The prior 31-run packet is [here](../grut_candidate1_two_object_20261019/persistence/REPORT.md). This brief is a construction target, not a claim that its nonlinear or apparatus tests have run.

## The new target

Keep the full scalar model a=1,b=.1 as the reference. The theory predicts a dense charged phase, its surface cost, its linear modes, and an unstable coherent dilute preparation. Check those with the evolution solver before expanding the physical interpretation. Store deviations as results; do not change the action to match each observable.

| Target | Prediction in the stated dimensionless convention | Controlled implementation |
|---|---|---|
| Increasing-charge stationary branch in flat space | Fcore -> 1.77827941, Score -> 0, ncore -> 1.77827941, Omega -> .562341325 | Supply a charge sequence and solve the full radial equations. Track core, radius and boundary error. Include stability checks; finite Q399 is not the asymptotic target. |
| Planar surface cost | sigma approximately .680013489 on the resolved branch | Independent full two-field interface solve. Check the first integral, tail slopes, alternate initial guesses and domain dependence. Least tension is not proved by the existing solve. |
| Large-Q energy surface correction | Ebar approximately .562341325 Qbar + 2.240448603 Qbar^(2/3) | Check leading scaling across increasing charge with a controlled large-radius hierarchy. Avoid treating the rough asymptotic crossing near Qbar134 as a calculated critical charge. |
| Dense bulk waves | nu_plus/minus squared = k squared + 3 Omega squared +/- sqrt(9 Omega^4 + 4 Omega squared k squared); nu_S squared = k squared + F squared - a | Initialize the uniform S=0, F=F0, Omega=Omega0 solution in a periodic box. Perturb one Fourier mode at a time at small amplitude; include all coupled fields. Compare frequencies, then refine space and time. Sound approaches c/sqrt(3); no claim this is light in vacuum. |
| Mixed-branch growth | F=.5, S=sqrt(.75), Omega=sqrt(.775). Unstable cutoff k=.5953450; sampled peak k=.4035, growth .1121396; low-k growth/k -> .412082 | Use an exactly stationary coherent uniform background, not the earlier nonstationary Gaussian cloud. Seed several wavelengths around and outside the band, record early growth, then test nonlinear saturation. Periodic box and mode availability are supplied choices. |
| Negative binding control | For a*b>=1, Ebar>=abs(Qbar) for every finite-energy state of the leading flat scalar model | Change b only in a separately labeled family control. Verify full matched energy; do not interpret missing core identification as the energy theorem. This is an altered-parameter control, not one of the old unchanged-coupling evolutions. |

No full nonlinear wave/growth experiment in this table has been executed by this packet. The matrix spectrum and interface calculation are reproducible in derive_and_check.py; the 10 author checks pass.

## Resolve the earlier m=1 ambiguity

Differentiate the full continuum trajectory under a global transverse translation and rotation. For every field and canonical momentum B:

```
translation amplitude: -d_rho B
tilt amplitude: rho*d_z B - z*d_rho B
```

Both multiply cos(azimuth). Compare evolved tangents with spatially differentiated background trajectories, using an interior norm and a declared boundary margin. Test convergence of the residual, not constancy of the chosen tangent norm. Rotation must act on momenta and velocity directions as well as positions. The current opposite-shift seed is not this full generator.

After these controls, align or project neutral motions with an explicit metric and report both raw and shape-response results. Retain parameter sensitivity. Then run finite-amplitude off-axis perturbations with backreaction in unrestricted 3D. Neither a clean symmetry control nor a homogeneous bulk spectrum proves nonlinear collision stability.

## Compose the capabilities

Only after the operating regime is established, build a trajectory that forms charge concentrations from the declared resources, brings them into an encounter through modeled preparation or internal controls, retains an acquired distinction, and lets a compatible reporter/controller act on that distinction. Include a missed encounter and a matched intervention that disables or scrambles the record without silently changing the available energy and access. An externally scripted force or freely added detector coupling would need an explicit enlarged model.

An internal control loop would show causal use of a record. It would not automatically establish subjective experience. Maintain the distinction between classical organization, a quantum state/measurement calculation, and calibrated tests against nature.

## Preserve the common model

The action, leading-sector approximation and benchmark are inherited from the original Candidate 1 packet. Gravity, quantum rules, 4D spacetime and Standard Model fields are adopted inputs. Full-world closure does not derive them. Do not import the earlier moving-channel or M1 memory theorems into this scalar system without a compatible construction. Do not use a higher-loop result to certify all-channel lifetime, or a small linear tangent to certify finite-amplitude 3D stability.

The decisive deliverable is a held-fixed prediction that survives a new controlled calculation or fails clearly, with the actual law, preparation, calibration and measured observable retained.
