# Finite-charge branch and explicit fragmentation checks

The same tree-level flat-space action and physical couplings were used throughout: a=1, b=0.1, g=0.1, lambda=0.01, lambda_X=0.001. No architecture or physical parameter was searched or adjusted. The reference profile is read from the preceding Candidate 1 directory without modifying it.

## Reproduce

Run in this directory, or give the full path to the script:

```sh
python branch_scan.py
python branch_scan.py --refine
```

Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. No packages installed. There were 112 successful scalar BVP attempts, including refinements and alternative initial guesses, plus six successful linear sensitivity BVP solves. Scalar BVP time recorded by the script sums to 3.88 seconds in this environment; eigenvalue and bookkeeping work are additional. No background calculation remains running.

The first numerical protocol was written before the scan to `branch_protocol.json`. Its initial turning minimizer reached the declared 16-call cap; this failure to meet its optimizer stopping rule is preserved in `branch_first_pass_results.json` and the execution log. A separate bounded sensitivity/domain-refinement protocol was written before the second pass to `branch_refinement_protocol.json`. No scalar BVP failed. All four independent alternate seeds converged to the same found frequency family. Two distinct nodeless profiles at a given eligible fixed charge were resolved and retained; these are the lower- and higher-frequency portions of the connected family.

## Found family and thresholds

Continuation covered Omega=0.80,0.81,...,0.99,0.995,0.9975,0.999. Domains grow as the tail length increases. The reference is Qbar=399.4384624223, Ebar=363.9229985517, Ebar/Qbar=0.9110865197.

| Feature on the found family | Omega | Qbar | Ebar/Qbar |
|---|---:|---:|---:|
| First binding against free X carriers on lower-frequency portion | 0.90654404078 | 176.60438832 | 1.00000000000 |
| Minimum stationary charge / turning point | 0.95980474869 | 148.73404211 | 1.01389288576 |

The turning point was obtained from a root of dQbar/dOmega using the differentiated scalar equations. Qbar changed by 1.65e-9 between L=80,tol=1e-9 and L=100,tol=3e-10. Binding Qbar changed by 8.36e-10 between L=60,tol=1e-9 and L=80,tol=3e-10. These are convergence comparisons, not rigorous interval error bounds; root tolerances were 2e-9 and 2e-10 in Omega, respectively.

For g=0.1, physical charge is Q=100 Qbar. The found stationary-charge minimum is therefore about 14873.4, and the free-particle binding threshold is about 17660.4. A minimum along this found spherical family does not establish that no other localized states exist below it.

At the reference charge, the second resolved fixed-charge state has Omega=0.998480387629 and Ebar=400.0081599313, above both the reference energy and the free-carrier energy. It is retained in `targets["1.0"]`.

Fixed-Q real radial Hessian checks on either side of the turning point give:

| Omega | Lowest eigenvalue, N=800,L=70 | Lowest eigenvalue, N=1200,L=70 |
|---|---:|---:|
| 0.95 | +0.03044218 | +0.03052071 |
| 0.97 | -0.02740401 | -0.02733581 |

The positive rank-one charge term was included. These signs support a radial stability change. These additional branch checks cover ell=0 only; the original Candidate 1 angular/phase checks are not silently extended to every scanned state.

## Resolved energy comparisons

Total charge in every comparison is Qstar. Energy is dimensionless Ebar; actual energy differences are (v/g) times these values. Each tabulated two-lump value uses the less energetic of the two resolved fixed-charge profiles for each component. Every combination of resolved branch choices is separately retained in `all_resolved_two_lump_pairs` and `all_resolved_branch_combinations`.

| Components as fractions of Qstar | Total Ebar | Excess over reference Ebar |
|---|---:|---:|
| All free same-sign X carriers | 399.43846242 | +35.51546387 |
| Two Q-balls: 0.375 + 0.625 | 392.43733471 | +28.51433616 |
| Two Q-balls: 0.4 + 0.6 | 393.26201569 | +29.33901713 |
| Two Q-balls: 0.45 + 0.55 | 394.27341535 | +30.35041680 |
| Two Q-balls: 0.5 + 0.5 | 394.59684111 | +30.67384256 |
| Two Q-balls 0.375 each + 0.25 free | 403.47437881 | +39.55138026 |
| Two Q-balls 0.4 each + 0.2 free | 402.23944861 | +38.31645006 |
| Two Q-balls 0.45 each + 0.1 free | 398.84044373 | +34.91744518 |
| One Q-ball 0.8 + 0.2 free | 379.02781471 | +15.10481616 |

All other resolved lump-plus-free combinations are retained in `all_resolved_lump_plus_free`. The largest coarse-to-fine fragment energy change was 4.02e-7; the largest refined charge residual was 1.43e-9. These numerical errors are much smaller than the listed excess energies.

Equal three-, four-, and five-Q-ball targets lie below the minimum charge of the found family and were left unresolved, not declared physically impossible. In fact, three times the found minimum charge is 446.20212632 > Qstar. Thus partitions of three or more positive-charge lumps from this particular family cannot meet Qstar; this is a statement restricted to the found family. Multi-component channels with free carriers were explicitly resolved above. Some small resolved lumps have E/Q>1, so free carriers at the same charge cost less than those lumps; the free-carrier alternatives remain separately included.

## Conditional comparison for all partitions drawn from the found envelope

For fixed couplings and a smooth stationary branch,

    dE/dQ = Omega.

Writing K for the positive gradient energy, the flat-space virial relation gives

    E - Omega*Q = 2*K/3 > 0,
    d(E/Q)/dQ = (Omega*Q - E)/Q^2 = -2*K/(3*Q^2) < 0.

The identity E-Omega*Q=2K/3 was checked on the refined profiles, with maximum absolute residual 2.21e-9. The lower-frequency branch has decreasing E/Q as charge increases toward Qstar. The sampled higher-frequency portion has E/Q>1. Free same-sign carriers have energy per unit charge 1.

Consequently, **conditional on components being drawn from these continuous found branches or free same-sign carriers, and on no unobserved lower-energy state intervening**, each component with 0<q<=Qstar has E(q)/q >= E(Qstar)/Qstar. For any positive-charge partition with sum q_i plus free charge q_free equal to Qstar,

    sum E(q_i) + q_free
      >= [E(Qstar)/Qstar] * (sum q_i + q_free)
      = E(Qstar).

This argument is conditional on the branch envelope; a finite scan does not prove its completeness. There is no uniform positive gap for arbitrarily small free charge: a Qstar-dq lump plus free dq approaches the reference energy as dq approaches zero, with first-order excess (1-Omega_star)*dq.

## Verdict and scope

Every resolved channel costs more energy than the reference profile, and the found lower-charge branch supports the conditional comparison above. This strengthens the reference state's finite-charge classical stability evidence within the found spherical branch envelope. It does not establish a global minimum over all fields or all partitions, exclude nonspherical or additional branches, analyze opposite-charge pair creation, establish nonlinear dynamical stability, or prove quantum stability. The action, couplings, and tree-level flat-space approximation are unchanged.

Raw evidence: `branch_results.json`, `branch_states.csv`, `branch_profiles.npz`, `branch_thresholds_refined.json`, `branch_execution.jsonl`. Protocols and the first-pass archive accompany them. No original Candidate 1 or coordination files were modified.
