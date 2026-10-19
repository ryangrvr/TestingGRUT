# Rotating scalar fluctuation spectrum: bounded numerical evidence

This pass uses the original four-scalar action and fixed matching, at the original reference frequency Omega=.8000000019763136. It calculates the tree-level flat-space rotating scalar generator, including the existing decoupled D fluctuation. It does not change the couplings, fit a new state, insert the homogeneous one-loop potential into the soliton, or calculate a quantum energy.

## Reproduce and frozen workload

From this directory:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python modes.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python modes.py --box-check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python modes.py --ell5-check
```

Each protocol was written before its calculations: `mode_protocol.json`, `mode_protocol_amendment_01.json`, and `mode_protocol_amendment_02.json`. The original results and execution record remain in `mode_initial_results.json` and `mode_initial_execution.jsonl`. The first amendment separates volume effects at exactly matched spacing and tightens one phase-zero test. The second fills ell=5 after the separately frozen numerical angular comparison indicated a conditional cutoff at ell=6. No unfavorable root was discarded as a physical result; symmetry-related numerical splits are retained and diagnosed explicitly.

Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0 were used, with one BLAS thread. The three workloads took 16.56, 6.49 and 3.13 seconds in this environment, respectively. No budget cutoff, stationary-solver failure or ARPACK nonconvergence occurred. No package was installed and no numerical job remains running.

The original profile is read without modification. Its SHA256 is recorded in the protocol. Dense complete finite-grid generator spectra cover ell=0,...,4 at (N,L)=(120,30),(180,30), and ell=5 at (180,30). Sparse checks cover ell=0,...,4 at (400,30),(800,30),(800,40), plus ell=0,...,5 at the exactly matched-spacing pair (599,30),(799,40), and ell=5 at (800,30). Each sparse calculation uses 18 eigenvalues at each of s=.08i,.18i,.05, tolerance1e-10 and maxiter1500. Summary duplicates within2e-7 are merged; the near-zero phase Jordan pair may therefore be represented once in a sparse summary. Dense spectra retain both roots. Returned raw shift values are also stored by the final source and in the ell=5 execution records.

## Operator, stationary discretization and charge

Time and radius are tau=g*v*t and rho=g*v*r. The perturbation is X=e^(-i Omega tau)[F+u+i w]v/sqrt(2), phi=v(S+eta). Transform q=(rho*u,rho*w,rho*eta), giving radial measure d rho. The finite-grid equations are

    qddot + G*qdot + K_ell*q = 0,
    A_ell = [[0,I],[-K_ell,-G]],
    G_uw=+2*Omega, G_wu=-2*Omega.

The K blocks are

    Lu = -d²/drho² + ell(ell+1)/rho² + S² + .3 F² - Omega²,
    Lw = -d²/drho² + ell(ell+1)/rho² + S² + .1 F² - Omega²,
    Leta = -d²/drho² + ell(ell+1)/rho² + 3S² - 1 + F²,
    M = 2FS,
    K = [[Lu,0,M],[0,Lw,0],[M,0,Leta]].

The D block is LD=-d²/drho²+ell(ell+1)/rho²+1+.01S². Its background is zero, so its quartic does not enter this linear problem. Eigenvalues s generate e^(s tau). Stable positive frequencies below use nu=Im(s)>0; this is the conjugate convention to e^(-i nu tau). Static Hessian eigenvalues are not squared rotating frequencies because G is nonzero.

All transformed perturbations obey Dirichlet conditions at rho=0,L. These are finite reflecting spherical boxes, not outgoing or infinite-volume scattering conditions. The origin condition and centrifugal discretization approximate the regular radial self-adjoint problem.

To avoid interpreting a sampled continuum profile's discrete residual as an instability, each grid first solves the stationary equations of the same centered-difference operator, at unchanged Omega, for p=rho F and z=rho(S-1). Newton iteration uses the original smooth profile as its seed. These discrete states approach the reference as spacing decreases; they are numerical approximations, not physical preparation changes. Examples at L30 are:

| N | h=L/(N+1) | Discrete Qbar |
|---|---:|---:|
| 120 | .2479338843 | 398.596525626 |
| 180 | .1657458564 | 399.062504484 |
| 400 | .0748129676 | 399.361905993 |
| 800 | .0374531835 | 399.419277406 |

The continuum reference has Qbar=399.438462422. At exactly h=.05, Qbar differs by only1.2e-10 between L30 and L40. Stationary maximum residuals and all Newton histories are retained.

For ell=0 the linear charge is proportional to

    Q1 = <p, 2*Omega*u_transformed - wdot_transformed>.

With v=(q,qdot), its coefficient vector c satisfies c^T A=0 when Lw p=0. The code verifies this identity. Every nonzero exact generator eigenmode consequently has Q1=0. For ell>0 the angular integral vanishes, giving Q1=0 directly. The rank-one fixed-charge term is used in a separate energy diagnostic; it is not arbitrarily inserted into the dynamical K matrix.

Across returned nonzero nonsymmetry modes, the largest normalized charge overlap was7.27e-13, and fine sparse values are around1e-14. Near-zero phase roots are excluded from this nonzero-mode diagnostic because division by a numerically split zero is ill-conditioned.

## Vacuum dispersion validation

The numerical vacuum generator reproduces the two charged branches sqrt(k_discrete²+1)−Omega and sqrt(k_discrete²+1)+Omega, the neutral eta branch sqrt(k_discrete²+2), and the D branch sqrt(k_discrete²+1.01). At N80,L30,ell=0 and2, maximum charged/eta differences from these exact discrete dispersions are1.15e-13; the D difference is below2.5e-15. Real parts of the vacuum generator eigenvalues are below1.8e-14.

The infinite-volume tree thresholds in these units are:

| Branch | Threshold nu |
|---|---:|
| Charged lower | .199999998024 |
| Charged upper | 1.800000001976 |
| eta | 1.414213562373 |
| D | 1.004987562112 |

Finite-box vacuum frequencies use the allowed nonzero radial k and therefore sit above their thresholds. These are the tree thresholds of this operator, not one-loop pole-mass thresholds.

## Phase and translation modes

The global phase mode is q=(0,p,0), qdot=0. The phase generator residual and the charge-left-null residual reach6.87e-15 and3.64e-15 in the tightened N180 test. The phase generalized vector constructed from the differentiated discrete stationary equations obeys A*v_generalized=v_phase, with residual6.23e-14 there. Its charge is −dQbar/dOmega, approximately+4093 on the finer grids, and it is excluded from the fixed-charge tangent space. Global phase is a collective coordinate, not an internal oscillator to put in a logarithm.

An apparent phase growth pair s≈±3.4022e-7 at the original N180 stopping tolerance was retained. One additional Newton refinement changes the stationary maximum residual from1.47e-11 to3.73e-14 and the pair becomes s≈±i1.6464e-8. Its phase overlap is effectively one. This tolerance dependence and the exact phase/charge identities diagnose a numerical Jordan splitting, not a converged positive-growth mode. Other nonzero frequencies are unchanged by this tightening.

Translations lie in ell=1 and have q proportional to (rho Fprime,0,rho Sprime), with multiplicity three. Their finite-difference frequency approaches zero:

| N, L30 | Translation nu |
|---|---:|
| 120 | .00704572194 |
| 180 | .00470278810 |
| 400 | .00212059226 |
| 800 | .00106141784 |

The overlap with the smooth translation vector tends to one. At h=.05 it is .99997361, and nu=.0014170633355 at L30 versus .0014170633339 at L40. The finite translation frequency is therefore dominated by radial discretization rather than box size. A fit on the finer L30 controls gives nu²≈.00080357*h² with intercept−7.1e-10, consistent with the exact continuum zero within the observed discretization scale. The small finite frequency is not claimed to be a physical internal mode. The associated collective drift and generalized zero directions require their usual separate treatment.

## Internal-frequency and growth search

No nonsymmetry internal frequency was resolved in the frozen window 1e-5<nu<.199998, for ell=0,...,5. This is a finite-box numerical statement, not an exclusion of an arbitrarily shallow state or a complete continuum spectral theorem.

At exactly matched h=.05, the smallest positive nonsymmetry frequencies are:

| ell | L30,N599 | L40,N799 |
|---|---:|---:|
| 0 | .20661186569 | .20354848812 |
| 1 | .21171086850 | .20643457717 |
| 2 | .21666015239 | .21017881045 |
| 3 | .22674909993 | .21514328231 |
| 4 | .23652700375 | .22070862612 |
| 5 | .24749970926 | .22698922830 |

They move toward the charged continuum threshold as the box enlarges. They must not be presented as six volume-independent internal excitation gaps. Radial refinement at fixed L is separately recorded in the raw spectra.

Apart from the identified phase numerical split, no growing root was resolved. The largest computed absolute real part among returned nonsymmetry roots was1.93e-14. Generator residuals across all returned modes reach1.57e-10, so the real-part figure is a computed value, not a rigorous1.93e-14 error bound on continuum growth. Dense searches are complete only for their finite matrices; refined sparse searches cover the declared shifts and returned modes rather than every fine-grid eigenvalue.

## Fixed-charge energy and higher angular momentum

The independent energy diagnostic gives a positive real-amplitude fixed-charge ell0 minimum .366966510879 at h=.05,L40. The corresponding uncorrected grand-canonical amplitude operator has a negative eigenvalue−.188833915084. That negative local direction is not itself a dynamical growing mode at fixed charge: minimizing the phase velocity subject to Q1=0 produces exactly the positive amplitude rank-one term4*Omega²|p><p|/||p||². The ell0 phase block has the symmetry zero and positive remaining eigenvalues.

At h=.05,L40, ell2 amplitude and phase minima are .285521978057 and .380710532650. For all ell>=2, each radial block differs from ell2 by the positive centrifugal term[ell(ell+1)−6]/rho², and D is positive. This is an algebraic all-higher-ell energy ordering; it does not require an unbounded numerical angular scan. Its continuum use remains conditional on the limiting background/operator and the corresponding nonnegative constrained-energy hypotheses. The companion `DERIVATION.md` supplies the more precise phase/translation ground-state identities and sufficient hypotheses. In particular, numerical positivity is not silently substituted for a certified continuum proof.

The separately owned angular comparison computes B=sup rho²[-lambda_min(Klocal−Klocal_vac)]_+ on the saved profile. Its optimized numerical value is28.0100296087 near rho4.24595; the conservative numerical comparison value is31. For real |nu|<.2, vacuum-pencil positivity together with ell(ell+1)>B would exclude a subgap mode. Since6*7=42>31, this condition begins at ell6, and ell5 was explicitly added above. `angular_bound_results.json` expressly marks the half-line bound as uncertified. Accordingly, finite ell0–5 searches plus this comparison provide conditional high-ell subgap evidence, not a complete continuum certification. This real-frequency subgap argument and the constrained-energy growth criterion are distinct.

## Evidence and limits

Raw records are `mode_results.json`, `mode_eigenvalues.csv`, and `mode_execution.jsonl`; protocols and the initial-result archive accompany them. `mode_internal_vectors.npz` contains no nonzero internal vectors because no mode met the internal-window criterion. The empty archive is expected, not a missing mode calculation. All symmetry roots and their overlaps remain in the numerical results and CSV, and the executable reconstructs the operators and consistent backgrounds.

The outcome supports the same classical saddle as a starting point for a future constrained quantum calculation. It neither evaluates a determinant nor supplies its continuum normalization, ultraviolet/Born subtraction, counterterms, collective-coordinate measure, fixed-Q quantum energy or lifetime. No frequency sum was performed. Finite boxes, radial grids, bounded sparse probes, the near-threshold window and the uncertified half-line angular bound remain explicit limitations. Gravity and Standard Model fluctuations are not included in this tree scalar calculation. Original candidate, classical-stability and previous quantum files were not modified.
