# Playing out the four next steps for GRUT

This follows the four next steps from "Can GRUT get gravity and particles?" through as calculations. Before running anything, I wrote the predictions and pass windows into `qed/PREREGISTRATION.json` (sha256 `d6808a17…b5a3`, checked after the runs). Misses are reported as misses, and no window was widened.

## Short answer

- **Interacting charged matter fits GRUT's detector-history framework.** In the lattice Schwinger model (1+1D QED), GRUT's two-time histories behave as its effective-description theorem predicts. All three preregistered history checks pass.
- **The spectrum check, which compares the computed spectrum with known continuum values, missed.** The zero-mass extrapolation lands 3.3% off the exact value against a 3% window. The massive holdout lands 9.9% off. Exact diagonalization at up to 24 sites is too small to reach the continuum to 3%. This is a limit of the method, not evidence against GRUT, but it still counts as a miss.
- **The light-cone requirement is much tighter than the photon-only ceiling.** If GRUT uses Wilson fermions, its electrons outrun its photons in every direction tested. The Crab Nebula's PeV electrons then push the spacing ceiling from about \(10^{-27}\) m down to about \(3\times10^{-32}\) m.
- **The chiral and spin-2 steps are stated as obligations, not solved.** No derivation was attempted. Below I write them in the form the GRUT program could adopt.

## Step 1: lattice QED in GRUT's framework

**Model.** The model uses Kogut–Susskind staggered fermions with U(1) links and open boundaries. Gauss's law is solved exactly, and θ=0. The Hamiltonian is \(W = x\sum(\sigma^+\sigma^-+{\rm h.c.}) + \tfrac{\mu}{2}\sum(1+(-1)^n\sigma^z) + \sum L_n^2\), with \(x=1/(ag)^2\) and \(\mu=2\sqrt{x}\,m/g\). The reference values come from the [exact massless result and Byrnes et al. (DMRG)](https://ar5iv.arxiv.org/html/hep-lat/0202014): vector binding energy 0.5641896 at \(m=0\) and 0.51918 at \(m/g=0.25\).

### Detector histories (P-C, P-D, P-E)

The instrument is GRUT's two-outcome Weyl instrument, \(K_\pm=(1\pm e^{i\varphi}W)/2\) with \(W=e^{icF}\). Here c=2, and F is a Gaussian-smeared electric field (σg=0.8 in a box with Lg≈5). Histories are computed at two times, out to T=6/g. The statistic is the maximum difference between the exact history probabilities and a quasi-free prediction built only from two-point functions.

| Lattice | \(m=0\) deviation | \(m/g=0.5\) deviation |
|---|---:|---:|
| N=10, x=4 | 0.00763 | 0.00788 |
| N=14, x=9 | 0.00551 | 0.00918 |
| N=20, x=16 | 0.00412 | 0.00895 |

- **P-C passes.** At \(m=0\), the deviation shrinks as the lattice gets finer. This is the expected behavior: the massless continuum field \(E=(g/\sqrt\pi)\phi\) is free, so the two-point summary should be enough. That is GRUT's effective-description criterion working on an interacting gauge theory.
- **P-D passes.** At \(m/g=0.5\), the deviation is larger on every lattice, though only barely at x=4. It does not shrink as the lattice gets finer, so the history framework detects the interaction rather than averaging it away.
- **P-E passes.** The first instrument disturbs the second: by about 0.04 at \(m=0\) and 0.012–0.015 at \(m/g=0.5\). Probabilities stay normalized to machine precision.
- **Caveat.** The two masses have different field fluctuation sizes (Var(cF) ≈ 0.24 vs 0.12). The comparison is therefore not perfectly matched, though the smaller fluctuations at \(m>0\) should make deviations harder to see, not easier.

### Spectrum (P-A, P-B)

**Recipe.** I fixed the recipe in `qed/extrapolate.py` using only the \(m=0\) data, before the holdout run finished.

1. At each x, fit \(\omega^2=M^2+b/N^2+c/N^4\) on N=14–24.
2. Extrapolate in \(1/\sqrt{x}\) (quadratic) over x=2–32.
3. Take the error band from the min–max spread over 24 fit variants.

| Prediction | Target | Result (variant spread) | Error | Window | Verdict |
|---|---:|---:|---:|---:|---|
| P-A, \(m=0\) | 0.5642 | 0.583 (0.561–0.641) | +3.3% | ±3% | Narrow miss |
| P-B holdout, \(m/g=0.25\) | 0.5192 | 0.571 (0.523–0.609) | +9.9% | ±3% | Miss |

The infinite-length values decrease steadily with x but are still far from the continuum at x=32. Examples: \(m=0\) runs 0.785 → 0.613; \(m/g=0.25\) runs 0.737 → 0.582. At large x, the bound state spans more lattice sites than a 24-site chain can hold. The qualitative trend is right, since binding drops when mass is added, but the size is not. The published values were reached with DMRG/MPS on chains hundreds of sites long. Matching them would take that tooling, which TestingGRUT does not yet have.

### Update: sealed DMRG rerun passes both

A DMRG/MPS solver was validated against ED (rel ≤ 2e-13), and its recipe was sealed before any production run (`mps/PREREGISTRATION_AMENDMENT_01.json`, commit `dd1e107`). It then reached chains of up to N=192 (physical length 24/g). P-A gives 0.5641906 against 0.5641896, and P-B (holdout) gives 0.519232 against 0.51918. Both are inside ±3%, with every fit variant within 0.25%. The ED misses above stay on record as the first test. See `RESULTS.md`.

## Step 2: the chiral obligation

The earlier probe showed that every local, Hermitian, translation-invariant lattice fermion in GRUT's style produces equal numbers of left- and right-handed species ([Nielsen–Ninomiya](https://en.wikipedia.org/wiki/Nielsen%E2%80%93Ninomiya_theorem)). Step 1 avoids the issue on purpose because QED is vector-like. The Standard Model is not, so GRUT has to name an escape route, in the same way the photon review named the plaquette-only choice:

- **Ginsparg–Wilson / overlap fermions:** keep a modified chiral symmetry exactly.
- **Domain-wall fermions:** chirality lives on a boundary of an extra dimension.
- **Mirror-fermion decoupling or symmetric mass generation:** give the doublers mass through interactions.
- **Non-lattice access structure:** drop one of the theorem's assumptions, which means saying which one.

The obligation reads: state which of these GRUT adopts, and show that its access or closure rules select it rather than import it. Nothing in the repository does this yet ([Book VII is an "absence map"](https://github.com/ryangrvr/TestingGRUT/tree/8e90c36ecf8f0986fbccb1097b327f37526959c4)).

## Step 3: universal spin-2 coupling as the target

The earlier light-bending test, calibrated on Newton's law, reduces the gravity problem to one statement. GRUT has to derive a massless spin-2 field that couples universally to the conserved stress tensor. A massless graviton is within 0.9σ of [Cassini](https://pubmed.ncbi.nlm.nih.gov/14508481/). A massive one gives 3/4 of GR's bending ([vDVZ](https://arxiv.org/html/hep-th/0011138)), and a scalar gives none; both are excluded by thousands of σ. `GR1_GRAVITY_VERDICT_01.md` already marks the minimal-stress coupling and the gapless sector as imports. [Weinberg–Witten](https://en.wikipedia.org/wiki/Weinberg%E2%80%93Witten_theorem) adds that the graviton cannot be a composite of the lattice's ordinary Lorentz-covariant conserved currents. It has to emerge without a conserved, covariant stress tensor at the lattice level. Alternatively, GRUT can say openly that it is a postulate.

## Step 4: spacing ceiling and shared light cone

`lightcone.py` expands both dispersion relations to \(O(k^2a^2)\), with \(v=1+C(ka)^2\):

- **Wilson fermion (r=1):** \(C_e=3/8-S/2\)
- **Plaquette photon:** \(C_\gamma=-S/8\)

Here \(S=\sum n_i^4\). Both formulas check numerically to about \(10^{-5}\).

- **The two particles do not share a light cone.** Electrons are faster than photons in 100% of sampled directions (equal only along the axes), and faster than light in 77.6% of them.
- **The ceiling tightens about 50,000-fold.** [LHAASO](https://arxiv.org/html/2204.02956v1) sees Crab electrons at 2.3 PeV or more with no vacuum Cherenkov cutoff. That requires \(a<3.0\times10^{-32}\) m along the body diagonal, or \(3.9\times10^{-32}\) m in the median superluminal direction. The photon-only ceiling from [GRB 221009A](https://ui.adsabs.harvard.edu/abs/2024PhRvL.133g1501C/abstract) was \(1.6\times10^{-27}\) m.
- **Caveats.** This is specific to Wilson fermions. A mass-dependent speed shift of order m·a would need tuning, and the Crab's direction relative to the lattice axes is unknown. A GRUT lattice could restore a common light cone by tuning, but that tuning is another import to declare.

## What this adds up to

| Step | Outcome |
|---|---|
| 1. Interacting charged matter | History checks pass (3/3). Spectrum checks miss (0/2) because of solver size. Needs DMRG/MPS. |
| 2. Chiral obligation | Stated, not solved. |
| 3. Universal spin-2 coupling | Stated as the target, not derived. |
| 4. Light cone | Shared light cone fails for Wilson fermions without tuning. Ceiling now about \(3\times10^{-32}\) m. |

GRUT's detector-history formalism holds up when interactions are switched on. What it lacks is a way to fix particle content and gravity's coupling from its own principles. A useful next build for TestingGRUT would be an MPS solver: rerun P-A and P-B against the same frozen windows, then add a GW/overlap fermion for step 2.

## Scope

- These are 1+1D and linearized calculations, written independently of the GRUT codebase. Nothing was committed to the GRUT repository.
- They have not been reviewed by the GRUT author or by outside experts.
- Files: `qed/PREREGISTRATION.json`, `schwinger.py`, `spectrum2.py`, `extrapolate.py`, `histories.py`, `lightcone.py`, and their JSON outputs.
