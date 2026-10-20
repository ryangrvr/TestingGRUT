# Candidate 1: audited exclusion and the independent-input boundary

20 October 2026. Same-author research audit; independent review still owed.
Baseline: `479f111117ea7162bd7888ff8f3792b2f6c45fdb`.

The unchanged 100 GeV, YQ=YB preparation remains excluded **within its stated
leading-order, flat-FRW, standard-BBN and observational-sensitivity assumptions**.
The energy/charge arithmetic survives an independent derivation from the action
and a separate implementation. The earlier description of the observed H0 as
entirely cosmology/CMB independent needs correction: its published inference
uses a Planck matter-density cosmography. A direct distance-envelope check
avoids that input and retains the benchmark rejection.

The action and available evidence do **not** independently fix the scalar
scale and primordial hidden charge yield. This audit therefore does not select
a new successful point, predict an abundance by fitting one, or establish
positive empirical identification. No action extension or production campaign
was undertaken. No unavailable Claude results were used.

## 1. Evidence baseline and claim–evidence table

Read [the original report](../grut_3d_empirical_20261020/REPORT.md) for its 29
nonlinear 3D runs, original-mode controls, spectral calculation and production
network. Their files and claims were preserved. The audit verified
65 manifest-indexed packet files, 3
preserved parents, 20 archives and all 190
raw members (37,339,295 bytes). Hash agreement is integrity,
not scientific replication.

| Claim | Evidence now available | Supported scope / remaining gap |
|---|---|---|
| Finite modeled objects endure disturbances | Old 26 finite-object/control cases; about 99.85% isolated and 97% pair charge retention | Declared finite 3D perturbation family through tau160; no cosmological lifetime. |
| Old transverse norm measures an internal breakup rate | Original replay has material translation/rotation residuals; direct core motion/shape controls | **Not established.** Chiefly surviving-core deflection over the tested interval. |
| Isolated branch is a local fixed-charge minimum | Old numerical Hessian/zero-mode refinement and analytic angular ordering | Leading isolated branch; no rigorous global/quantum lifetime certificate. |
| Charge has a leading energy floor | Current/potential/square-completion derivation, saturation and negative controls | Canonical classical proper matter energy; not arbitrary ADM or quantum stress. |
| 100 GeV, r=1 preparation meets observations | Audited floor plus frozen direct six-maser distance-envelope comparison | **Excluded within stated conditional scope.** No retuning. |
| A converged relic population and individual survival history exist | Old threshold/grid-dependent network; no new population computation | **Not established.** Not needed by this total-charge test. |
| Scalar scale and hidden yield are independently known | Model parameters and initial charge are free supplied inputs | **Missing.** BBN calibrates baryons, not equality to hidden charge. |
| A new point is observationally compatible | Only necessary conditional product constraints derived | No new point selected; no sufficient viable region or full cosmology fit. |
| Observations favor the new mechanism over alternatives | Generic massive charged sectors share this gravitational inequality | **No such evidence.** Adopted GR/SM success is not scalar-sector identification. |
| An independent researcher reproduced these results | Same-author alternate derivation/code and integrity checks | **No independent scientific replication.** Checklist supplied. |

## 2. Exclusion audit

For `X=(q1+i q2)/sqrt(2)`, the specified current gives
`nQ=q2 dot(q1)-q1 dot(q2)` and `nQ=omega f^2` for the rotating solution.
Physical and dimensionless quantities obey

\[
Q=\bar Q/g_X^2,\quad E=(v/g_X)\bar E,\quad
E/Q=m_X\bar E/\bar Q.
\]

Directly minimizing the potential ratio and completing the kinetic squares yields

\[
m_*=v\min(g_X,(\lambda_\phi\lambda_X)^{1/4}),\qquad
\rho_h\ge m_*|n_Q|.
\]

Here `mu=m*/mX=0.562341325190` and `m*=56.23413252 GeV`.
Positive spatial gradients and D energy do not weaken the bound. It holds on
an instantaneous FRW slice at the adopted leading order; it does not assume
stationary cores, a formation efficiency, a cold hidden EOS or surviving
individual composites. Gravitational binding is not included as positive
matter energy; translating to mean cosmic gravity requires the hypotheses
in [ASSUMPTIONS.md](ASSUMPTIONS.md).

Conserved mean net charge and visible comoving entropy carry the supplied
initial yield forward. At present `r=abs(YQ0)/YB0` gives

\[
\omega_{h,\min}=\frac{m_*r}{\bar m_b}\omega_b,
\qquad H_0\ge100\sqrt{\omega_b+\omega_{h,\min}}.
\]

The audit obtains `omega_h,min=1.3387080801` and
`H0>=116.66353672 km/s/Mpc`.
Three-width-low BBN density and the old 1% upward baryon-mass allowance give
`H0>=113.25181752`. These reproduce the old results
to roundoff. The SI/entropy route gives
`YB=8.6758138605e-11` and
`nB=2.5083693055e-07 /cm^3`. The entropy/temperature conventions cancel
from the main density ratio; they do not calibrate a hidden asymmetry.

Although omega_b is written Omega_b*h^2, its denominator rho_c/h^2 uses the
fixed unit H100=100 km/s/Mpc. BBN's physical baryon-density calibration therefore
does not consume the maser H0 being audited. Standard-BBN assumptions and the
photon-temperature history remain conditions on that calibration.

The local inequality is not an ADM-mass theorem or an all-orders quantum bound.
Flat FRW, nonnegative other mean energy/Lambda, standard BBN, preserved late
baryons, no dilution of the specified initial yield, and negligible strong
binding/collapse/backreaction are load-bearing. The old short simulation and
small initial hidden fraction do not derive this complete cosmic history.
In particular a lower energy bound cannot prove an upper hidden energy at BBN.
See [DERIVATION.md](DERIVATION.md) for the complete proof and its domain.

## 3. Correcting the observational inference and freezing one target

Pesce et al. 2020 Eq. 1 uses flat LambdaCDM with Omega_m=.315 from Planck to
infer the familiar combined H0. The original floor's independence from hidden
pressure does not make that observed estimator independent of cosmography.
[CORRECTION.md](CORRECTION.md) records this qualification. The benchmark's
arithmetic is retained, not quietly changed.

[FEASIBILITY.md](FEASIBILITY.md) compares eight targets, including CMB, BBN,
compact-object lensing, ordinary scattering, the ideal D reporter and B-L
tests. The selected target is a conservative upper angular-distance envelope.
No population or cosmological composite lifetime is needed. At each nearby z,
the same charge floor and conserved baryon density imply

\[
H(z)\ge H_{\rm floor}(1+z)^{3/2},\qquad
D_A(z)\le\frac{2c}{H_{\rm floor}(1+z)}
\left[1-(1+z)^{-1/2}\right].
\]

This uses neither a fitted hidden EOS nor the Planck Omega_m prior. It still
uses ordinary GR/photons, flat representative FRW light paths, standard-BBN
calibration and the astrophysical disk/velocity interpretation. The CMB-frame
velocity convention and photon-temperature history remain explicit inputs.

The [public source freeze](https://github.com/ryangrvr/TestingGRUT/commit/eb5a1084cb8e9eaf4bef8c6cfa1f9cbd8b93d153) specifies source, all six Table 1
transcriptions, calibration, estimator, two velocity allowances and decision
rules before new numerical evaluation. [PROTOCOL.json](PROTOCOL.json),
[FREEZE.json](FREEZE.json) and [PUBLIC_FREEZE.json](PUBLIC_FREEZE.json) retain
exact identities. The analyst already knew the table and old rejection:
**retrospective, not blinded or prospective**. The direct distances and old
H0 inference are two uses of the same masers, not independent evidence.

## 4. Frozen comparison and conditional region

The primary sensitivity corners use BBN mean minus three reported widths,
1% higher mb, distance median minus three lower percentile widths, upper
systemic-velocity widths, and radial peculiar velocity -750 km/s (three times
250 rms). Exact multiplicative Doppler conversion is used. The predeclared
stress case doubles peculiar rms to 500, allowing -1500 km/s. Such corners
are not simultaneous confidence bounds or a full likelihood.

All distances below are Mpc. A row without a corner mismatch has only passed
this necessary envelope check, not established a viable model.

| Galaxy | Distance median | Lower sensitivity corner | Model upper, 250 rms | Model upper, 500 rms | Primary mismatch |
|---|---:|---:|---:|---:|---|
| UGC 3789 | 51.50 | 39.50 | 35.20 | 41.58 | yes |
| NGC 6264 | 132.10 | 81.10 | 90.97 | 97.00 | no |
| NGC 6323 | 109.40 | 40.40 | 72.08 | 78.23 | no |
| NGC 5765b | 112.20 | 96.90 | 77.84 | 83.95 | yes |
| CGCG 074-064 | 87.60 | 66.00 | 67.03 | 73.21 | no |
| NGC 4258 | 7.58 | 7.25 | 12.55 | 19.07 | no |

UGC 3789 and NGC 5765b remain incompatible at primary corners. NGC 5765b still
has `96.90 > 83.95` Mpc in the doubled-velocity stress case. The original
100 GeV, r=1 preparation remains the tested point; no mass/yield or dataset
selection was changed to recover agreement.

For fixed dimensionless couplings these data imply the necessary conditions

\[
m_X\,|Y_{Q0}/Y_{B0}|\lesssim 63.9344\ {\rm GeV}
\quad\text{(declared primary corners)},
\]

or `74.6429 GeV` with the doubled-velocity
stress prescription. These are conditional summary-based product bounds,
not calibrated confidence limits, a sufficient allowed cosmology, a selected
replacement preparation or an independently predicted abundance. Actual
energy above the floor can tighten them.

## 5. Independent inputs, accessibility and uncertainty

[PARAMETER_PROVENANCE.md](PARAMETER_PROVENANCE.md) answers the deciding question:
the model does not supply a measured scalar scale or an asymmetry-generation
history. At fixed scalar couplings, rescaling v rescales masses/lengths while
leaving the dimensionless flat dynamics unchanged. Exact U(1) conserves a
given charge; it does not choose its yield or generate net charge from a
neutral closed state. BBN calibrates YB; choosing YQ=YB adds a hypothesis.
No leading ordinary portal ties the hidden vacuum scale to the Higgs scale.

The allowed empty hidden classical vacuum also reproduces the absence of
additional classical hidden stress. A null gravity test cannot identify these
fields, and a generic charged massive component can share the same bound.
There are adopted framework postulates with separate falsifiable consequences
(including B-L), but no independently fixed discriminating **new scalar-sector**
prediction is available from this packet. This qualification is the intended
scope of the machine-readable global-model decision, not a claim that no
postulate of the action could ever be tested.

As an accessibility check only, the old reference object at 100 GeV has mass
`6.48752e-21 kg`, radius
`1.00704e-17 m`, compactness
`9.56806e-31` and
gravitational-lensing wave parameter for a 1 eV photon
`9.76599e-41`.
It is not a justified geometric-optics microlensing target. This says nothing
about an uncomputed high-mass tail. An event rate needs a population. The
ideal D probe has no supplied ordinary source/detector chain.

There are 22 new algebra/unit/implementation checks and
5 recomputation/integrity checks. Analytic versus
128-interval Simpson distance quadrature differs by
`6.661e-15` relative;
float versus 60-digit Decimal differs by at most
`1.110e-16`. This does not control uncomputed
EFT, cosmic-history, binding, disk-model or posterior-tail uncertainties.
BBN's reported density width alone gives about
`0.9404` km/s/Mpc
linearized Hfloor uncertainty. See [UNCERTAINTY_BUDGET.md](UNCERTAINTY_BUDGET.md).
No Gaussian significance or new raw maser/CMB likelihood is reported.

## 6. Decision, reproduction and stopping boundary

The five evidence categories remain separate:

1. **Internal verification:** the leading bound, unit conversion and numerical
   envelope are checked, and the frozen baseline/raw identities are preserved.
2. **Constraint/rejection:** the specified preparation is excluded within the
   stated physical hypotheses and observational-summary sensitivity treatment.
3. **Compatibility:** a necessary product region is derived; no new complete
   compatible preparation or relic solution is established.
4. **Evidence favoring the mechanism:** none supplied by this comparison.
5. **Independent reproduction:** none; same-author checking is not that result.

No new action version is adopted. The outcome for positive empirical validation
is: **the new sector cannot yet make a discriminating independently fixed
prediction because identified scale/yield inputs and physical access/history
calculations are missing.** The next dependency is an independent calibration
or approved, specified source of those inputs, followed by a genuinely reserved
discriminating observable. More local production runs do not fix that problem.

An additional portal, inflaton/asymmetry interaction or changed symmetry must
be a separately versioned extension with parameter provenance and owner approval
before becoming the empirical target. This audit stops at that boundary.

Run from this directory in a checkout containing the public source history:

```bash
python audit.py audit
python audit.py predict
python audit.py compare
python audit.py verify
python report.py
python quality.py
```

Only the Python standard library is required. [ENVIRONMENT.json](ENVIRONMENT.json)
records the execution environment. Source/data identities, numeric diagnostics
and every comparison row are in the JSON/CSV files; `PACKET_MANIFEST.json`
indexes the final packet. Old raw data remain in their unchanged archives.
The 12,000 state draws have a fixed seed and recorded aggregate diagnostics.
No new nonlinear field history exists or is implied.

[SOURCE_PROVENANCE.md](SOURCE_PROVENANCE.md) identifies original sources and
roles; [DEVELOPMENT_NOTES.md](DEVELOPMENT_NOTES.md) retains exposure, failures
and limits. The independent review still owed is specified in
[REVIEW_CHECKLIST.md](REVIEW_CHECKLIST.md). Continue from
[START_HERE.md](START_HERE.md); do not treat integrity counts as discovery.
