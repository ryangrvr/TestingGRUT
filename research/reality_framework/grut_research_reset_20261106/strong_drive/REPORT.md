# Propagation changes production on the frozen strong pulse

**Decision.** The cell-local reduction overproduces both fermion number
and net energy transfer on the reserved order-one spatial drive.
Including propagation gives **25.95% less final number** and **20.26%
less net work**, using the local values as denominators. The differences
exceed the registered numerical allowances by factors of 22.0 and
10.0. Close the cell-local source route for this tested regime.

This is a controlled comparison of a quadratic Majorana field on a
prescribed flat-space background with one periodic spatial direction.
The comparison leaves H55's actual background, the renormalized
instantaneous source, the interacting transport transition, and the
relic open.

[Registered pulse and rule](PROTOCOL.md) ·
[Machine-readable decision](results/decision.json) ·
[Independent observable screen](OBSERVABLES.md)

## Matched evolution and normalization

In units of the fermion mass, the free block for momentum $(r,0,p_z)$
and the driven mass residual are

$$
H_0=\tau_z\otimes(r\sigma_x+p_z\sigma_z)+\tau_x\otimes I_2,\qquad
V=(\cos\theta-1)\tau_x\otimes I_2-\sin\theta\,\tau_y\otimes I_2.
$$

The transverse azimuth has been rotated away. The propagating arm
retains the spinor at each $p_z=\kappa+jk$ with
$-k/2\le\kappa<k/2$ and $-L\le j\le L$. The sector measure is
$d\Pi=r\,dr\,d\kappa/(2\pi)^2$. In the open regulator $\cos(kz)$
connects adjacent bands with matrix element $1/2$. Diagonalize it
with a sine transform, apply the exponential of the projected
cosine in the mass term, and transform back. The cyclic check
changes only this finite-band prescription. The local arm uses
the same free spinor and measure, evolving independent homogeneous
copies at uniform phases of the spatial wave. Both drives end in
the same free Hamiltonian.

For the occupied spinor frame $U$, the covariance is $C=UU^\dagger$.
Outgoing particle and hole traces enter with $1/2$ for one Majorana
field, giving the free energy density difference
$\Delta\rho=\frac12\int d\Pi\,\mathrm{Tr}[H_0(C-C_0)]$.
The second energy calculation uses the work density

$$
w=\frac12\int dt\,d\Pi\,\mathrm{Tr}[(C-C_0)\partial_t H].
$$

The subtracted term integrates to zero on the closed pulse because
the initial vacuum covariance $C_0$ is time independent and the
initial and final Hamiltonians coincide. Evaluating the explicit
time derivative at split midpoints checks net transfer independently
of the outgoing projection. This subtraction defines the reported
finite-grid transient power.

## Qualification before the reserved pulse

The tangent to the actual mass kick at zero drive connects neighboring
bands. Each matrix element uses the sum of the free incoming hole and
outgoing particle energies in the midpoint time transform. This
exercises pair kinematics and normalization together.

| Gaussian width $M\tau$ | $k/M$ | Propagating weak $w/(\theta_0^2 M^4)$ | Local weak $w/(\theta_0^2 M^4)$ |
| ---: | ---: | ---: | ---: |
| 12 | 1.0 | 0.4796655858 | 0.8102005805 |
| 12 | 1.6 | 0.0004684312 | 0.8102005805 |
| 24 | 1.0 | 0.9646477704 | 1.6176511903 |
| 24 | 1.6 | $3.46632\times10^{-9}$ | 1.6176511903 |

The first qualification missed the narrow local reference by 0.61%.
Refining only the weak-test quadrature from 64 to 80 Bloch nodes and
12 to 24 transverse nodes per panel passed the unchanged gate. The
worst final relative error against the archived references was
$3.45\times10^{-7}$. A finite-amplitude weak control matched the
tangent to $3.85\times10^{-6}$ relative or better. A dense, unsplit
DOP853 evolution matched the split sector's number and energy to
$1.50\times10^{-5}$ relative. The static vacuum resolved no signal
above roundoff. The homogeneous arms agreed to $4\times10^{-12}$
absolute on the control grid.

[First gate output](results/qualification_first.json) ·
[Passing gate output](results/qualification.json) ·
[Weak invariant references](../results/checks.json)

## Reserved pulse and refinement

The inherited pulse has amplitude 1.2, half-width $12/M$, carrier
frequency $2.4M$, and wavenumber $1.6M$.

| Grid change | Propagating $n/M^3$ | Propagating $w/M^4$ | Local $n/M^3$ | Local $w/M^4$ |
| --- | ---: | ---: | ---: | ---: |
| Base | 0.055677688 | 0.133075395 | 0.075636098 | 0.167961169 |
| Timestep | 0.055671776 | 0.133062868 | 0.075637481 | 0.167954937 |
| Momentum spacing and phase quadrature | 0.056178296 | 0.134728940 | 0.075873931 | 0.168996542 |
| Longitudinal cutoff | 0.055681033 | 0.133095221 | 0.075636099 | 0.167961172 |
| Transverse cutoff | 0.055741605 | 0.133386832 | 0.075713614 | 0.168337018 |
| Combined fine | 0.056237321 | 0.135043162 | 0.075948113 | 0.169356013 |
| Fine cyclic regulator | 0.056237321 | 0.135043163 | 0.075948113 | 0.169356013 |

The measured allowances sum the one-axis shifts, crossed-grid
remainders, and regulator shifts in the two arms. Work also
includes the fine budget residuals. The allowances are
$8.98\times10^{-4}M^3$ in number and $3.43\times10^{-3}M^4$
in work, or 1.18% and 2.03% of the larger result. The gaps are
$1.97\times10^{-2}M^3$ and $3.43\times10^{-2}M^4$.
Momentum quadrature is the largest measured numerical contribution.
The fine work–outgoing energy discrepancies are $1.45\times10^{-5}$
and $1.40\times10^{-5}$ relative in the propagating and local arms.
These are measured refinement envelopes for one pulse and two
high-momentum prescriptions.

Two diagnostics were computed after applying the decision rule.

- The linear-response work ratio for the *same compact envelope* is
  $w_{\rm prop}/w_{\rm local}=0.1122$. At amplitude 1.2 the measured
  ratio is 0.7974. The compact envelope supplies above-threshold
  frequencies and the nonlinear evolution adds channels and
  interference. The change in ratio measures the breakdown of
  linear extrapolation. Attribution to individual harmonics
  would require a channel-resolved calculation.
- The propagating positive and negative work integrals are
  $+2.9667M^4$ and $-2.8317M^4$. With the cyclic regulator the
  positive integral moves by $0.01268M^4$ and net work by
  $5.11\times10^{-10}M^4$. At matching sampled times the largest
  shift is $0.00550M^5$ in power and $0.00237M^4$ in cumulative work.
  Final energy closure cannot certify the work history. A live
  scalar application needs a controlled subtraction for the
  transient source and a convergence test on that source.

## Replay

From this directory, use Python with NumPy and SciPy. The saved
run used NumPy 2.3.5 and SciPy 1.17.0.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python qualify.py
for case in base time spacing longitudinal transverse fine regulator; do
    python propagation.py "$case"
done
python decide.py
```

Per-grid JSON preserves momentum bins, sampled power and cumulative
work, the outgoing pressure trace divided by three, signed work,
and number and energy ledgers. Qualification and the base run
replayed to byte-identical JSON in a fresh output directory.
