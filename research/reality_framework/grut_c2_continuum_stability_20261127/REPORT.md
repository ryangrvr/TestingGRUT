# C2 fixed-charge continuum stability

**Verdict A.** All three original Candidate 1 stationary configurations underlying C1 have a certified classically linearly stable fixed-charge scalar fluctuation spectrum on the continuum half-line, modulo phase and translation symmetry modes. This uses the leading flat scalar dynamics of C1. The proof takes the exact stationary profiles from a contraction of a validated bordered map and uses interval spectral bounds on the continuum forms. [Certificate method](CERTIFICATE_METHOD.md) · [Amendment](AMENDMENT_01.md) · [Raw results](RESULTS.json)

## Exact saddles and radial modes

The table reports outward frequency and charge-slope intervals. The profile radius bounds the norm with components $\|(u,w)\|_2$, $\|(u,w)''\|_2/16$ and frequency. The number $g$ bounds the radial spectrum of $L_+$ away from zero, with one eigenvalue below $-g$ and the rest at or above $g$. The displayed gap bounds are rounded down.

| Charge | Exact $\Omega$ enclosure | Profile radius upper | $g$ lower | $d\bar Q/d\Omega$ enclosure |
| ---: | --- | ---: | ---: | --- |
| 39,944 | [0.79999952449, 0.79999972818] | $1.019\times10^{-7}$ | 0.03121578 | [−326.52232, −324.92872] |
| 39,943 | [0.80000196660, 0.80000217233] | $1.029\times10^{-7}$ | 0.03121543 | [−326.51768, −324.90815] |
| 19,972 | [0.88507412869, 0.88507426173] | $6.652\times10^{-8}$ | 0.03122765 | [−99.23759, −98.77926] |

Here $\bar Q=g_X^2Q/(4\pi)=\Omega I$ and $I=\int u^2\,d\rho$. The certified slopes invert the bordered system and, together with the radial index, imply zero negative index and zero kernel of $L_Q=L_++4\Omega^2|e\rangle\langle e|$. The positivity conclusion uses the exact selected continuum profiles. The largest stationarity defect is $1.338\times10^{-9}$ in the transformed $L^2$ norm, the largest bordered inverse bound is 38.460, and the largest Newton contraction bound is 0.001425. The exterior analytic tail has decay rate at least $1/3$ in $\rho$, with $\|(u,w)\|$ at the joining radius at most $3.086\times10^{-7}$ across the three enclosures. The lowest exterior spectral floor exceeds 0.21664. [Detailed enclosures](RESULTS.json)

## Angular and symmetry sectors

| Sector | $Q=39,944$ | $Q=39,943$ | $Q=19,972$ |
| --- | --- | --- | --- |
| Amplitude $\ell=0$ at fixed charge | Strictly positive | Strictly positive | Strictly positive |
| Phase $\ell=0$ | Phase mode and gap ≥0.01559078 | Phase mode and gap ≥0.01559043 | Phase mode and gap ≥0.01560265 |
| Amplitude $\ell=1$ | Translation modes and gap ≥0.01559078 | Translation modes and gap ≥0.01559043 | Translation modes and gap ≥0.01560265 |
| Phase $\ell\ge1$ and amplitude $\ell\ge2$ | Positive by ordering and exterior floor | Positive by ordering and exterior floor | Positive by ordering and exterior floor |
| Reporter, all $\ell$ | ≥1 | ≥1 | ≥1 |
| Reduced gyroscopic evolution | Linear energy stability | Linear energy stability | Linear energy stability |

The exact phase and translation eigenfunctions identify the kernels. The positive-shift index envelopes bound the number of low eigenvalues by one in each radial symmetry channel. Centrifugal ordering gives strict positivity for each higher sector and uniform control from its first member. The skew velocity term preserves the positive reduced quadratic energy when phase and rest momentum are fixed and the symmetry coordinates are quotiented.

The unchecked first interval LDL recurrence reached a zero-containing interval at block 123. Its frozen result is retained. The bounded corrected pass used the same mesh, exterior conditions, shifts and precision, with an a posteriori interval residual and inverse test. The largest Neumann product upper bound was $1.051\times10^{-37}$. Its exact replay and independent proof scrutiny are documented. [Negative result](INITIAL_ATTEMPT.json) · [Reproduction](REPRODUCTION.md) · [Internal critique](INTERNAL_CRITICAL_REVIEW.md)

C2 supplies the local scalar saddle premise for C1 at these three charges. C1's quantum energy allowances, nonlinear orbital stability, higher-loop physics, gravitational dynamics, lifetimes, and formation remain separate. This certificate warrants the bounded empirical-feasibility assessment of the prepared parent. External replication is pending. [External RAI review packet](RAI_REVIEW_PACKET.md)
