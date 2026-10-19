# Preparation, evolution, records and measured output

All core values below come from the preregistered parameter point in
[protocol.json](protocol.json). Time is in \(\hbar/J\); no conversion to
seconds or fit to nature is available. Exact finite-system spectral
evolution is compared with Krylov evolution of an independently assembled
Hamiltonian. This is algorithmic cross-checking by the same worker.

## An exact regime before numerical extrapolation

At \(\epsilon=0\), define the domain-wall operator

\[
D=\frac12\sum_i(1-\tau_i^z\tau_{i+1}^z).
\]

A permitted memory flip has opposite neighbors. It moves a domain wall
without creating or destroying one. All other Hamiltonian terms commute
with \(D\), so \([H_{M1}(0),D]=0\). The independent matrix calculation
returns an exactly zero commutator at \(L=3,4,5\).

In the zero-wall subspace,

\[
U(t)(|\psi\rangle_\sigma\otimes|r\cdots r\rangle_\tau)
=e^{iKLt}\,e^{-itH_\sigma^{(r)}}|\psi\rangle_\sigma
\otimes|r\cdots r\rangle_\tau,
\]

\[
H_\sigma^{(r)}=J(1+\alpha r)\sum_i\sigma_i^z\sigma_{i+1}^z
+h_x\sum_i\sigma_i^x+h_z\sum_i\sigma_i^z.
\]

Thus a perfectly retained classical record can coexist with nontrivial,
record-dependent matter dynamics. This exact result does not force those
matter dynamics to have a small or zero chaos exponent.

It also exposes a formation cost. An initially all-\(+x\) memory has
zero-wall probability \(2^{1-L}\). When \(\epsilon=0\), unitary evolution
cannot increase this probability: the desired sector is invariant.
Perfect protection is not a mechanism for generating the protected state.

At nonzero \(\epsilon\), the first protected-record terms are

\[
\langle\tau_0^z(t)\rangle_r
=r[1-2\Delta^2\epsilon^2t^2+O(t^3)].
\]

For the specified matter readout preparation,

\[
\left.\frac{d\langle\sigma_1^y\rangle_r}{dt}\right|_0
=4J(1+\alpha r)+2h_z=\begin{cases}7,&r=+1,\\3,&r=-1.\end{cases}
\]

For a maximally mixed matter register and uniform memory, the adjacent
matter squared commutator begins as

\[
\frac12\langle[\sigma_0^x(t),\sigma_1^x]^\dagger
[\sigma_0^x(t),\sigma_1^x]\rangle
=8J^2(1+\alpha r)^2t^2+O(t^3),
\]

with coefficients 18 and 2. Double-commutator/first-commutator matrix
calculations reproduce all three formulas. The final readout sign and
probability were not inserted as a boundary condition.

## Fixed-law numerical regime

| Core run | Minimum memory success through \(t=10\) | Memory success at \(t=40\) | Best matter readout success |
|---|---:|---:|---:|
| M1, L=3 | 0.994999 | 0.985006 | 0.806447 |
| M1, L=4 | 0.987510 | 0.935065 | 0.778468 |
| M1, L=5 | 0.989132 | 0.940399 | 0.778460 |
| Saturation removed, L=4 | 0.621452 | 0.709936 | 0.741417 |
| Feedback removed, L=4 | 0.996241 | 0.994449 | 0.500000 |
| Exactly frozen, L=4 | 1.000000 | 1.000000 | 0.778577 |

For M1 at \(L=4\), the best registered readout time is 0.58. The two
matter-Y means are -0.491561535 and +0.622309749. Hence

\[
P_Y=\frac12+\frac{|-0.491561535-0.622309749|}{4}
=0.778467821.
\]

The simultaneous memory success is 0.999312115. Both \(L=4\) and \(L=5\)
use this same readout time. Sampling is explicit; the reported minima are
not rigorous all-continuous-time enclosures.

Removing feedback leaves a persistent memory but no matter readout signal.
Removing saturation leaves an interacting readout but substantially erases
the record. Both parts of the proposed mechanism have a tested operational
role. The exactly frozen control is retained as an important limiting
case, not offered as autonomous formation.

The core's largest norm deviation is \(1.29\times10^{-14}\); sampled total
energy deviations are below \(3.07\times10^{-14}J\). Independent
Kronecker assembly and Krylov evolution reproduce the main retention trace
within \(2.24\times10^{-13}\) and the readout means within
\(5.2\times10^{-15}\).

## Energy, measurement and preparation are not free

At \(L=4\), the ground energy is -10.8731790514. Both retention
preparations have mean energy 0.2 and excitation energy 11.0731790514
above that ground state. Their energy variances differ; they are not the
same state in every observable.

The readout preparations have energies 1.55 and -0.45. The nonselective
matter-Y instrument is
\(\mathcal D_Y(\rho)=(\rho+Y\rho Y)/2\). At \(t=0.58\), its changes
in system energy are -1.2010144303 and -0.8997375490, respectively.
Negative values mean energy leaves this modeled system during the
instrument; they do not imply a zero-cost measurement. An energy-conserving
detector completion must receive that energy. Detector reset is not modeled.

## Attempts that did not work

[challenge_protocol.json](challenge_protocol.json) was written before the
stronger attacks. They are not removed by the core's 13 passed conditions.

- If only one memory cell contains the bit, local guessing success at
  \(t=10\) falls to 0.801481 (L=4) and 0.751472 (L=5). This is not a
  register of independently protected cells.
- The bare matter-to-memory writing preparation reaches only 0.553352
  and 0.526376 maximum success over the tested window. It fails the
  registered write-and-hold requirement.
- Adding an explicitly different, unconstrained flip term
  \(-\kappa\sum_i\tau_i^x\) lowers minimum success through \(t=10\)
  to 0.875748 at \(\kappa=0.1\) and 0.764940 at \(\kappa=0.3\).
  The saturation rule is a substantive supplied law, not symmetry-free
  universal protection.

The separately named M1-W completion adds a modeled source/carrier and
keeps the core law and its failed tests unchanged. Its complete numerical
record, refinements and scattering qualifications are in
`writer_results_initial.json` and `writer_results_refined.json`, summarized
in [REPORT.md](REPORT.md).
The apparatus writes from prepared resources; this does not upgrade the
failed generic-formation claim.

## Equilibrium and scrambling diagnostics

At \(\beta=1\), the normalized connected core memory autocorrelation
after projecting out identity and total energy is 0.538351 at \(t=10\).
Its finite-size degenerate spectral plateau is 0.458812. Neither number is
the prepared-bit success or an infinite-volume localization proof.

The adjacent matter squared commutator reaches 1.41102 within \([0,4]\).
It demonstrates information spreading while the different prepared-record
experiment retains memory; the preparations and observables are kept
separate. It does not establish an MSS Lyapunov exponent.

The raw registered regulated OTOCs use Pauli operators with nonzero Gibbs
means. A separately recorded centered diagnostic removes those means;
its initial factorization defect is still 0.607085. There is no controlled
large-system exponential window. No \(\lambda_L/(2\pi T)\) is reported.
