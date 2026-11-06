# Instant preheating feasibility at a separate published point

**The vacuum source is explicit. The unprotected tree history fails
this analytic screen.** The scalar loop term becomes comparable to
the nominal kinetic energy at \(N\simeq2.97\). This happens well
before the formal radiation equality at \(N\simeq11.27\).
The result concerns the declared minimal realization and finite
renormalization choice. The preflight used a separate published
input and leaves the vector claim ledger at conditional transfer.
[Input protocol](PROTOCOL.md) · [Claim ledger](../CLAIM_LEDGER.md)

## Fixed mechanism and input

Use the nonoscillatory Peebles–Vilenkin example studied by Haro and
Aresté Saló, with their illustrative \(g=10^{-5}\),
\(\Gamma/\bar M=10^{-10}\) and \(H_0/\bar M=1.44\times10^{-6}\).
The source takes the post-crossing tree tail to be shallow in energy
and subdominant to kinetic energy. Here \(N=0\) and \(a=1\) at the
field crossing. The mediator \(\chi\) is a real, conformally coupled
scalar, with mass \(g|\phi|\), initially in its adiabatic vacuum.
Add the Yukawa vertex \(-y\chi\bar\psi\psi\), with one massless
Dirac daughter in vacuum. The published \(\Gamma\) sets the width
at \(m_\chi=g\bar M\). The width then grows with the actual mass.
This gives \(y=\sqrt{8\pi\Gamma/(g\bar M)}=0.01585\).
[Haro and Aresté Saló, Sec. 4.3](https://arxiv.org/abs/2004.11843) ·
[Felder, Kofman and Linde, Eqs. 3, 4 and 10](https://arxiv.org/abs/hep-ph/9812289)

Nonadiabatic passage through \(m_\chi=0\) creates excitations from
vacuum. The flat-crossing estimate has
\(k_*^2=g|\dot\phi_0|=g\sqrt6\bar M H_0\),
\(n_k=e^{-\pi k^2/k_*^2}\), and
\(n_{\chi,0}=k_*^3/(8\pi^3)\). This estimate gives
\(k_*/H_0=4.12\), or \(H_0/k_*=0.242\), and a mean
1.18 mediator quanta per Hubble sphere at the crossing. The
0.1 screening threshold for a controlled flat burst is missed.
A quantitative production test would use expanding-space modes
through a smooth crossing. The Yukawa channel supplies relativistic
daughter energy even with zero daughter occupancy. The screened
daughter is a freely streaming Dirac species.

## Test of the tree expansion

Under the kinetic reference solution \(H=H_0e^{-3N}\),
\(\phi=\sqrt6\bar M N\), and
\(\Gamma(N)=\Gamma_0\sqrt6 N\). The rest-frame decay optical depth is
\[
 \tau(N)=\frac{\sqrt6\,\Gamma_0}{9H_0}
 \bigl[e^{3N}(3N-1)+1\bigr].
\]
The optical depth accumulated before the typical momentum drops
below the mediator mass is \(5.05\times10^{-6}\).
At median decay the typical momentum-to-mass ratio is 0.00504.
These check the nonrelativistic decay estimate within the tree model.

The comoving radiation energy is
\[
 C_r=n_{\chi,0}\int_0^\infty
 m_\chi[N(\tau)]e^{N(\tau)}e^{-\tau}\,d\tau .
\]
After decay the kinetic and daughter densities dilute as \(a^{-6}\)
and \(a^{-4}\). Their formal equality is
\(N_{\rm eq}=\tfrac12\ln(3\bar M^2H_0^2/C_r)\).
Using the source's separate constant-mass, sudden-decay approximation
recovers its quoted scale \(2.18\times10^7\) GeV. The mass-dependent
decay quadrature gives a different tree-level estimate.

| Diagnostic | Calculated value |
|---|---:|
| Median decay in e-folds | 2.832 |
| \(C_r/(3\bar M^2 H_0^2)\) | \(1.64\times10^{-10}\) |
| Formal tree equality in e-folds | 11.267 |
| Thermal-equivalent energy scale at formal equality | \(7.23\times10^7\) GeV |
| E-fold where the scalar loop term matches the tree kinetic density | 2.975 |
| Formal tree radiation-to-kinetic ratio there | \(3.24\times10^{-8}\) |

The temperature uses \(g_*=106.75\) solely to express an energy
density. The specified freely streaming daughter has a nonthermal
spectrum. Visible thermalization would require its own interactions.

The real mediator gives the one-loop term
\[
 V_1(\phi)=\frac{g^4\phi^4}{64\pi^2}
 \left[\ln\!\frac{g^2\phi^2}{\mu^2}-\frac32\right].
\]
In the declared \(\overline{\rm MS}\) screen at fixed \(\mu=g\bar M\)
with the uncompensated runaway tree tail, \(V_1/\rho_\phi=1\) at
\(N=2.975\). There \(m_\chi/H\simeq3.8\times10^5\), so the mediator
is heavy on the curvature scale, and only 67.5% of the tree
mediators have decayed. Extrapolating the same tree trajectory
to formal equality gives \(V_1/\rho_\phi\simeq1.73\times10^{24}\).
The loop estimate tests the finite effective action chosen here.
A symmetry or explicitly matched compensating potential defines
a separate realization. [Real-scalar one-loop term, Eq. 3.2](https://arxiv.org/pdf/hep-ph/0111209)

## Perturbations, cost and decision

Using the source's slow-roll values and negative-field potential
gives \(m_\chi/H_p\simeq9.91\) at the pivot. The mediator is heavy
there. A common inflaton time shift moves the crossing and all
local decay clocks together. It supplies the long-wavelength
adiabatic solution for the stated vacuum and constant couplings.
An independently generated variation of the comoving daughter
energy would give \(\delta N=\tfrac14\delta\ln C_r\) on a late
uniform-radiation-density slice. The actual spectrum of that
variation requires the expanding-space production and decay state.

The present screen is two one-dimensional quadratures and scalar
roots. With a radiatively controlled potential specified, a first
production convergence test could use about \(10^3\) momentum modes
and \(10^4\) time steps, of order \(10^7\) mode updates. Retaining
two complex state variables at every such step would use roughly
320 MB in double precision. These are planning estimates.
Transport to a visible bath and the independent scalar response
require a specified daughter sector before a reliable cost estimate.

The flat-burst approximation and the unprotected tree trajectory
both miss their control checks. Stop this realization at analytic
preflight. The formal equality scale is an extrapolation of the
tree reference past its control limit. The review target is the
input choice and the radiative-control comparison.

```sh
python screen.py
```

[Computed values](results.json)
