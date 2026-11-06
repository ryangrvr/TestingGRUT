# Source matching and reciprocal energy

The real classical amplitude $A$ multiplies $\cos(kz)$ in the phase of
the mass of one Majorana fermion. All densities are per spatial volume.
At finite grid use the evolved covariance $C$ from the free in-vacuum
$C_0$, including the coherence between free positive- and negative-energy
states. With $d\Pi$ the Bloch and transverse phase-space measure,

$$
F_\Lambda[A](t)=\frac12\int d\Pi\,\mathrm{Tr}\big[C\,\partial_A H(A)\big],
\qquad
E_\Lambda(t)=\frac12\int d\Pi\,\mathrm{Tr}\big[CH(A)-C_0H(0)\big].
$$

The half is the Majorana factor for the doubled spinor. Unitary evolution
implies $\dot E_\Lambda=\dot A F_\Lambda$. The force *on* the scalar is
the negative of the conjugate force $F_\Lambda$. In the previous closed
pulse the reported work used $C-C_0$ under the trace for the force.
The removed term is the derivative of the explicit function
$\frac12\int d\Pi\,\mathrm{Tr}[C_0(H(A)-H(0))]$. It integrates to zero
net work when $A$ comes back to zero. Using it in a live scalar equation
would add a cutoff-dependent scalar potential. This step restores
the full force and subtracts a constant reference energy.

## Divergent and finite terms

In Euclidean notation write the mass as $S+i\gamma_5 P$ with
$S=M\cos\theta$ and $P=M\sin\theta$. The even determinant is obtained
from the squared operator $-\partial^2+X$ with
$X=S^2+P^2-\gamma^\mu\partial_\mu(S+i\gamma_5P)$. In flat space
the relevant spin traces are $\mathrm{tr}\,X=4M^2$ and

$$
\mathrm{tr}\,X^2=4M^4+4\big((\partial S)^2+(\partial P)^2\big)
=4M^4+4M^2(\partial\theta)^2 .
$$

The local one-loop divergences of this continuum fermion problem
therefore give a constant and a kinetic counterterm, up to total
derivatives. There are no gauge fields or curvature here. Higher
derivative operators with independently specified scalar EFT coefficients
are separate physical inputs. The kinetic renormalization is consistent
with the axial backreaction literature in its own setting.
[Adshead, Liu and Lozanov, Eq. 4.34](https://arxiv.org/abs/2203.09370) ·
[Adiabatic Yukawa renormalization](https://doi.org/10.1103/PhysRevD.95.105003)

Set $K_R(0)=K'_R(0)=0$ with respect to the invariant
$s=(\omega+i0)^2-k^2$. This specifies the split between the fermion
response and the renormalized scalar quadratic action. For positive
frequency, the absorptive part is

$$
\operatorname{Im}K_R(s)=-\frac{M^2 s}{32\pi}
\sqrt{1-\frac{4M^2}{s}}\,\Theta(s-4M^2).
$$

It follows from the two-body spectral density, including the spatial
mean square $1/2$ of the cosine and the Majorana half. The twice
subtracted dispersion relation gives

$$
K_R(s)=\frac{M^2s}{32\pi^2}\int_0^1 dx\,
       \log\left(1-\frac{s}{M^2}x(1-x)\right).
$$

Below threshold its imaginary part vanishes. The inverse transform
acts on the history of $A$ as a retarded linear force. The benchmark
checks its integrated work against the independent two-body integral.
[Previous spatial threshold check](../AUDIT.md#what-the-checks-eliminated)

## Matching on the finite one-mode grid

The full finite-cutoff linear response contains a contact term from
the second derivative of the mass vertex. For pair frequencies
$\nu=E_p+E_q$ and nonnegative matrix-element weights $g_\nu$,

$$
K_\Lambda(\omega)=C_\Lambda+
       \sum_\nu \frac{g_\nu\,\nu}{(\omega+i0)^2-\nu^2},\qquad
C_\Lambda=\frac12\int d\Pi\,\mathrm{Tr}[C_0\,\partial_A^2 H(0)] .
$$

For constant spatial phase the static spectral term cancels the contact
term exactly. At fixed $k>0$ the finite grid need not preserve boosts.
Match its static value and temporal slope to the same invariant kernel:

$$
a=K_R(-k^2)-K_\Lambda(0),\qquad
b=K'_R(-k^2)-\left.\frac{\partial K_\Lambda}{\partial\omega^2}\right|_0,\qquad
F_{\rm match}=F_\Lambda+aA-b\ddot A .
$$

Here $a$ has dimensions $M^4$ and $b$ dimensions $M^2$. The matching
is performed independently for each spatial mode tested. It fixes a
time-local quadratic correction for the reduced model. There is also
a diagnostic that subtracts the finite-grid linear tangent of the
actual integrator and replaces it with the continuum retarded result.
The two prescriptions should agree when both the continuum limit
and the numerical integration are controlled.

For the midpoint split the sampled static susceptibility differs
from the continuous-time finite-grid value by

$$
K_{\Lambda,\Delta t}(0)-K_\Lambda(0)=
\sum_\nu g_\nu\left[\frac{1}{\nu}
-\frac{\Delta t}{2}\cot\left(\frac{\nu\Delta t}{2}\right)\right]
=\frac{\Delta t^2}{12}\sum_\nu g_\nu\nu+O(\Delta t^4).
$$

The leading high-momentum integral in three dimensions scales as
the fourth power of a common cutoff. This explains why increasing
the cutoff at fixed timestep can make the force worse while the
outgoing energy remains stable. The second protocol switches to a
symmetric fourth-order composition and measures at complete step
endpoints. The new tangent uses all signed substeps of that same
composition.
[Yoshida composition](https://doi.org/10.1016/0375-9601(90)90092-3)

## Feedback energy identity

For the registered scalar inertia $Z$, frequency parameter $\Omega$,
and external force $J$, the matched equation is

$$
(Z-b)\ddot A+(Z\Omega^2+a)A=J-F_\Lambda.
$$

The energy ledger, up to the common vacuum constant, is

$$
E_{\rm tot}=E_\Lambda+
   \frac12(Z-b)\dot A^2+\frac12(Z\Omega^2+a)A^2,\qquad
\dot E_{\rm tot}=J\,\dot A .
$$

Thus the same matching enters acceleration and scalar energy. Both
counterterms are fixed from the vacuum linear reference before the
live pulse. A finite-grid calculation must still demonstrate numerical
convergence and energy balance. The tested scalar mode has a supplied
quadratic action and preparation. Other spatial harmonics and the
stress components needed for a metric evolution remain separate.
