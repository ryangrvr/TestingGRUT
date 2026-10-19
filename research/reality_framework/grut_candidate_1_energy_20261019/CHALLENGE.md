# Challenger audit: proposed fixed-charge one-loop energy

19 October 2026. Same original Candidate 1 action, benchmark and matching. This file audits the rotating determinant, Born subtraction/restoration, charge projection and numerical remainder as they are actually supplied. It does not assume that a formal determinant or a finite numerical remainder already constitutes a physical fixed-charge quantum energy. No lifetime is calculated here.

## Immediate gates

- A physical fixed-charge prescription must resolve the phase/charge Jordan pair and the unconstrained negative static amplitude direction. A positive determinant at nonzero Euclidean frequency is not, by itself, a convergent unfixed-charge Gaussian definition.
- Vacuum subtraction must preserve the complete original vacuum oscillator count. Removing a phase coordinate from the numerator does not license deleting an arbitrary vacuum oscillator.
- Zero-mode factors have integrable infrared logarithms. Dividing them out without restoring their normalization and integral changes the result.
- The computed remainder and the restored diagrams must use the same Hessian potential, normalization, profile and spatial cutoff. A compact-support proof is not a proof for the original untruncated exponential tail.
- The UV and high-angular bounds below are conditional bounds in terms of their supplied Vmax and profile integrals. A numerically sampled supremum is not an interval-certified bound on the exact continuum profile.
- A one-loop Routhian at the reference frequency becomes a fixed-charge energy correction only through a consistent Legendre/stationarity argument. A classical energy and a pole threshold evaluated at inconsistent orders cannot establish quantum binding.

## Independent laboratory Hessian and harmonic check

Let the real charged components follow the background direction

\[
n(t)=(\cos\Omega t,-\sin\Omega t),\qquad
X=(x+iy)/\sqrt2.
\]

The laboratory real-field Hessian perturbations, relative to the vacuum masses squared `(2,1,1,1.01)` in units `(gv)^2`, are

\[
A=3(S^2-1)+F^2,\quad B=S^2-1+2bF^2,\quad C=bF^2,
\quad D=2FS,\quad E=.01(S^2-1).
\]

The charged block is `B I_2 + C M(t)`, where

\[
M(t)=\begin{pmatrix}\cos2\Omega t&-\sin2\Omega t\\-\sin2\Omega t&-\cos2\Omega t\end{pmatrix},
\]

and the organization/charged off-diagonal block is `D n(t)`. These are the full original Hessian entries, rather than the single static U'' potential used in simpler Q-ball approximations.

Expanding `Gamma1=(1/2) Tr log(H0+V)` gives `(1/2) Tr G0 V` and `-(1/4) Tr G0 V G0 V`. With

\[
\mathcal A(m^2)=\frac{m^2[\log(m^2/\mu^2)-1]}{16\pi^2},\qquad
B_0(s;a,b)=-\int_0^1dx\log\frac{xa+(1-x)b-x(1-x)s}{\mu^2},
\]

the first graph is exactly

\[
K_1^{(1)}=\tfrac12\{\mathcal A(2)\!\int A+2\mathcal A(1)\!\int B+\mathcal A(1.01)\!\int E\}.
\]

For the second graph, `M_+` and `M_-` have `tr(M_+ M_-)=1`; summing both harmonics gives `2 Ctilde^2`. The vector coefficients have `n_+ dot n_-=1/2`; both ordered field paths and both harmonics give `2 Dtilde^2`. Thus the proposed second graph passes the factor and sign check:

\[
\begin{aligned}
K_1^{(2)}=-\frac1{64\pi^2}\int\frac{d^3k}{(2\pi)^3}\{&\widetilde A^2B_0(-k^2;2,2)
+2\widetilde B^2B_0(-k^2;1,1)\\
&+2\widetilde C^2B_0(4\Omega^2-k^2;1,1)
+2\widetilde D^2B_0(\Omega^2-k^2;2,1)\\
&+\widetilde E^2B_0(-k^2;1.01,1.01)\}.
\end{aligned}
\]

The signs of the temporal invariants are positive: the laboratory harmonics carry real frequencies 2Omega and Omega. Treating them as static or inserting negative squared frequencies would restore the wrong diagrams. The charged C channel has maximal timelike invariant `4 Omega^2=2.56<4`; the mixed D channel has `Omega^2=.64<(sqrt2+1)^2`. The static channels are spacelike. Consequently these graph kernels remain below their two-particle cuts; the denominator check must still be executed by the implementation.

All cross terms absent from this formula vanish by species tracing or time-harmonic conservation. The UV divergent second-order local contraction is proportional to `A^2+2B^2+2C^2+2D^2+E^2=tr(V^2)`, independent of Omega. This agrees with locality of the scalar-potential counterterms. The first two graphs alone are not the complete one-loop soliton energy.

The graph outputs have units `K1/(gv)` as dimensionless numbers. The physical scalar-loop correction is `gv*K1`, while the classical energy is `(v/g)*Hbar`. Omitting their relative factor `g^2` would misidentify the loop order.

## Charge, contour and zero-mode normalization

With the earlier conventions, the Euclidean operator for Fourier convention `exp(i zeta tau)` is

\[
P(\zeta)=L+\zeta^2 I-2\Omega\zeta J.
\]

The other Fourier convention reverses the skew term and exchanges the charged conjugate branches, without changing the determinant. In a constant unitary charged basis the vacuum resolvents have squared masses `1+(zeta +/- i Omega)^2`. The Euclidean matrix is not a real positive symmetric stiffness matrix. Convergence and branch choices cannot be inferred just by treating it as one.

The fixed-frequency saddle has an unconstrained negative static amplitude direction. At fixed charge, eliminating the parallel phase velocity adds the positive `4 Omega^2 |F><F|/I` term in the **reduced** energy, while the original full gyroscopic operator remains unchanged. This is precisely the sector that a correct charge-projected canonical construction must address.

The supplied finite canonical toy is a useful direct check. Its exact polynomial identity is `det(Pfull)=zeta^2 det(Pred)`, with the phase coordinate removed and its Routh correction included only in the reduced operator. Its reduced stiffness is positive, reconstructed nonzero modes obey the full generator and charge constraint, and the complete-vacuum-normalized logdet integral agrees with half the nonzero-frequency difference to about `5e-16`. The unfixed charge curvature is negative, so the same toy explicitly rejects a naive unconstrained real Gaussian trace. This is finite-dimensional algebra/normalization evidence, not a continuum existence or convergence proof.

For an exact saddle the phase sector contributes `2 log zeta` at small positive zeta; the translation sector contributes the same per spatial zero coordinate, with multiplicity three. These logarithms are integrable. An analytic infrared subtraction must be added back. For example,

\[
\int_0^\infty d\zeta\log\frac{\zeta^2}{\zeta^2+a^2}=-\pi a.
\]

The energy prefactor `1/(2 pi)` makes this `-a/2` per zero mode. Deleting this factor without restoration therefore shifts the energy by `a/2`. Dividing by bare zeta powers also changes the ultraviolet degree. Collective-coordinate measures may give subextensive `log T/T` normalization factors, but that statement requires the fixed-charge/rest-frame construction; an arbitrary prime on a determinant is not the argument.

At strict one-loop order, the stationary Legendre conversion would give `E1|Q=K1(Omega_class)`: the frequency-shift term cancels because `K0'(Omega_class)+Q=0`. It does not require computing the O(hbar) profile shift to insert it back into the stationary classical functional. This cancellation is conditional on an actual correctly normalized K1 and a consistent definition of Q, not a replacement for either.

The reference continuous charge also needs explicit interpretation. With unit elementary U(1) charge, `Qref=Qbar/g^2≈39943.84624` is not an exact charge eigenvalue. A continuous semiclassical E(Q) interpolant is legitimate if labeled. For an integer sector N, the expansion about this reference additionally contains `omega_class*(N-Qref)` at the one-loop energy scale. Calling the floating reference an exact integer-charge eigenstate would be inaccurate.

## UV bound, including complex charged resolvents

Write `t=zeta^2+1-Omega^2=zeta^2+.36`, `Vmax=sup ||V||op`, and `I2=int d^3r tr(V^2)`. In the circular charged basis the free resolvents are `(h_l+t +/- 2iOmega zeta)^-1`. The radial free heat kernel is positive, and

\[
\left|\int_0^\infty ds\,e^{-(t\pm2i\Omega\zeta)s}e^{-s h_\ell}(r,r')\right|
\le G_{t,\ell}(r,r')\le G_{0,\ell}(r,r').
\]

The organization and reporter channels have still larger positive real masses and satisfy the same domination. This proves the needed complex-kernel estimate without incorrectly assuming real modified Bessel functions. The constant unitary change of basis preserves operator and Hilbert-Schmidt norms.

In three dimensions,

\[
\|G_0V\|\le V_{\max}/t=:r,\qquad
\|G_0V\|_{HS}^2\le\frac{I_2}{8\pi\sqrt t}.
\]

For `r<1`, Schatten inequalities bound every trace of order n>=3, giving

\[
|\operatorname{Tr}\log(1+G_0V)-\operatorname{Tr}G_0V
+\tfrac12\operatorname{Tr}(G_0V)^2|
\le\frac{r}{3(1-r)}\frac{I_2}{8\pi\sqrt t}.
\]

The proposed full-angular frequency-tail bound is therefore correct:

\[
|K_{\rm tail}|\le\frac{V_{\max}I_2}{48\pi^2}
\int_Z^\infty\frac{d\zeta}{(\zeta^2+.36-V_{\max})\sqrt{\zeta^2+.36}},
\quad Z^2+.36>V_{\max}.
\]

This controls the first-two-Born-subtracted remainder, not the separately restored graph tails or any error in Vmax/I2.

## High-angular bound and its support assumption

For ell>1 the massless radial Green kernel is

\[
g_\ell(r,r')=\frac{r_<^{\ell+1}r_>^{-\ell}}{2\ell+1}.
\]

If V is genuinely supported on r<=R, its Schur bound and squared-kernel integral give

\[
\|G_{0,\ell}V\|\le\frac{V_{\max}R^2}{(\ell-1)(\ell+2)},
\]

\[
\|G_{0,\ell}V\|_{HS}^2\le
\frac{2\int_0^R dr\,r^3\operatorname{tr}V^2}
{(2\ell+1)(2\ell+3)(2\ell-1)}.
\]

The constants in both proposed estimates pass the direct integration check. Their logarithm bound requires the displayed norm ratio below one. After multiplying by angular degeneracy, the order>=3 angular remainder is summable. This massless bound alone is independent of frequency; to obtain a finite energy bound, apply it on a finite frequency interval and use the separate UV bound above beyond that interval.

These inequalities do not bound replacement of the actual exponential tail by V=0. Spatial-cutoff variation supports numerical convergence but does not automatically supply a rigorous absolute interval. A rigorous full-profile result additionally requires profile enclosures/tail estimates and consistent zero-mode control. Approximate symmetry zero shifts can be especially important in the infrared.

## Primary-method check

Graham's [hep-th/9912286](https://arxiv.org/pdf/hep-th/9912286) explains the equivalence between Born subtractions and the corresponding Feynman diagrams, with counterterms fixed before comparing backgrounds. That method supports the subtraction/restoration principle, not importing its different model's finite on-shell counterterms into this MSbar benchmark.

Graham's [hep-th/0105009](https://arxiv.org/pdf/hep-th/0105009) uses a reduced single-potential scattering equation and on-shell conditions. Those formulas are not the coupled phi/X phase/amplitude operator here. No numerical correction or stability conclusion from that model is transferred to Candidate 1.

## Actual-output status and rejection rules

The supplied graph protocol and canonical toy source/results have been read. The harmonic signs/multiplicities, toy charge reduction and proposed UV/angular-bound constants pass the independent checks above. The continuum determinant derivation, actual spectral remainder, restored graph outputs and integrating report are still to be audited as supplied; no finite-Q quantum-energy verdict is presumed.

A converged growing mode, an unjustified absolute value/complex-part clipping, a mismatched Born subtraction, an un-restored zero-mode factor, an unbounded important tail, or an unproven charge-contour conversion would reject the corresponding physical-energy claim. Partial pieces should then be reported as partial pieces, without fitting finite counterterms, moving the matching or claiming a lifetime.

## Integrating-root follow-up on completed outputs

The preceding audit is preserved as the challenger's actual scope at the time.
The following is the integrating root's subsequent review, not a claim that
the challenger or an outside reviewer audited the completed calculation.

The completed continuum output contains 2362 successful point records and
complete six/eight-point quadratures through ell40. Every stored Jost phase
is below `6e-14`; the real-part extraction is not concealing a sampled
imaginary instability. The stricter 1312-point replay also passes its phase
gate. A separate finite-difference determinant/Born calculation converges at
two physical points, while a fresh canonical toy and vacuum/kernel controls
pass. These checks support the implemented normalization, not an exact
continuous-frequency stability theorem.

The original radial tolerance leaves an appreciable aggregate remainder bias:
the stricter domain/origin/tolerance replay changes the energy coefficient by
`+.00851301935`, despite a much smaller original quadrature change. Both
outputs remain archived. The `.02` empirical numerical allowance was declared
before that replay; it is not represented as a rigorous global numerical bound.

The root derived an additional weighted Birman-Schwinger angular estimate,
recorded separately with its proof, protocol and results. Its sharper angular
bound is `.05042782`, in addition to the `.03827594` frequency-tail bound.
The original, weaker compact-support estimates remain in the audit results.
Profile moment and supremum inputs enclose the padded numerical spline; they
do not interval-certify the exact nonlinear soliton. Spatial-cutoff/profile
errors retain numerical, not theorem-level, control.

The fixed-integer-charge energy expansion is checked by a new direct BVP;
the like-order free-X threshold uses the previously calculated strict pole
mass. The resulting negative free-carrier comparison is supported within the
reported numerical/saddle/loop qualifications. It does not settle residual
soliton emission, quantum fission, a lifetime, full continuum stability,
gravitational vacuum matching, higher-loop control or empirical viability.
Those gates remain open; they were not removed because this channel binds.
