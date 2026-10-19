# Vacuum particle thresholds from momentum-dependent two-point functions

This calculation executes the next vacuum-threshold obligation at the original scalar matching. It includes all new scalar species, their tadpole-induced vacuum displacement and momentum-dependent scalar bubbles at strict one-loop order. It excludes Standard Model/gravitational loops and does not calculate a soliton quantum energy.

## Definitions and sign

Use canonical fields \((\phi,x,y,D)\), \(X=(x+iy)/\sqrt2\). The matching vacuum is \((v,0,0,0)\) at tree order. Write the positive timelike invariant as \(s=-p_E^2\), and define

\[
A(t)=\frac{t}{16\pi^2}\left[\ln(t/\mu^2)-1\right],
\quad B_0(s;a,b)=-\int_0^1 dz\,\ln\frac{za+(1-z)b-z(1-z)s}{\mu^2}.
\]

At matching \(a=M_{\phi,0}^2=0.02v^2\), \(b=M_{X,0}^2=0.01v^2\), \(M_{D,0}^2=0.0101v^2\), and the mixed cubic vertex is \(g_3=2cv=0.02v\). The inverse charged Euclidean propagator, continued below the cut, is

\[
\Gamma^{(2)}_{xx}(-s)=-s+b+2cv\,\delta v+cA(a)+4\lambda_XA(b)
+\alpha A(M_{D,0}^2)-\frac{g_3^2}{16\pi^2}B_0(s;a,b)+O(\hbar^2).
\]

The tadpole displacement is \(\delta v=-V_1'(v)/M_{\phi,0}^2\), with the same \(\overline{\mathrm{MS}}\) coefficients as the preceding report. It is not fitted to cancel a physical effect. The factor \(4\lambda_X\) includes both real X modes. The mixed \(\phi,x\) loop has two ordered internal contractions, canceling the usual one-half factor.

The bubble sign follows directly from the second variation of \(\tfrac12\mathrm{Tr}\ln H\): the two insertions give \(-\tfrac12\mathrm{Tr}(H^{-1}\delta H H^{-1}\delta H)\) in the inverse-propagator variation. In standard scalar self-energy notation, the same result follows from the negative \(B_{SS}=-B\) convention in Martin's equations (3.1)–(3.3), cited below. There is no positive-sign ambiguity after the inverse-propagator convention is fixed.

At strict first order, evaluate the self-energy at \(s=b\). Evaluating at \(s=0\) instead reconstructs the zero-momentum curvature. Therefore

\[
M_{X,\mathrm{pole}}^2-M_{X,\mathrm{curv}}^2
=-\frac{g_3^2}{16\pi^2}[B_0(b;a,b)-B_0(0;a,b)].
\]

At these exact benchmark masses,

\[
B_0(0)=1-2\ln2,\qquad B_0(b)-B_0(0)=1+\ln2-\frac\pi2.
\]

This independent analytic check fixes the sign and value without relying on the earlier effective-potential code.

## Charged threshold result

The separate integral reconstruction gives

\[
\delta v/v=1.21937398077\times10^{-4},\quad
M_{X,\mathrm{curv}}^2/v^2=0.01000277530641414,
\]

\[
\boxed{M_{X,\mathrm{pole}}^2/v^2=0.01000246538808106},\qquad
\boxed{M_{X,\mathrm{pole}}/v=0.1000123269404053+O(\hbar^2)}.
\]

The second boxed value uses the **strict first-order mass expansion**, \(\sqrt b+\delta M^2/(2\sqrt b)\). The square root of the truncated squared result is separately \(0.1000123261807317v\); their difference is induced higher-order content, not a two-loop calculation. The pole-minus-curvature squared shift is \(-3.09918333078\times10^{-7}v^2\).

The mixed pair threshold is \((\sqrt a+\sqrt b)^2=0.0582842712475v^2\), well above the evaluated shell \(s=0.01v^2\). The Feynman denominator is positive, so this scalar one-loop pole has no absorptive part. That statement is about the elementary vacuum particle, not a soliton lifetime. Its residue in the stated field convention is positive, \(Z_X=0.9999665711643\) at this order.

## All scalar species from full tensor contractions

As a second calculation, `neutral_poles.py` forms every cubic vertex \(g_{ijk}=\Lambda_{ijk\phi}v\) from the complete original real-field quartic tensor. For each species,

\[
\Pi_i(s)=\frac12\sum_k\Lambda_{iikk}A(M_{k,0}^2)
-\frac1{32\pi^2}\sum_{k,l}g_{ikl}^2 B_0(s;M_{k,0}^2,M_{l,0}^2),
\]

and the tree-Hessian displacement \(g_{ii\phi}\delta v\) is added once. All ordered pairs are counted explicitly. Both X components agree, and every zero-momentum curvature agrees with the earlier effective-potential calculation to the prescribed numerical tolerance.

| Particle | Tree mass squared / \(v^2\) | Strict one-loop pole squared / \(v^2\) | Strict one-loop pole mass / \(v\) |
|---|---:|---:|---:|
| \(\phi\) | 0.020000000000 | 0.0200095687556 | 0.141455186897 |
| X, both real components | 0.010000000000 | 0.0100024653881 | 0.100012326940 |
| D | 0.010100000000 | 0.0100998306010 | 0.100497913420 |

Every nonzero scalar cubic channel lies below its pair cut at the external tree shell. These are vacuum particle thresholds in the pure-new-scalar approximation. The earlier classical rotating-mode spectrum is not automatically loop-corrected by replacing its continuum thresholds with these numbers: the profile, conjugate frequency, kinetic terms and constrained quantum energy also require consistent corrections.

## Tests, running and limitations

The X-specific route passes ten checks. Adaptive Feynman-parameter quadrature and 32/64/128/256-node Gauss quadrature agree to a maximum \(5.0\times10^{-16}\) in B0. The analytic momentum derivative agrees with central differences to best relative \(8.6\times10^{-10}\). The full tensor route passes six further checks and independently reproduces X, while adding \(\phi,D\).

Running all ten already specified coefficients from the same matching point gives strict X pole squared values

| \(\mu/\mu_0\) | \(M_{X,\mathrm{pole}}^2/v^2\) |
|---:|---:|
| 1/2 | 0.0100024649954 |
| 1 | 0.0100024653881 |
| 2 | 0.0100024669387 |

The residual spread \(1.94\times10^{-9}v^2\) is an uncomputed-higher-order diagnostic, not a rigorous uncertainty. Solving the truncated inverse propagator implicitly yields \(s/v^2=0.0100024653056670\), differing from the strict answer by \(-8.24\times10^{-11}\); again this is induced higher-order content.

No finite tadpole, mass or gravitational vacuum constant was tuned. The previous raw scalar vacuum-energy shift remains; these particle-threshold calculations do not impose a zero gravitating vacuum density. Comparing a corrected vacuum pole to an **uncorrected classical** object energy would mix perturbative orders and cannot prove quantum binding. The renormalized fixed-charge quantum energy remains required.

## Reproduction and primary reference

Frozen inputs are in `pole_protocol.json` and `neutral_pole_protocol.json`; raw results and execution records are adjacent. Run `pole.py` before `neutral_poles.py`, using the same sibling layout as the repository. The initial X output hit a NumPy-boolean serialization error; the correction is recorded in `pole_correction_record.json`, with no physical formula or threshold change.

Stephen P. Martin, *Two-loop scalar self-energies in a general renormalizable theory at leading order in gauge couplings*, [hep-ph/0312092](https://arxiv.org/pdf/hep-ph/0312092), especially the one-loop scalar part of equations (3.1)–(3.3) and the propagator pole convention in section I. Only the one-loop scalar formulas are used here; no two-loop, gauge, fermionic or supersymmetric spectrum is imported.
