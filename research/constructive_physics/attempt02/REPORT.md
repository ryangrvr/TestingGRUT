# CP-02: a protected interaction relation and its fundamental boundary

**AUTHOR RESULT — conditional geometric EFT survives; fundamental GRUT law not
established.** One authorized construction, D003, frozen at
**d88a1a0bb0d4edea54b075437f0edd0921381df2** after disclosed analytic design and
before executable scientific checks. Starting handoff:
466c3d9769c3bc53c10662d2b01146e194fc7208.

CP-02 derives a relation between two independently definable observables:
the soft scalar self-scattering coefficient kappa and the long-range coefficient
C_7 of a scalar-induced force between matter sources:

\[
\boxed{C_7=\frac{3}{8\pi^3}\kappa^2.}
\]

Its load-bearing premise is that scalar motion and matter inhabit the same
induced geometry of a timelike hypersurface in flat five-dimensional spacetime.
The relation follows from that physical symmetry and ordinary quantum theory,
not from fitting a force to scattering. Quantum corrections generate higher
operators but cannot independently change these leading coefficients while
the stated symmetry, spectrum and matching assumptions hold.

This is not a novel restriction over the full conventional theory with those
same premises: ordinary brane EFT already imposes it. The construction also
cannot recover Einstein gravity about its planar vacuum. No experimentally
established scalar, extra dimension, new GRUT law or theory of reality is claimed.

## 1. The law, its physical price and correspondence domain

The frozen LAW.md gives the complete EFT domain and counterterm family. Here
is its lowest-derivative action, with eta_AB=diag(-1,1,1,1,1):

\[
\begin{split}
\gamma_{\mu\nu}&=\eta_{AB}\partial_\mu X^A\partial_\nu X^B,\\
S_0&=-T\int d^4x\sqrt{-\gamma}
-\frac12\int d^4x\sqrt{-\gamma}
 \left(\gamma^{\mu\nu}\partial_\mu\chi\partial_\nu\chi+\mu^2\chi^2\right).
\end{split}
\]

X^A embeds a four-dimensional surface. In static gauge X^A=(x^mu,Y), the
canonical motion field is pi=sqrt(T_R)Y:

\[
\gamma_{\mu\nu}=\eta_{\mu\nu}+\frac{\partial_\mu\pi\partial_\nu\pi}{T_R}.
\]

T_R is renormalized, positive and independently measurable. This construction
does not inherit CP-01's potential, massive scalar or quartic/mass relation.
The motion field is massless because ambient translations forbid a Y potential.

The quantum EFT includes every local invariant made from the induced metric,
extrinsic curvature K_mu nu, covariant derivatives and matter fields, with the
frozen reflection symmetries. Coefficients of R, R^2, Ricci-squared, tr(K^4),
chi^2 R, matter self-interactions and higher invariants are not silently removed.
R=K^2-K_mu nu K^mu nu by the Gauss equation. Redundant operators can be traded
by integration by parts and perturbative field redefinitions, not discarded
to force a result. This is a derivative expansion, not a renormalizable
two-parameter fundamental action.

| Structure | Status and price |
|---|---|
| Flat ambient dimension, signature and one embedding normal | POSTULATED physical inventory |
| Actual ambient ISO(4,1), worldvolume covariance, normal reflection | POSTULATED selection/protection premise |
| Matter and probe coupling to the same induced metric | POSTULATED universal coupling within the chosen sector |
| Quantum preparations, probabilities, local quantization and composition | POSTULATED established operational structure |
| T_R, matter masses, cutoff and all finite matching data | SUPPLIED parameters; unbounded operator inventory at arbitrary order |
| Vacuum, small-data domain, source history and measurement procedures | SUPPLIED boundary/operational premises |
| Leading self-scattering and source-coupling relation | DERIVED conditional consequences |
| Physical necessity of the geometry; UV completion; observers and consciousness | NOT DERIVED |

No formal Appendix-B/L_0 price, compression credit or Stage-3 score is issued.
A short geometric formula does not erase the supplied dimensional, symmetry,
matter or matching information. There is one CP-02 analytic candidate draft;
CP-01 had two. The lifetime constructive design log therefore records three
analytic drafts and two selected mechanisms. This is not Stage-3 Card 3.

For energies much smaller than T_R^(1/4), the cutoff and relevant heavy gaps,
the flat-vacuum kinetic theory reduces to a canonically normalized massless
scalar and ordinary relativistic matter. Ordinary quantum mechanics is used
to quantize them; its Hilbert space or Born rule is not derived. Einstein
gravity is not supplied or recovered. Actual finite-G apparatus predictions
would need additional standard gravitational terms and calibration; none is
certified here.

## 2. Dynamics, physical symmetry and stability

For the displayed leading action define matter stress by
delta S_m=(1/2)int sqrt(-gamma) T_m^mu nu delta gamma_mu nu. Variation gives

\[
\partial_\mu\!\left[\sqrt{-\gamma}
 (T\gamma^{\mu\nu}-T_m^{\mu\nu})\partial_\nu X_A\right]=0,
\qquad (\Box_\gamma-\mu^2)\chi=0.
\]

Tangential equations are the diffeomorphism identities, including matter
stress conservation on its equations. The normal equation is
(T gamma^mu nu-T_m^mu nu)K_mu nu=0. Higher operators add their ordinary
Euler-Lagrange variations and are treated perturbatively. Translation and
Lorentz invariance supply conserved total ambient Noether momentum and angular
momentum. Matter can exchange energy/momentum with surface motion; the matter
stress alone need not be conserved relative to the background eta.

The static-gauge nonlinear boost is

\[
\delta Y=b^\mu(x_\mu+Y\partial_\mu Y),\qquad
\delta\chi=\xi^\mu\partial_\mu\chi,\qquad \xi^\mu=Yb^\mu.
\]

It induces delta gamma=L_xi gamma. Consequently the matter action transforms
as an induced scalar density. This is an active constraint relating physically
different embeddings, not independent active local-unitary invariance. It
does not imply frozen local marginals, the obstruction that killed Card 2.
Passive changes of coordinates are separately handled by covariance.

For a first-gradient scalar density f(z), z=(partial Y)^2, the boost variation
is a total derivative plus
[2(1+z)f'(z)-f(z)] b^mu partial_mu Y. Modulo an irrelevant additive
constant, invariance imposes

\[
2(1+z)f'(z)=f(z),\qquad f(z)=-T\sqrt{1+z}.
\]

The additive constant is a vacuum subtraction in static gauge, not a new
interaction. A constant-gradient embedding has no extrinsic curvature:
no other invariant can alter this first-gradient sector. Expanding in the
canonical field,

\[
\mathcal L_\pi=-T_R-\frac12(\partial\pi)^2
+\frac{[(\partial\pi)^2]^2}{8T_R}
-\frac{[(\partial\pi)^2]^3}{16T_R^2}+\cdots.
\]

For the pure leading scalar, with momentum P conjugate to pi,

\[
\mathcal H=T_R
\sqrt{\left(1+\frac{(\nabla\pi)^2}{T_R}\right)
      \left(1+\frac{P^2}{T_R}\right)}
\]

is positive; subtract its vacuum T_R for excitation energy. A homogeneous
background dot Y=v, |v|<1, has characteristic speed squared 1-v^2.
Induced-null tangent vectors satisfy eta_mu nu v^mu v^nu=-(partial Y dot v)^2,
so they are not outside the ambient Minkowski cone. These are controlled
leading-sector checks. They do not establish stability or causality for every
large-curvature background or arbitrary resummation of Wilson coefficients.
Regular small data admit a local hyperbolic evolution; global graph existence
is unproved and unclaimed.

There is no Hamiltonian depending on how a density matrix was decomposed into
preparations. At a regulated, ordinary quantum level,
E(rho)=U rho U^dagger is affine on mixtures. Perturbative local QFT supplies its
usual commutator/causality framework in the weak domain. An all-scale
nonperturbative quantum completion is not proved by this observation.

## 3. Deriving the joint observable

Four differentiations of the quartic vertex give the massless tree amplitude

\[
\mathcal M_4=\frac{s^2+t^2+u^2}{4T_R},\qquad s+t+u=0.
\]

Here eta=(-,+,+,+), s=-(p_1+p_2)^2 and all momenta are incoming.
The coefficient kappa is defined by scaling all external momenta to zero at
fixed scattering angle after canonical normalization:

\[
\kappa=\lim_{\varepsilon\to0}
\frac{\mathcal M_4(\varepsilon p_i)}
     {s(\varepsilon)^2+t(\varepsilon)^2+u(\varepsilon)^2}
=\frac1{4T_R}.
\]

This limit isolates a coefficient; it is not a claim that arbitrary
finite-energy amplitudes have no corrections.

Expanding the same induced matter metric yields

\[
S_{\pi m}=\frac1{2T_R}\int d^4x\,
T_m^{\mu\nu}\partial_\mu\pi\partial_\nu\pi+\cdots.
\]

For a static probe S_i=-M_i int dt sqrt(1-dot pi^2/T_R).
After Euclidean rotation its leading vertex is
C_i int d tau (partial_tau pi)^2, with C_i=M_i/(2T_R).
There is no one-pi source vertex and no classical static scalar field profile.
The first connected force comes from exchange of two virtual massless scalars.

The free Euclidean propagator is

\[
D(\tau,r)=\frac1{4\pi^2(\tau^2+r^2)}.
\]

Wick contraction gives

\[
\left\langle(\partial_\tau\pi)^2_{1}
(\partial_\tau\pi)^2_{2}\right\rangle_c
=2[D''(\tau,r)]^2
=\frac{(3\tau^2-r^2)^2}{2\pi^4(\tau^2+r^2)^6}.
\]

The connected vacuum energy is negative:

\[
V_\pi(r)=-C_1C_2\int_{-\infty}^{\infty}d\tau\,2[D''(\tau,r)]^2
=-\frac{3M_1M_2}{128\pi^3 T_R^2r^7}.
\]

For completeness, after tau=r x the exact moments over the real line are
int dx x^(2k)/(1+x^2)^6 = pi {63,7,3}/256 for k=0,1,2.
Their combination 9 I_2-6 I_1+I_0=3pi/16 gives the displayed coefficient.
Independent numerical improper quadrature checks that integral.

Eliminating the independently calibrated T_R yields

\[
C_7=-\lim_{r\to\infty}\frac{r^7 V_\pi(r)}{M_1M_2}
=\frac3{128\pi^3T_R^2}
=\frac3{8\pi^3}\kappa^2.
\]

The force coefficient is defined at leading order in weak-source strength.
The r->infinity limit follows that expansion in the frozen flat vacuum.
It is not the limit of a real apparatus's total potential with other forces
left unsubtracted. It is not a tautology: the observables come from separate
scattering and source-response operations; an ordinary 4D EFT can vary their
ratio.

## 4. Quantum protection, counterterms and RG

The protection statement is conditional on an exact, anomaly-free realization
of the postulated ambient symmetry and symmetry-preserving matching. In
geometric variables the quantum effective action obeys its Ward identity.
Its local zero-curvature, zero-matter density is again a volume term. Heavy
matter can renormalize its tension; defining pi=sqrt(T_R)Y renormalizes
self-interactions and the induced matter metric together. Normal reflection
forbids odd-K terms; curvature and derivative invariants do not change the
leading four-point p^4 coefficient or leading zero-transfer probe vertex.

Light-pi loops carry additional powers of momentum. They give nonlocal logs
at higher order, not a new leading p^4 coefficient. Insertions of geometric
higher operators also vanish more rapidly in the soft/long-distance limits
for finite matching coefficients and the stated gapped matter spectrum.
No additional massless species are admitted. Thus the leading relation is
stable in this perturbative EFT, although finite-energy observables remain
unselected.

The self-loop power-counting statement uses dimensional regularization.
In a cutoff description power-law matching corrections can renormalize
the tension and its derivative tower. The Ward identity ties that tower
and the canonical normalization together; it does not say T itself never
renormalizes. The observable relation is expressed in physical low-energy
coefficients, not in an unrenormalized parameter.

This is an all-order structural Ward argument with supplied symmetry, not an
all-order cancellation deduced from one successful loop. For the selected
bosonic inventory a covariant regulator and counterterms can preserve the
embedding symmetry perturbatively; no chiral gauge anomaly is introduced.
This does not prove an anomaly-free unknown UV theory, gravitational
completion or validity when extra light fields or explicit breaking appear.

### A concrete one-loop calculation

In the physical s-channel write t=-s(1-z)/2, u=-s(1+z)/2. Then

\[
\mathcal M_{\rm tree}(s,z)=\frac{s^2(3+z^2)}{8T_R}.
\]

With identical intermediate scalars,
2 Im M_s=(1/2)int dPhi_2 M_in M_out and dPhi_2=dOmega/(32pi^2).
The angular moment
int dOmega(3+a^2)(3+b^2)=4pi(166+2z^2)/15 gives

\[
\operatorname{Im}\mathcal M_s^{(1)}
=\frac{s^4(83+z^2)}{15360\pi T_R^2}.
\]

An independent partial-wave calculation, using
M=32pi sum_l(2l+1)a_l P_l for identical bosons, has
a_0=5s^2/(384pi T_R), a_2=s^2/(1920pi T_R) and reproduces the cut.

In d=4-2epsilon the logarithmic contribution is

\[
\mathcal M^{(1)}_{\log}=
-\frac1{15360\pi^2 T_R^2}\sum_{\rm cyc}
\left[83s^4+s^2(t-u)^2\right]
\log\frac{-s-i0}{\mu_R^2}.
\]

Massless tadpoles are scaleless in this scheme. The identity
sum_cyc s^2(t-u)^2=s^4+t^4+u^4 fixes the UV pole:

\[
\mathcal M_{\rm div}^{(1)}
=\frac7{1280\pi^2 T_R^2\epsilon}(s^4+t^4+u^4).
\]

For the physical on-shell local coefficient
M_local^(8)=c_8(mu_R)(s^4+t^4+u^4)/T_R^2, this requires

\[
\delta c_8=-\frac7{1280\pi^2\epsilon},\qquad
\frac{dc_8}{d\log\mu_R}=-\frac7{640\pi^2}.
\]

The code explicitly verifies cancellation of renormalization-scale dependence
against the logs. The beta function applies to this on-shell coefficient in
the displayed pure-pi order, not to an individual off-shell operator in every
basis. Higher-order/matter terms alter their appropriate higher projections.
The finite value c_8(mu_0) remains independently matchable. Small loops do not
bound it.

At generic nonzero s,t,u the derivative-coupled massless bubble has no separate
infrared pole in this calculation. Its p^8 log(p^2) term divided by the
leading p^4 polynomial vanishes in the fixed-angle soft limit. This does
not certify every inclusive experimental cross-section: real unresolved
emission, acceptance and any added massless sector must be treated at the
order used by an actual experiment. In particular, finite-G graviton
infrared effects have not been computed. Local analytic contact terms in
the source potential do not change its separated-source r^-7 coefficient;
higher derivative terms change subleading tails and finite-distance errors.

**Source normalization issue preserved:** arXiv:2010.05913v1 Eq.(4.20) prints
7/10240 instead. Our explicitly normalized unitarity result agrees with
our 4! contraction of the same paper's geometric Eq.(4.17): its pole density
7 tr(K^4)/(320pi^2 epsilon) has contact contraction
(s^4+t^4+u^4)/(4T_R^2), giving 7/1280. Independently, its N=1,
L->infinity Eq.(5.20), summed over crossings, reduces to 7/1280 as well.
The executable exact polynomial reduction records the factor-eight mismatch.
We have not obtained author confirmation or checked an erratum; this is an
explicit external-review obligation, not a silently corrected quotation.

### Massive matter cannot be omitted from the inventory

For the minimally coupled real chi on a background gamma, the Euclidean
determinant is (1/2)Tr log(-Box_gamma+mu^2). Modulo total derivatives its pole is

\[
\Gamma_{\chi,\rm pole}=
-\frac1{32\pi^2\epsilon}\int d^4x\sqrt{\gamma}
\left[\frac{\mu^4}{2}-\frac{\mu^2R}{6}
+\frac{R^2}{72}+\frac{R_{\mu\nu\rho\sigma}^2-R_{\mu\nu}^2}{180}\right].
\]

This follows by expanding exp(-mu^2 t) in the standard scalar heat kernel.
It requires tension, induced-curvature and curvature-squared counterterms.
Its sign is in the Euclidean determinant convention; it is not the Lorentzian
scattering sign above. Allowed nonminimal matter terms require their own
matching and change higher curvature coefficients. This calculation
demonstrates why an exact finite operator family is not claimed. It preserves
the common metric rather than generating an independent leading disformal
coefficient. Full mixed-loop, high-field and UV calculations have not been done.

In minimal subtraction the displayed massive field gives
beta_T=mu^4/(32pi^2). Its finite zero-curvature contribution is
mu^4[log(mu^2/mu_R^2)-3/2]/(64pi^2); its scale derivative cancels beta_T.
After that volume contribution is included and the field recanonicalized,
T_R denotes the physical vacuum effective tension in the soft prediction.
Using a running bare/MS tension in the tree formula without this matching
would give an artificial scale dependence. The code checks the exact
scale-cancellation coefficients.

## 5. The strongest conventional comparison

### What the extra physical premise actually removes

Ordinary four-dimensional Lorentz, shift and reflection invariance permits
the same scalar self-scattering coefficient with an independent matter metric

\[
\widetilde\gamma_{\mu\nu}
=\eta_{\mu\nu}+\frac{b}{T_R}\partial_\mu\pi\partial_\nu\pi.
\]

In the same weak domain it has C_7(b)=b^2 C_7(1) and unchanged kappa.
The declared b=0,1,2 cases yield force ratios 0,1,4 with the same low-energy
field inventory, masses, vacuum and self-scattering calibration. These are
ordinary local perturbative EFTs; positive kinetic spectra and weak coupling
do not force b=1. Their full UV completions are not established, just as
CP-02's is not.

b=2 is a genuine conventional counterexample to a relation asserted solely
from four-dimensional shift/Lorentz quantum physics. It **does not** share the
additional nonlinear five-dimensional physical symmetry. That premise is
exactly what buys the relation; omitting this distinction would be an unfair
comparison.

### No exclusion beyond the complete same-symmetry baseline

Include the full actual CP-02 symmetries and shared induced geometry in the
baseline. Its conventional EFT already has the same Ward argument and
leading relation. No admissible same-premise model violating that asymptotic
relation has been exhibited. No new GRUT restriction over this full baseline
is established. A model need not have all Wilson coefficients fixed to make
a leading observable relation rigid, but the rigidity here is established
brane physics under an extra, unverified physical hypothesis.

### A same-symmetry counterexample to finite-energy rigidity

Add zeta int sqrt(-gamma)R^2, an allowed geometric term. At four-pi order and
on shell,

\[
\Delta\mathcal M_\zeta
=\frac{\zeta}{2T_R^2}(s^4+t^4+u^4).
\]

At a right-angle massless collision with each particle energy E,
s=4E^2 and t=u=-2E^2:

\[
\mathcal M_0=\frac{6E^4}{T_R},\qquad
\frac{\Delta\mathcal M_\zeta}{\mathcal M_0}
=24\zeta\frac{E^4}{T_R}.
\]

For the predeclared zeta=1 controls this is 3/20000 at E/T_R^(1/4)=1/20
and 3/1250 at 1/10. These exact conditional differences are software
benchmarks, not observations. The allowed R term produces an independent
p^6 amplitude already, so R^2 is not the only uncertainty in a real
finite-energy forecast. We match the other terms equal only to isolate
this explicit witness; we do not set them universally to zero.

The EFT remains perturbative in these small fixtures; higher operators are
not resummed into new light ghosts. A positive R^2 coefficient can also be
obtained by integrating out a healthy heavy field linearly coupled to R;
this illustrates an ordinary matching origin, not an addition to CP-02's law
or evidence for a particular UV completion. Both comparators share the
leading observables but differ at finite energy. A demand for exact
finite-energy rigidity is therefore refuted even by the full symmetry.

## 6. Attempting to recover gravity: precise failure

About Y=0 the metric perturbation is

\[
\gamma_{\mu\nu}-\eta_{\mu\nu}
=T_R^{-1}\partial_\mu\pi\partial_\nu\pi,
\qquad \delta\gamma_{\mu\nu}|_{\pi=0}=0.
\]

There is one transverse massless scalar and no linear independent spin-2
metric mode. Coupling to matter begins at pi^2, so no tree single-particle
exchange creates a Newtonian static 1/r potential.

An allowed induced Einstein-Hilbert term cannot repair that vacuum:

\[
\left.\sqrt{-\gamma}R\right|_{\pi^2}
=\frac1{T_R}\left[(\Box\pi)^2-
(\partial_\mu\partial_\nu\pi)^2\right].
\]

Its spacetime integral is a boundary term. The Fourier quadratic kernel
(k^2)^2-k_mu k_nu k^mu k^nu vanishes identically. It does not create a
new propagating graviton. The nonzero scalar quantum force is r^-7,
not r^-1.

For a nonplanar background a normal displacement gives a metric variation
proportional to K_mu nu delta Y, plus a tangential diffeomorphism. There is
still only one normal field. The tensor shape is supplied by the chosen
embedding curvature; it does not provide two independent Lorentz-invariant
spin-2 polarizations around the flat vacuum. Curved-background effects are
not a rescue of the required observed flat-space gravity limit.

**Failed stronger claim:** the selected codimension-one induced geometry,
even with an induced R term, generates Einstein gravity about its declared
vacuum. **Counterexample/obstruction:** its exact linearized field inventory
and vanishing EH quadratic kernel. Adding a bulk metric, more embedding
normals or a new principle would change the law. None is introduced here.

## 7. Operational test and uncertainty

A prospective conditional test first measures the massless scalar's dispersion
and self-scattering at several low energies/angles, extracting kappa with
independently bounded higher-derivative coefficients. Separately measure source
masses and the connected scalar-induced potential at several large separations
in a matched vacuum. Form the dimensionless ratio

\[
Q=\frac{8\pi^3 C_7}{3\kappa^2}.
\]

CP-02 predicts Q=1; the b=2 ordinary 4D comparator predicts Q=4.
Neither T_R nor the matter coupling is refitted using that force holdout.
Mass calibration and any other field exchange/background are independent.
The amplitude sign can be fixed by interference or otherwise calibrated;
its square alone suffices for Q.

A nonzero discrepancy beyond separately certified uncertainties in soft
extrapolation, source size, recoil, extra light fields, vacuum preparation,
ordinary electromagnetic/gravitational backgrounds and apparatus calibration
would refute the **combined CP-02 hypotheses**. It would not establish a
different fundamental theory automatically.

Writing R_Q=log Q, an independent small-error budget is
sigma_RQ^2 approximately (sigma_C7/C7)^2+4(sigma_kappa/kappa)^2;
correlated calibration requires its covariance term, not this formula.
A two-sided z-sigma discrimination from b=2 would require
z sigma_RQ plus a certified bound on total log-systematic error below log 4.
No measurements or bounds achieving that inequality are supplied. At finite
E,r unrestricted matching coefficients can destroy predictive power.
No numerical cutoff or naturalness estimate is treated as a certified error.

There is no detected CP-02 massless scalar or apparatus with both calibrated
scalar self-scattering and isolated r^-7 force. Existing constraints on
disformal interactions are model- and spectrum-dependent; historical values
are not a current apparatus certificate. Experimental accessibility is
**NOT ESTABLISHED**. A genuinely asymptotic limit is not an experiment at
infinite distance with infinite resolution. Finite-G contributions need a
controlled extension and comparator before a laboratory prediction.

## 8. Novelty and primary-source audit

This packet does real calculations, but its mechanism is not original GRUT
physics. The assembled leading relation is an inference combining standard
DBI self-interactions and induced-metric disformal interactions; no priority
claim is made.

| Primary source | Material actually inspected | Role |
|---|---|---|
| Goon, Hinterbichler, Trodden, arXiv:1103.5745 | Brane construction and symmetry sections | Established geometric/symmetry framework |
| Goon et al., arXiv:1606.02295 | Nonrenormalization discussion for derivative theories | Established protection context; not a new GRUT identity |
| Goon, Melville, Noller, arXiv:2010.05913v1 | Explicit PDF/text: Secs.4.1–4.2, Eqs.4.17–4.20, Sec.5.2 Eqs.5.19–5.20, App.D | Covariant brane loops and explicit normalization cross-check |
| Brax, Burrage, arXiv:1407.1861 | Secs.2 and 4–6, Eqs.2.5 and 6.1, Sec.4 footnote | Published disformal force; footnote already states M^4=2 f^4 |
| Vassilevich, arXiv:hep-th/0306138 | Scalar heat-kernel formula, Eq.4.28 | Established matter determinant coefficients |

URLs and access scopes are recorded in SOURCES.json. We do not claim complete
page-by-page reviews of those papers. The Brax–Burrage convention
(1/M^4)T_m^mu nu partial_mu pi partial_nu pi has M^4=2T_R;
its -3M_1M_2/(32pi^3 M^8 r^7) matches our coefficient exactly.
That source already makes the brane-coupling connection. The present
elimination into kappa is not evidence of a new fundamental relation.

Ordinary QFT, nonlinear realizations of spacetime symmetry, brane EFT and
Casimir-Polder-type virtual exchange cover the result. Nothing here derives a
preferred measurement, conscious observer, persistent partition, or new
gravitational back-reaction rule. Renaming the fields would not change that.

## 9. Computation and verdict

Run:

~~~bash
python research/constructive_physics/attempt02/check.py
~~~

Verified output:

~~~text
CP02_GEOMETRIC_EFT_CHECKS_PASS
groups=21 max_numeric_error=5.329e-15
CONDITIONAL_EFT_ONLY; NO_NEW_FUNDAMENTAL_GRUT_LAW; NO_NEWTONIAN_GRAVITY_LIMIT
~~~

The 21 assertion/reporting groups include freeze verification, exact vertex
permutations, Ward recursion, improper-integral moments, countermodels, two
independent unitarity routes, whole-polynomial crossing identities, a primary
source normalization audit, RG scale cancellation, heat-kernel inventory,
leading stability, the exact gravity obstruction and the ordinary mixture
boundary. There are 23 numerical comparisons. Not every group is an
independent physical test: the heat-kernel inventory uses the established
coefficient formula, and the mixture statement explicitly uses supplied
ordinary quantization. Requirements pin NumPy 2.3.5 and SciPy 1.17.0;
the recorded Python is 3.12.14. All results are software/author evidence.
The execution log also discloses routine solver settings, two Fourier diagnostic
vectors and the RG reference c_8=0 that were set in code before execution but
not individually enumerated in the frozen design ledger. They are recorded
after execution without backdating. No new candidate or parameter-survivor
selection occurred.

**Derived:** the common geometry's leading interaction coefficients, quantum
two-source force, asymptotic joint relation, pure-pi one-loop cut/pole/RG,
same-symmetry finite-energy countermodel and flat-vacuum gravity obstruction.

**Postulated:** the extra geometry and physical symmetry, shared matter
coupling, quantum operational rules, field inventory, matching and vacuum.

**Falsified:** exact finite-energy rigidity under the full allowed geometric
EFT, and a Newton/Einstein-gravity limit generated by this single embedding
field about its frozen vacuum. The hypothesis that nature has this sector
has not been experimentally tested.

**Author verdict: CONDITIONAL ONLY.** CP-02 is alive as a controlled
symmetry-based EFT for its declared sector. It fails promotion to an accurate
fundamental theory of reality: physical necessity, novelty beyond the complete
same-premise conventional baseline, observed gravitational correspondence
and experimental support are missing. There is no independent validation,
banking or kit-certified pass.

The consequential result for GRUT is specific: quantum consistency can preserve
a leading observable relation while leaving infinitely many higher matching
directions free. What remains unjustified is why nature realizes the protecting
geometry and universal matter coupling, and how observed spin-2 gravity would
arise. This packet does not formulate a replacement answer.

One CP-02 attempt is accounted for and reaches this substantive boundary.
CP-03 and Card 3 remain unopened; the broad GRUT mission remains active.
Independent review must check the Ward/matching scope, source-force
normalization, UV-pole discrepancy and gravity argument. PR #12 remains
draft; all previous frozen packets remain preserved.

## Packet files and preservation

New files, all within research/constructive_physics/attempt02:
AUTHORIZATION.md, DESIGN_LOG.json, LAW.md, FREEZE.json, REPORT.md, check.py,
results.json, requirements.txt, SOURCES.json, EXECUTION_LOG.json,
DISPOSITION.json and MANIFEST.json. The first four were published in the
pre-execution freeze; the remaining eight complete the scientific packet.
Only the current research README.md, REORIENTATION_STATE.json and CHECKS.md
are updated outside this directory. The extended handoff is not altered.

All 101 previously manifested scientific files and the CP-01 freeze's
additional top-level files were verified byte-identical. The CP-02 law and
design freeze are also unchanged. The implementation ran twice: the first
complete run and a second run after adding the massive-tension
scale-cancellation check; no law or parameter search was restarted.
