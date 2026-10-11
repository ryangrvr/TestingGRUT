# CP-03: gravity from a scalar determinant — worked construction and obstruction

**Author verdict: NOT YET ESTABLISHED as physical gravity or a new GRUT law.**

A scalar composite mechanism can reproduce the massless gravitational exchange
tensor without an Einstein-Hilbert term at the start. We independently reproduce
the leading loop normalization, two-helicity residue, common scalar coupling
and conditional Newton potential. This is a published mechanism, not a GRUT
discovery. Its constraints are also derived below, rather than inferred from
the existence of a pole.

The selected clockless action has an exact real-classical signature obstruction.
The formal quantum calculation requires a Lorentzian continuation, a measure
and UV physics that this construction does not establish. Its supplied
reparametrization invariance does not independently generate the full Einstein
constraint dynamics or the physical observable algebra. Those are substantive
failures of the requested fundamental-origin claim.

Freeze: `d197e6ac8ad91e8d14ab82eca8c34b53add15017`. One selected mechanism,
D004, CP-03. The freeze followed disclosed analytic derivation and literature
inspection and preceded executable scientific checks. It was not blind law
selection. The owner's uploaded note about constraints and subsystem structure
was incorporated before freeze. No change to the selected action followed it.

## 1. THE LAW

The full frozen specification is [LAW.md](LAW.md). In four dimensions, with
natural units and N real scalars,

\[
 Q_{\mu\nu}=\sum_{a=1}^{N}\partial_\mu\phi^a\partial_\nu\phi^a,
 \qquad V=V_0+\frac{m^2}{2}\sum_a(\phi^a)^2,
 \qquad S=\int d^4x\,\frac{\sqrt{|\det Q|}}{V}.
\]

The differentiable real branch has rank Q=4 and V nonzero. Variations are
compactly supported. The construction does not prove global existence across
rank-changing configurations. Positive V0 and m, a positive field-space metric,
four dimensions, the manifold and scalar inventory are supplied. No metric,
clock fields, graviton or Einstein-Hilbert term appears in this action.

This is the clockless determinant mechanism described in Carone's
[2019 review, Section IV](https://arxiv.org/pdf/1911.13166). The scalar
composite-gravity calculation was introduced by Carone, Erlich and Vaman in
[2016](https://arxiv.org/pdf/1610.08521); its important sign corrections and
cubic/self-coupling analysis are in
[Carone, Claringbold and Vaman](https://arxiv.org/pdf/1710.09367v2).
This packet reconstructs and tests that mechanism; it does not claim its
invention or a proof of its complete quantum theory.

The formal quantum organization uses a large-N expansion, covariant
dimensional regularization D=4-epsilon, normal-ordering/tadpole conditions and
an expansion about a flat metric. Adding and subtracting the same tensor is an
algebraic organization of the scalar action, not evidence that real Lorentzian
spacetime has been generated. The constrained measure and an unconstrained
auxiliary-metric integral need not coincide at the quantum level.

**Supplied assumptions and information price.** The following ten categories
are load-bearing: (1) manifold, differentiability and dimension; (2) real scalar
inventory and positive target metric; (3) the determinant action; (4) the
particular potential and its parameters; (5) reparametrization invariance;
(6) ordinary quantum preparations/probabilities; (7) large N and the prescribed
loop organization; (8) regulator and measure assumptions; (9) flat-vacuum,
tadpole and masslessness conditions; (10) Lorentzian reconstruction and
asymptotic states when interpreting a scattering amplitude. Items 8 and 10 are
not completed mathematical constructions. A count of categories is not a bit
price. No Stage-3 L0/Appendix-B score or compression credit is claimed. This
scope does not attempt Stage 3's observer/carrier/interface chain.

## 2. THE DERIVATION

### Algebraic metric and the exact signature obstruction

For classical elimination only, introduce an auxiliary metric with action

\[
 S_{\rm aux}=\int\sqrt{|g|}\left(\frac12g^{\mu\nu}Q_{\mu\nu}-V\right)d^Dx.
\]

Its metric variation gives

\[
 T^{\rm full}_{\mu\nu}=Q_{\mu\nu}-g_{\mu\nu}
 \left(\frac12 g^{\rho\sigma}Q_{\rho\sigma}-V\right)=0.
\]

Take its trace, solve the scalar factor and obtain, for D>2,

\[
 g_{\mu\nu}=\frac{D/2-1}{V}Q_{\mu\nu}.
\]

Substitution yields the frozen determinant action at D=4. **DERIVED:** an
algebraic composite metric. **POSTULATED:** the variational law producing it.

For every real tangent vector v,

\[
 v^\mu Q_{\mu\nu}v^\nu
 =\sum_a\left(v^\mu\partial_\mu\phi^a\right)^2\ge0.
\]

At rank four this is positive definite. Multiplication by a nonzero scalar
can make it positive or negative definite, never Lorentzian. Inertia is also
unchanged by a real coordinate transformation. Thus:

\[
 \boxed{\text{The frozen positive-Gram action has no nondegenerate real
 Lorentzian classical metric branch.}}
\]

This is a proof for the entire stated real classical family, not a sample
of gradient matrices. It is not a proof that a Euclidean quantum theory could
never admit a Lorentzian reconstruction. That requires an appropriate measure,
positivity and continuation, none proved here. A renormalized composite
operator is not automatically the same object as a pointwise classical Gram
matrix. We do not promote the classical theorem into a quantum no-go result.

### Constraints actually produced by the action

Write G_mu,a=partial_mu phi^a, L=sqrt(det Q)/V and use a coordinate x0 only as a
canonical bookkeeping parameter. Then

\[
 p_a=L(Q^{-1})^{0\nu}\partial_\nu\phi^a,
 \qquad p_a\partial_\mu\phi^a=L\delta^0_\mu.
\]

Consequently the canonical Hamiltonian vanishes and the three momentum
constraints hold. With q_ij=partial_i phi dot partial_j phi, U(phi)=1/V(phi),

\[
 H_i=p_a\partial_i\phi^a=0,
 \qquad H_0=\frac12\left[p_ap_a-U^2\det q\right]=0.
\]

These follow from the selected action. They are inherited from its supplied
reparametrization symmetry; that symmetry has not itself emerged from a more
primitive dynamics. On a regular branch the velocity Hessian is

\[
 \frac{\partial p_a}{\partial\dot\phi^b}
 =L(Q^{-1})^{00}\left[\delta_{ab}
       -(G^TQ^{-1}G)_{ab}\right].
\]

The bracket is a target-space normal projector of rank N-4. Four gauge
directions are removed. N=4 is locally devoid of ordinary embedding degrees
of freedom on this regular branch; N=6 having two classical embedding degrees
does not make them graviton helicities. Those are different statements.

For the canonical Poisson bracket, smearing H0 with density-weight-minus-one
functions n,m gives exactly

\[
 \{H_0[n],H_0[m]\}
 =D\!\left[-U^2\det(q)q^{ij}
             (n\partial_jm-m\partial_jn)\right].
\]

The derivation uses

\[
 \frac{\partial H_0}{\partial p_a}=p_a,
 \quad
 \frac{\partial H_0}{\partial(\partial_i\phi^a)}
 =-U^2\det(q)q^{ij}\partial_j\phi^a.
\]

Terms without derivatives of the smearing functions cancel by antisymmetry,
including terms from U(phi). The remaining coefficient is precisely the
spatial constraint. With H=H0/[U sqrt(det q)], its normalized weak bracket is

\[
 \{H[N],H[M]\}\approx
 -D\!\left[q^{ij}(N\partial_jM-M\partial_jN)\right].
\]

The negative sign is the Euclidean normal-deformation sign for this positive
target metric and these generator conventions. It is consistent with the
signature theorem. **DERIVED:** a worldvolume deformation algebra.
**NOT DERIVED:** Lorentzian ADM dynamics for an independent spatial metric,
its quantum constraint algebra, or its physical Hilbert space. An embedding
theory can realize a deformation algebra without being Einstein gravity.

### Leading interaction and independently calculated scalar loop

Expand the determinant and potential. At four dimensions, the interacting
part can be written

\[
 \mathcal L_{\rm int}
 =-\frac{1}{4V_0}\left(T_{\mu\nu}T^{\mu\nu}-\frac12T^2\right)
       +O(V_0^{-2}).
\]

Here T is the free matter stress tensor, **not** the full constrained stress
tensor. We check the equality by two independent exact routes: direct
permutation expansion of the determinant and contraction of the stress tensor.

For equal-mass scalar states, the stress insertion is

\[
 E_{\mu\nu}(p_1,p_2)=-(p_{1\mu}p_{2\nu}+p_{2\mu}p_{1\nu})
       +\eta_{\mu\nu}(p_1\cdot p_2+m^2),
 \qquad q^\mu E_{\mu\nu}=0.
\]

Introduce a Feynman parameter x, z=x(1-x), b=x-1/2. After Wick rotation and
momentum shift, the Euclidean loop insertion becomes

\[
 E_{\mu\nu}=2k_\mu k_\nu-\delta_{\mu\nu}(k^2+m^2)
 +2b(k_\mu Q_\nu+Q_\mu k_\nu-\delta_{\mu\nu}k\cdot Q)
 +z(\delta_{\mu\nu}Q^2-2Q_\mu Q_\nu),
\]

with denominator [k^2+Delta]^2 and Delta=m^2+zQ^2. The radial integrals are

\[
 J_r=\int\frac{d^Dk}{(2\pi)^D}\frac{(k^2)^r}{(k^2+\Delta)^2}
 =\frac{\Gamma(r+D/2)\Gamma(2-r-D/2)}
        {(4\pi)^{D/2}\Gamma(D/2)}\Delta^{D/2+r-2}.
\]

Their leading Laurent terms in D=4-epsilon, with common
1/(16 pi^2 epsilon) removed, are J0=2, J1=-4Delta, J2=6Delta^2. Angular
averages give delta/D and the three pairings/[D(D+2)]. Integrating x exactly,
the stress bubble is

\[
 B_{\mu\nu|\rho\sigma}
 =\frac{1}{16\pi^2\epsilon}
 \left[\left(m^4+\frac{m^2Q^2}{3}\right)
 (\delta_{\mu\rho}\delta_{\nu\sigma}
 +\delta_{\mu\sigma}\delta_{\nu\rho}
 -\delta_{\mu\nu}\delta_{\rho\sigma})
 -\frac{m^2}{3}J_{\mu\nu|\rho\sigma}(Q)\right]+cdots,
\]

where J has the four terms delta_mu,rho Q_nu Q_sigma and their symmetric
partners, minus 2 delta_rho,sigma Q_mu Q_nu and
2 delta_mu,nu Q_rho Q_sigma, plus delta_mu,nu delta_rho,sigma Q^2.
The ellipsis includes higher momentum terms and finite regulator terms, not
an exact finite-momentum completion.

Contract with Pi=I_sym-delta tensor delta/2. For conserved external sources,
Pi B has transverse coefficients 2m^4 and 2m^2Q^2/3; its remaining longitudinal
terms vanish in the source contraction. The executable script builds the
insertion polynomial, integrates its monomials, performs the parameter
integral, then contracts it. It does not insert the final kernel as its input.

Analytic continuation Q^2=-q^2 and the corrected loop sign yield

\[
 K=-\lambda\left(1-\frac{Dq^2}{12m^2}\right)I_{\rm sym}+O(q^4),
 \qquad
 \lambda=\frac{N(D/2-1)\Gamma(-D/2)m^D}
                {2V_0(4\pi)^{D/2}}.
\]

The leading four-dimensional residue is independently checked. The general-D
closed form is a comparison to the corrected published calculation, not an
independent evaluation of all its finite parts.

### Composite pole, physical helicities and matter coupling

Bubble resummation gives A=A0+KA, with A0=-Pi/(2V0). Its conserved-source
denominator is

\[
 1+\lambda-\lambda\frac{Dq^2}{12m^2}+O(q^4).
\]

At the tuned value lambda=-1,

\[
 V_0=-\frac{Nm^4}{32\pi^2\epsilon}+O(\epsilon^0),
 \quad M_P^2=-\frac{Nm^2}{48\pi^2\epsilon}+O(\epsilon^0),
 \quad A_{\rm pole}=-\frac{I_{\rm sym}-\eta\otimes\eta/2}{M_P^2q^2}.
\]

These are leading Laurent statements. Positive M_P^2 in this prescription
requires epsilon<0. Finite epsilon is a proxy for unspecified covariant UV
physics. Taking epsilon to zero does not leave a finite predictive Planck
scale without further input. The correct critical sign is lambda=-1; the
uncorrected original paper had the other sign.

For null q=(1,0,0,1), a conserved symmetric source has six independent
components. Its polarization residue obeys exactly

\[
 T_{\mu\nu}T^{\mu\nu}-\frac12T^2
 =\frac12(T_{11}-T_{22})^2+2T_{12}^2.
\]

This positive quadratic form has rank two. It identifies the two radiative
helicity components in the conditional massless exchange. It does not prove
positivity of the whole interacting Hilbert space or absence of every other
pole. A timelike pure spatial-trace source instead has residue -3/2 in our
normalization. Promoting this numerator unchanged into a massive propagator
would require a separate analysis; it cannot be used as proof of a ghost at
the exact massless point.

For multiple scalar masses, the total stress square fixes the mixed terms.
Let S4=sum_i Ni mi^4 and S2=sum_i Ni mi^2. Then the leading tuned quantities
are V0=-S4/(32 pi^2 epsilon), M_P^2=-S2/(48 pi^2 epsilon), and the same pole
couples to each species. The supplied single determinant/shared constraint
does the selection work. This is not a derivation of the observed fermion,
gauge-field and scalar inventory.

In a hypothetical valid Lorentzian interpretation, static source vertices
E00=2M^2 and the ordinary Born normalization give

\[
 \widetilde V(\mathbf k)=-\frac{M_1M_2}{2M_P^2\mathbf k^2},
 \quad V(r)=-\frac{M_1M_2}{8\pi M_P^2r},
 \quad G_N=\frac{1}{8\pi M_P^2}.
\]

**DERIVED conditionally:** the Newtonian exchange law. A Lorentzian theory and
the requisite asymptotic states have not been earned by the Gram action.
The published cubic-graviton comparison is useful supporting literature, but
we have not independently recomputed its six-scalar amplitude or proved full
nonlinear Einstein equations. The heat-kernel R/6 coefficient gives an
independent normalization control for M_P^2 under the usual auxiliary-field
loop prescription; it does not establish quantum equivalence of measures.

### The constraints and subsystem question in the uploaded note

The note correctly asks for more than a graviton-shaped pole. In this model
bare phi^a(x) transforms under coordinate gauge transformations and is not a
physical local observable. A tensor product over coordinate sites therefore
cannot simply be asserted to describe independently accessible physical
subsystems. The constraint reduction and a relational observable construction
remain necessary.

Two published results make the obligations sharper. Marolf's
[kinematic nonlocality theorem](https://arxiv.org/pdf/1409.2509) requires a
Hamiltonian which is a gravitational boundary flux, matching time evolution
and a compatible bulk/boundary locality map; with local gauge-invariant
kinematics those assumptions freeze interior evolution. Our calculation has
not generated that nonlinear boundary Hamiltonian or its map to microscopic
observables. A linear pole alone does not establish the theorem's premises.

Donnelly and Giddings'
[dressing theorem](https://arxiv.org/pdf/1607.01025) concerns perturbatively
diffeomorphism-invariant operators around flat gravity, including their
dependence on asymptotic geometry when they carry Poincare charges. This
packet does not construct those dressings from the scalar theory. It cannot
claim to have solved physical subsystem identity, independent local algebras
or observer access merely by exhibiting T_full=0.

Hojman-Kuchar-Teitelboim is a conditional recovery of geometrodynamics with
substantial canonical/geometric premises. Our worldvolume bracket makes the
boundary particularly concrete: deformation consistency by itself does not
select Einstein gravity. Here its spatial metric is a positive pullback
metric, rather than an independent ADM metric with derived gravitational
kinetic and curvature terms.

Weinberg-Witten also cannot be dismissed just by renaming T. The full
coordinate-invariant action has vanishing constrained stress, whereas the
gauge/background-organized free matter stress sources the conditional pole.
The distinction avoids directly inserting that free tensor into the theorem
as the total physical covariant energy tensor. It does not replace a
construction of charges, physical states and a unitary quantum theory.

## 3. THE RESTRICTION

What survives is conditional massless exchange with two radiative helicities
and a common scalar coupling. The masslessness condition is supplied as a
critical/vacuum condition, and reparametrization symmetry is supplied by the
action. They are not independently predicted facts about nature.

There is **no demonstrated extra restriction beyond the complete conventional
gravitational baseline**. In particular, the leading apparent relation
M_P^2=-S2/(48 pi^2 epsilon) is not an independently testable prediction when
epsilon represents unspecified UV physics and is adjusted to the observed
gravitational strength. A physical UV specification and matching prescription
would be additional premises. They are absent here.

The exact exclusions we do establish are mathematical: the frozen real Gram
branch cannot be Lorentzian, and its classical embedding constraints do not
by themselves constitute a derivation of independent Einstein metric dynamics.
No new quantum universal theorem about all emergent gravity follows.

## 4. THE COUNTERMODEL

The strongest conventional comparator is ordinary Einstein gravity plus the
same massive scalar inventory and all allowed covariant EFT operators:

\[
 S_{\rm std}=\int\sqrt{|g|}\left[
 \frac{M_P^2}{2}R-\Lambda+
 \sum_a\left(\frac12g^{\mu\nu}\partial_\mu\phi^a\partial_\nu\phi^a
       -\frac12m_a^2(\phi^a)^2\right)+\sum_\alpha c_\alpha\mathcal O_\alpha
 \right]d^4x,
\]

with curvature-sign convention chosen to give positive graviton residue.
This theory has the same two-helicity pole, universal scalar vertex and Newton
potential. Matching its low-energy parameters to the reconstructed pole
reproduces every claimed soft observable here. The corresponding matter loop
is also available in conventional EFT. The computed soft image is therefore
not separated from that baseline by a new allowed region.

Curvature-squared and higher operators are not prohibited by coordinate
invariance. Their finite coefficients are not fixed merely by calculating a
divergence. If the microscopic determinant model and its UV measure were fully
specified they could potentially predict matching data; those predictions are
not obtained here. Omitting these coefficients from the comparator would
manufacture novelty.

An unequal-species-coupling model would violate the shared-total-stress premise
that generated our relation. It is not offered as a same-premise comparator.
Detuning the flat-space bubble denominator is a sensitivity diagnostic, not a
certified consistent alternative vacuum. This matters because the original
clock-field family develops subtleties, including a singular static-gauge
limit and massive/non-Pauli-Fierz ghost concerns, discussed by its author.
We did not select that family as a repair.

## 5. THE TEST

Run from the repository root:

```bash
python research/constructive_physics/attempt03/check.py
```

Output:

```text
CP03_CONDITIONAL_COMPOSITE_POLE_CHECKS_PASS
groups=10
NOT_YET_ESTABLISHED_AS_PHYSICAL_GRAVITY_OR_NEW_GRUT_LAW
EXCLUDED_FOR_FROZEN_POSITIVE_GRAM_ACTION
```

[results.json](results.json) stores the exact rational tensor cells, source
residue matrix, fixed classical fixtures, tuning values and floating errors.
The independent bubble calculation checks 108 transverse tensor cells.
Three Hessian finite-difference controls have maximum absolute error
1.483e-10. The resummation's largest relative floating discrepancy is
1.104e-13. Neither number is an apparatus error bound.

At the diagnostic N=200, m=1, epsilon=-0.01, the leading Laurent expressions
give V0=63.3257397765 and M_P^2=42.2171598510 in the stated software units.
The full general-D expressions differ from those leading values by about
1.22% and 1.46%. This is explicitly recorded; finite regulator corrections are
not treated as machine error or zero physics.

There were two script executions: the initial nine groups passed; a final
explicit Newton/Born normalization group implementing the frozen
correspondence obligation was then added and all ten groups passed. No
candidate, parameter range or physical action was changed. No failed
scientific fixture was repaired into a survivor. Fixed gradient recipes and
derivative step are routine implementation choices disclosed in code; they
were not individually listed in the action freeze.

**Physical falsification boundary:** a genuinely identified scalar sector,
independently known inventory/masses and physically specified UV completion
could be compared against gravitational and scattering observations. We have
no such physical calibration or unique extra observable forecast. Newtonian
gravity or helicity measurements alone would not distinguish this microscopic
origin from conventional GR. The present calculation is consequently a
mathematical construction/control, not an experiment proposal certified for
testing new GRUT physics.

## 6. THE VERDICT

**NOT YET ESTABLISHED.** Retain the conditional composite-pole derivation and
the exact signature/constraint results. Reject promotion to a working physical
theory of reality. Specifically:

- Real Lorentzian classical emergence from the frozen positive-Gram action is
  excluded. That branch is not rescued by changing the target signature.
- The formal large-N spin-2 pole, its helicity residue and scalar coupling are
  reproduced; full quantum Lorentzian reconstruction remains unproved.
- A worldvolume constraint algebra is derived, but the quantum physical
  observable algebra, gravitational boundary Hamiltonian and full nonlinear
  Einstein dynamics are not generated by this packet.
- Critical/vacuum conditions and regulator physics remain supplied. No new
  joint-observable exclusion beyond complete conventional GR/EFT is obtained.

Evidence grade: author analytic proofs for the stated classical statements;
independent author rederivation of the leading loop pole and exact tensor
algebra; numerical corroboration of specified implementation controls. No
external review, experimental validation, complete UV theory or independent
quantum constraint proof.

The outstanding obligations are specific mathematical/physical constructions,
not additional routine software fixtures. No alternative law is formulated to
fill them in this attempt. CP-03 is separately accounted; CP-04 is not opened.
Stage 3 stays inactive and Card 3 is untouched. PR #12 remains draft.
All 117 unique files referenced by the existing scientific/routing manifests,
plus the three additionally frozen CP-01 program files, are verified unchanged.
The broader GRUT mission remains active; this result does not establish its
desired theory.
