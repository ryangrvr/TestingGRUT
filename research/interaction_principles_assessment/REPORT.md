# Physical origin of interactions: established derivations and one deeper constraint

2026-10-10. Starting development revision:
`45bfafbacc86994df5f34c96b4e25b40fff0f265`.

**Verdict: conditional interaction-selection mechanisms exist; an ultimate origin
or universal restriction beyond established physics is NOT YET ESTABLISHED here.**
The concrete deeper lead is a published axion/gravity positivity argument, with a
precisely identified physical assumption. This assessment responds to the owner's
new request. It introduces no GRUT generating law, card score, law search or new
foundations queue. Approach B stays frozen; Card 2 stays unopened.

## 1. What established physics actually derives

Weinberg's soft-particle argument assumes Lorentz-invariant scattering, its pole
structure and existing massless particles of the specified spin [S1]. For hard
future-directed on-shell momenta p_i, incoming/outgoing signs eta_i, soft outgoing
null momentum q and polarization epsilon, the leading emission factor has form

\[
\mathcal M_{n+1}\sim\mathcal M_n\sum_i\eta_i g_i
\frac{p_i^{\mu_1}\cdots p_i^{\mu_s}
\epsilon_{\mu_1\cdots\mu_s}}{p_i\cdot q}.
\]

Decoupling a gauge polarization gives

\[
s=1:\ \sum_i\eta_i g_i=0,\qquad
s=2:\ \sum_i\eta_i g_i p_i^\mu=0.
\]

The first is charge conservation. The second, over sufficiently rich scattering
connecting species, forces a common gravitational coupling by momentum
conservation. Neither establishes the existence of those particles or fixes the
common strength. These are conditional standard-physics restrictions.

Benincasa/Cachazo's four-particle test additionally assumes constructibility:
appropriate complex-momentum deformations reconstruct the amplitude without an
unknown boundary term. In its stated massless four-dimensional class it forces
spin-1 cubic couplings to obey a Lie-algebra Jacobi identity [S2]. Passing one
four-particle test is not proof of an all-orders theory. Higher-derivative,
nonconstructible interactions are outside that uniqueness argument.

Our exact controls expose both selection and surviving freedom. For A+B -> A+B,
the disclosed on-shell momenta give a rank-one two-species gravitational Ward
matrix: its null direction is (1,1). All common coupling values satisfy it, and
independently conserved electric charges remain free. Scaled SU(2) structure
constants satisfy Jacobi, whereas a disclosed antisymmetric tensor fails it.
See [generated coverage](CONTROL_TABLES.md). These finite checks illustrate the
equations, not a proof of an unrestricted S-matrix or the origin of all forces.

## 2. A proposed deeper connection, with its premise boundary

The electric Weak Gravity Conjecture compares charged states with gravitational
extremality. In four-dimensional Einstein-Maxwell theory without long-range
scalars, its familiar threshold is m <= sqrt(2)e|q| Mbar, where
Mbar=(8 pi G_N)^(-1/2), the Maxwell kinetic term is -F^2/4 and physical charge is
e q [S3,S4]. The claim is existential over the relevant full spectrum. A missing
light state is not a full-spectrum refutation; low-energy EFT consistency is not
a certified quantum-gravity completion.

Recent restricted results do not remove that distinction. A 2026 bosonic-string
account derives a sublattice result inside its specified string construction [S5].
Gravitational positivity arguments also have limits: the massless graviton pole
obstructs ordinary forward-limit reasoning, and a dispersive analysis allows
small black-hole-WGC violations within its assumptions [S7,S8]. Neither an
EFT example nor such a bound certifies a full quantum-gravity countermodel.

The most concrete connection examined is Di Ubaldo et al.'s 2026 argument [S6]:
gravitational path-integral norm positivity, together with included stable axion
wormholes and controlled higher connected contributions, requires sufficiently
strong shift-breaking instanton effects. For their unit-charge, single-axion
AdS_{d+1} case,

\[
S_{\rm inst}\le {1\over2f_a}
\sqrt{\frac{\pi(d-1)}{8dG_N}}.                 \tag{1}
\]

The paper explicitly leaves higher-boundary suppression and actual restoration
of all positivity conditions unresolved. It also discusses exclusion of the
wormhole saddles as an alternative. This is a conditional published argument,
not a verified universal theorem or a GRUT discovery.

## 3. Exact independent control of the decisive moment step

Here is the mathematical implication we checked, separately from any gravity
calculation. Let X be a nonnegative random variable with finite first three
moments m_j=E[X^j]. Set mu=m_1>0, v=m_2-mu^2 and
kappa_3=m_3-3mu m_2+2mu^3. For any real u,w,

\[
E[X(u+wX)^2]\ge0.
\]

Thus the shifted moment matrix is positive semidefinite:

\[
H_1=\begin{pmatrix}m_1&m_2\\m_2&m_3\end{pmatrix}\succeq0,
\quad
\det H_1=\mu\kappa_3+\mu^2v-v^2\ge0.
\]

Consequently,

\[
\boxed{\kappa_3\ge\frac{v^2}{\mu}-\mu v.}       \tag{2}
\]

This is ordinary moment positivity. To turn it into a variance prohibition,
extra control of kappa_3 is necessary. For example, if kappa_3 <= C v for a fixed
finite C>=0, then either v=0 or

\[
v\le\mu^2+\mu C.                             \tag{3}
\]

More generally, along an approach with bounded mu and v->infinity, a uniformly
valid kappa_3=o(v^2) conflicts with (2). Divide the determinant inequality by v^2:
the left term mu kappa_3/v^2 tends to zero, whereas 1-mu^2/v tends to one.
Suppression at each fixed parameter away from a singular boundary does not
establish the needed joint-limit control. That uniformity must come from the
physical saddle expansion, not the finite algebra fixtures.

### Hostile positive control

For every mu>0 and every finite v>=0, choose

\[
b=\mu+v/\mu,\qquad
P(X=b)=\mu/b,\qquad P(X=0)=1-\mu/b.
\]

Directly, E[X]=mu, Var(X)=v and
kappa_3=v^2/mu-mu v. Equality holds in (2) for every v. All moment matrices are
positive: each is a positive weighted sum of outer products. A strictly positive
variant sets a lower value a in (0,mu), upper value b=mu+v/(mu-a), and upper
probability (mu-a)/(b-a); it also permits arbitrary v at the same mean.

These controls defeat an inference from positivity alone to bounded variance.
They **do not** refute the conditional gravitational argument: they deliberately
allow the large higher moments that its suppression assumption excludes.
No wormhole geometry, quantum-gravity completion or accessible apparatus is
constructed by these distributions. The numerical examples and exact principal
minors are emitted into `results.json` and `CONTROL_TABLES.md`.

The code also checks a determinant-free geometric flux proxy. With r=exp(-C)
and t=exp(x), its two tails converge exactly when r<t<1/r, reproducing the
strip |x|<C for sum_m exp(-C|m|+mx). This is not an evaluation of a gravitational
one-loop determinant or a certificate for the full path integral.

## 4. What must be supplied, and what would count as a distinction

The established derivations supply spacetime, quantum scattering, particle
content, helicity and analytic/factorization assumptions. The gravitational
proposal additionally supplies an Einstein/axion EFT, periodicity and f_a,
the path-integral gluing prescription and admissible saddles. They do not derive
GRUT's persistent identity/interface pipeline, a nonzero force, the observed
particle inventory or all interaction strengths.

Our inference: positive Hilbert-space norms are already a requirement of quantum
physics. A new consequence can arise when a physical prescription for gravity
connects them to interaction parameters. Its novelty relative to the baseline
must reside in that prescription or its justified consequence, not in renaming
positivity. It cannot earn a GRUT exclusion by opposing a baseline that was
allowed unphysical negative norms in the first place.

A complete comparison would need a single conventional quantum-gravity model
with the same operational data and resources, plus independently justified
control of the relevant higher contributions. Neither is established here.
Nor is (1) currently an apparatus forecast: no measured complete instanton
inventory or independently certified actions/prefactors is supplied. A potential
harmonic's amplitude alone does not separate its instanton action from its
prefactor or rule out other contributing sectors. No experiment is specified as
feasible merely from the inequality.

If these ideas were converted into a novel GRUT generating law, supplied physical
structures would require actual frozen L_0 pricing, authorization and the draft
ledger. No such conversion or pricing is done in this assessment.

## 5. Reproduction and verdict

Run `bash research/interaction_principles_assessment/reproduce.sh` from the
repository root. Python's standard library suffices. Ward, Jacobi, moment and
geometric controls use exact rational arithmetic. The two displayed coefficient
evaluations use floating point and are not interval certificates. The suite
does not compute a physical path integral, reproduce the complete string proof,
fit measurements or establish an experimental result.

**Established:** conditional consistency can force interaction form and coupling
relations. **Concrete deeper lead:** global quantum-gravity consistency may force
nonperturbative interactions, with a published quantitative argument.
**Not established:** the required physical uniformity, a universal origin of
interactions, or a new GRUT restriction beyond the full conventional baseline.
The mathematical controls are author-produced evidence pending independent
review. Card 1 remains author-killed pending review; Card 2 is unopened and two
slots remain. Existing frozen packets and the card ledger are untouched.

Primary sources and exact read scopes are in `SOURCES.json`; the search is bounded,
not a proof that no stronger result or principle exists.
