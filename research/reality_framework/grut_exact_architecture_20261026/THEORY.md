# GRUT: an explicit architecture for a closed physical world

26 October 2026, America/Chicago. Codex / AInstein.
This is a specified working architecture and two new mathematical constructions.
It is not a claim to a completed fundamental theory or an observed new sector.

The most defensible concrete proposal is **a constrained quantum description
of one universe, whose low-energy dynamics are a covariant field theory and
whose observations are correlations recorded by physical systems inside it**.
The dark scalar construction is one candidate matter sector in that theory.

The complete mathematical datum is

\[
\mathfrak T=(\mathcal K,\{\widehat C_A\},S_{\rm EFT},
\vartheta,\mathcal B,\{E_{a|T=t}\}).
\]

Here K is a kinematic state space, C_A are the constraints, S_EFT and its
matching data vartheta specify dynamics, B specifies a physical state or
preparation measure, and E specifies physical observation effects.
An actual world has a state; the space of possible states can still have
dimension greater than one. Nothing here identifies those two statements.
In a complete model these entries must be compatible: the constraints and
physical inner product must come from a consistent quantization of the
action. The finite construction below illustrates the dictionary; it has
not been derived from the gravitational field action.

## 1. Physical content and spacetime

Use ordinary Standard Model matter, one metric g, a real binding scalar phi,
three complex scalars X,Y,Sigma, and a hidden gauge connection A. Their hidden
charges are 1,3,6. In a definite compact branch choose spatial topology T^3.
This topology is a supplied hypothesis, not a consequence of internal baths.
Local physics is studied on a globally hyperbolic spacetime with Cauchy slices.

The organizing condition is that clocks, controls, apparatus and reservoirs
are all physical subsystems. It does not mean that every reservoir is thermal
or every reduced evolution is Markovian. In gauge/gravitational theories
subsystem algebras and dressings require care; an exact tensor factorization
of arbitrary physical spatial regions is not being assumed.

The old reporter D is not required as a new fundamental degree of freedom
in this core architecture. Ordinary matter can carry records. Earlier D-probe
calculations remain conditional probe results, not an automatically integrated
observer model. This choice changes no pinned D=0 classical benchmark and
claims no calculation of D quantum effects.

## 2. One local action

In signature (−+++), use

\[
\begin{aligned}
S_{\rm EFT}=\int_Md^4x\sqrt{-g}\bigg[&\frac12F(\phi,z,H_{\rm SM})R
 +\mathcal L_{\rm SM}-\frac14F_{X\mu\nu}F_X^{\mu\nu}
 -\frac{\epsilon_{BX}}2B_{\mu\nu}F_X^{\mu\nu}\\
 &-\frac12(\nabla\phi)^2-\sum_{z=X,Y,\Sigma}|D z|^2
 -V_{\rm new}(\phi,z,H_{\rm SM})
 +\mathcal L_{\rm higher}\bigg]+S_{\rm boundary}.
\end{aligned}
\]

D_mu z=(nabla_mu−i e q_z A_mu)z. F contains the reduced Planck mass squared
and allowed nonminimal couplings. L_higher contains curvature counterterms
and higher-dimensional operators, with a stated cutoff and matching
coefficients. S_boundary supplies the usual variational terms for a chosen
Cauchy slab; compact space has no spatial boundary but can still have initial
and final slices. [ACTION.md](ACTION.md) gives the scalar operator basis.

At the classical matching reference the tested matter potential is

\[
\begin{aligned}
V_0={}&\rho_0+\frac{\lambda_\phi}{4}(\phi^2-v^2)^2
 +g_X^2\phi^2|X|^2+\lambda_X|X|^4\\
 &+M_Y^2|Y|^2+\lambda_Y|Y|^4
 +\lambda_\Sigma(|\Sigma|^2-h^2/2)^2\\
 &+(\kappa Y^*X^3+\mu\Sigma^*Y^2+\mathrm{h.c.}).
\end{aligned}
\]

Additional allowed coefficients cannot be declared zero at every scale.
They are matching conditions. The present light target and heavy ratios are
listed in the parent [INPUTS.json](../grut_heavy_vector_budget_20261025/INPUTS.json).
They are supplied numbers. Renormalized rho_0 is an independent gravitational
vacuum coefficient; its small observed value is not derived by this action.

A concrete covariant representative of the old inflationary −10 H^2|X|^2
term is the standard nonminimal coupling xi_X R|X|^2 in the potential, with
xi_X=−5/6. It gives −10 H^2 in ideal de Sitter, zero in ideal radiation and
−2.5 H^2 in an ideal matter era. Its removal follows the metric evolution
instead of an externally timed switch. This is an additional specified
coupling, not a selection by closure. A real inflaton/reheating history and
the phase holder are still required; no new cosmological trajectory is run.

## 3. Physical states and internal time

The intended canonical description imposes momentum, Hamiltonian and gauge
constraints. Schematically,

\[
\widehat{\mathcal H}_\perp(x)|\Psi\rangle=0,
\quad\widehat{\mathcal H}_i(x)|\Psi\rangle=0,
\quad\widehat{\mathcal G}_a(x)|\Psi\rangle=0.
\]

This does not state a measured additive total energy of zero. A consistent
quantum constraint algebra, physical inner product and regulator for the
full gravitational problem are not constructed here. Einstein gravity is
adopted in the EFT; microscopic geometry and gravity have not been derived.

Observations are conditional on an internal clock T. In a specified physical
description, probabilities have the form

\[
p(a|T=t)=\operatorname{Tr}_{\rm phys}(\rho_{\rm phys}E_{a|t}).
\]

The Born rule is supplied. The effects must be physical and normalized, not
just arbitrary clock projectors that fail the constraints. We now give a
fully executable finite realization of this dictionary in [RECORDS.md](RECORDS.md).
It does not complete the continuum gravitational dictionary.

## 4. An exact small world with a clock and two records

Take eight clock readings, one matter bit and two record bits. A single
time-independent constraint checks an internal sequence of gates: make a
superposition, copy the pointer bit twice, apply a phase and reverse all
operations. Every control is encoded in that constraint. There is no
external intervention during this modeled history.

For each supplied input state psi_0, the stationary state is

\[
|\Psi\rangle=\frac1{\sqrt8}\sum_{t=0}^{7}|t\rangle W_t|\psi_0\rangle,
\qquad C|\Psi\rangle=0.
\]

For the chosen blank-record input, the records share one bit of pointer
information at readings 3–5. Their local probabilities are one half, while
the whole retains coherence. At reading 7 the information has been erased
and the matter bit is pure again. The global state is stationary throughout;
the conditional records change. This is exact finite quantum mechanics with
a supplied circuit and input, not a spacetime-emergence result.

The same constraint permits every branch probability p between zero and one.
All those history states are pure; purity cannot choose among them. The
maximum-entropy physical state instead gives the identity conditional state
at every reading and no new pointer correlations. Closure plus maximum
ignorance therefore supplies no record-producing preparation in this model.

An extra state filter can choose the blank-record example. Its compatibility
is explicitly checked, but choosing the filter is an additional state law.
The Hilbert space remains eight-dimensional; selecting one actual state is
not a claim that only one possible state exists.

## 5. Persistent matter and cosmic history

The scalar mechanism makes local binding conditions dynamical: a concentration
of X changes phi, lowering its mass cost and allowing a bound configuration.
After Sigma breaks the gauge group, charge-six interactions can transfer the
approximate X charge while the total gauge charge, including Sigma, remains
zero on a compact slice. A residual Z6 does not alone prove the lifetime of
a large object; the specified leakage channels still have to be calculated.

The present charge-transfer example is conditional on its supplied phase,
preparation and reheating assumptions. Its abundance is being tested by
Claude's kinetic calculation. Neither this finite clock nor the new curvature
representative supplies the missing phase or formation kernel. The low-entropy
initial condition remains in B. No special final entropy or record count has
been inserted to fix the cosmological constant.

## 6. How the architecture becomes predictive

The required chain is fixed action/matching data and preparation → formation
and charge distribution → decay/survival history → final populations and
profiles → detector or cosmological distributions. The same parameters must
be used through the chain. The existing multi-target recoil-shape comparison
family is retained; a full experimental likelihood and held-out test remain
to be specified. Conditional internal checks do not establish an observed
new particle or an empirically distinguishing signal.

The specifically missing GRUT step is a rule that selects some dynamics or
states, or gives another operational consequence beyond the adopted physics.
This architecture says exactly where it must enter: S_EFT, B or the physical
observation rule. It cannot be replaced by the statement that the bath is
internal. The exact finite state-filter example shows what an added
restriction does mathematically, without pretending to have derived the
cosmological one.

Eight mathematical groups pass for the new constructions. See
[REPORT.md](REPORT.md), [HISTORY_RESULTS.json](HISTORY_RESULTS.json),
[OPERATOR_RESULTS.json](OPERATOR_RESULTS.json) and [REPRODUCE.md](REPRODUCE.md).
