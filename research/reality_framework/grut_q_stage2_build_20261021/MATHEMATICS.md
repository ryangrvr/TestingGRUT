# Stage 2 construction mathematics

21 October 2026. All computations are author calculations in the supplied
classical models; primary literature and read scopes are in SOURCES.md.
The four planned threads are completed here at their stated approximations.
The metastable version is an explicit additional potential/preparation
choice, kept separate from the symmetric baseline and the dissipative model.

## 1. A conserved vacuum-change scenario

Start with the fixed-mu q action and an explicit minimally coupled scalar s:

\[
S=\int\sqrt{-g}\,[K(q)R/2-\epsilon(q)
 -(\partial s)^2/2-V(s)]+S_d+S_r.
\]

Where K is linear and epsilon is quadratic about the selected branch,
eliminating q gives the exactly specified local representative

\[
S_g={M^2\over2}\int\sqrt{-g}[R+R^2/(6m_q^2)-2\Lambda],\qquad
F\equiv K/M^2=1+R/(3m_q^2),\quad c=M^2\Lambda.
\]

This is a supplied local choice of the free functions, not a derivation of
them from closure. The q equation is still
epsilon'(q)-K'(q)R/2=mu, with constant mu. Its selected value and c are
matching/state data. M is the reference gravitational normalization.
The underlying scalar-tensor/f(R) equivalence is due to KV and established
scalar-tensor physics; see the previous response packet for its derivation.

In numerical units M=1 let x=F-1 and v=Fdot. For
\(V(s)=\lambda_s(s^2-v_s^2)^2/4\), the closed homogeneous system is

\[
\begin{split}
\dot N&=H,&\dot x&=v,\\
\dot v&=-3Hv-m_q^2x+(\rho_d+4V-\dot s^2+4c)/3,\\
\dot H&=m_q^2x/2-2H^2,\\
\ddot s&=-3H\dot s-V_s,&
\rho_d&=\rho_{d0}e^{-3N},\quad\rho_r=\rho_{r0}e^{-4N}.
\end{split}
\]

The scalar trace is rho_s-3p_s=4V-sdot². Omitting its kinetic term would
give a different gravitational response during the transition.
The exact Friedmann constraint and propagation identity are

\[
C=3(1+x)H^2+3Hv-\rho_d-\rho_r-\dot s^2/2-V-c
 -3m_q^2x^2/4=0,\qquad \dot C=-4HC.
\]

The matter identity is dot rho_s+3H sdot²=0. Potential energy becomes
kinetic energy and contributes to expansion work. No external time-dependent
potential is removed without a receiver, and no globally conserved GR
total energy is assumed.

The symmetric scenario starts near the unstable maximum and rolls to
one minimum. It is a homogeneous symmetry-breaking roll with a supplied
seed, not a calculation of thermal nucleation or the actual SM crossover.
At c=0 its late-time oscillations redshift and the tested background
decelerates. At supplied c>0 it approaches an accelerating background.
The q potential alone need not vanish on the latter background:

\[
H^2=c/3,\qquad x=4c/(3m_q^2),\qquad
U_q-c=4c^2/(3m_q^2).
\]

That static deformation is part of the same de Sitter solution, not an
additional independent residual: K changes consistently in Friedmann's
equation, leaving H²=c/3 exactly.

### Separate formulation used for verification

Let beta=sqrt(2/3), F=exp(beta chi), g_E=F g_J and d tau_E=sqrt(F)dt_J.
The Einstein-frame scalar metric is diag(1,F^-1), with potential

\[
W(\chi,s)={c+3m_q^2(F-1)^2/4+V(s)\over F^2}.
\]

Writing u_chi=dchi/dtau_E and u_s=ds/dtau_E, the separate verifier uses

\[
\begin{split}
3H_E^2&=u_\chi^2/2+u_s^2/(2F)+W+\rho_d^E+\rho_r^E,\\
u_\chi'&=-3H_Eu_\chi-\beta u_s^2/(2F)-W_\chi+\beta\rho_d^E/2,\\
u_s'&=-3H_Eu_s+\beta u_\chi u_s-V_s/F,\\
N_J'&=H_E-\beta u_\chi/2,\quad t_J'=F^{-1/2},\\
\rho_d^E&=\rho_{d0}e^{-3N_J}/F^2,\quad
\rho_r^E=\rho_{r0}e^{-4N_J}/F^2.
\end{split}
\]

Primes in this paragraph denote Einstein proper time. In particular the
field-space kinetic connection and dust coupling are included. Translating
back gives H_J=sqrt(F)(H_E-beta u_chi/2). The verifier does not import
build.py. It compares the two solutions at the same saved Jordan times.

## 2. The gravitational-response bound includes derivative terms

In the analytic stable heavy branch with no independently excited large
q condensate,

\[
U_q-c\simeq{(\rho_d+4c)^2\over12M^2m_q^2}.
\]

This is not the full correction to expansion. Substitute
F=1+R/(3m_q²), R_0=(rho_d+4c)/M² and the GR leading solution into the exact
Friedmann equation. To first order in m_q^-2,

\[
3M^2H^2=\rho_d+\rho_r+c+
{R_0^2/12-H_0\dot R_0-R_0H_0^2\over m_q^2}\,M^2
+O(m_q^{-4}).
\]

Here H_0 means the leading GR solution at that time, not today's measured
Hubble rate. For conserved dust and radiation and constant c this becomes

\[
\delta\rho_{\rm bg}={1\over M^2m_q^2}
\left({3\rho_d^2\over4}+{2\rho_d\rho_r\over3}
-{4c\rho_r\over3}\right)+O(m_q^{-4}).
\]

In dust domination with c=0 the full source correction is **nine times**
the potential-only contribution. Its leading effective density scales as
a^-6; it does not behave as constant dark energy. A homogeneous massive
free q oscillation instead has a dust-like averaged energy, with its
amplitude supplied by preparation. Neither result sets a positive constant.

The same unscreened mode has relative Yukawa strength 1/3. Lee et al.'s
historical 2020 exclusion plot was viewed directly. Its alpha=1/3 crossing
is only visually bracketed at 45--60 micrometers, corresponding to
3.29--4.39 meV; this bracket is a reading estimate, not a new confidence
interval or likelihood analysis. The exact quoted alpha=1 endpoint cannot
be substituted for the alpha=1/3 crossing.

At a uniform 4 meV mass, the full dust correction using the listed present
reference densities is 4.21e-62 times rho_Lambda. Extrapolating the formula
to equality with rho_Lambda would require m_q about 0.57H_0, beyond the
heavy-field approximation. It diagnoses the need for a nonadiabatic light
mode, not a valid fitted solution. Screening requires the actual nonlinear
functions, source environments and perturbation calculation; it is not
included in this bound. These statements do not constrain an arbitrary
free integration constant c, an added vacuum phase gap, or every f(R) model.

## 3. Several matter vacua and conditional Higgs/top relations

Take q-independent matter minima with energies V_i and one convex q branch.
For a Minkowski solution in each phase with the same mu,

\[
\epsilon'(q_i)=\mu,\qquad\epsilon(q_i)+V_i-\mu q_i=0.
\]

Strict convexity makes q_i the same on that branch. Requiring every phase
to be Minkowski then implies V_i=V_j. But QV1 requires matching the eventual
vacuum; it does not require every other local matter minimum to be flat.
For epsilon=q²/2+1/2, mu=1, the phases V_1=0 and V_2=.01 give gravitating
energies 0 and .01. This elementary counterexample leaves QV1 intact while
refuting an automatic multiple-point requirement. If epsilon has distinct
branches or matter depends on q, even the same-branch implication needs
re-examination.

A useful positive possibility remains: if the final true vacuum is matched
to zero but the current phase is metastable, then its gravitating residual
is approximately Delta V=V_false-V_true. Its sign is positive. Its magnitude
requires the potential and state, not merely the closure condition.

For an additional SM-desert assumption, a high-field Higgs minimum h_* with
V_eff approximately lambda_H(h)h^4/4 is nearly degenerate with the EW phase
only when lambda_H(h_*) and beta_lambda_H(h_*) are nearly zero. This is the
published multiple-point criticality condition, not new GRUT physics.
Ignoring bottom and tau, the one-loop relation at lambda_H=0 is

\[
y_t(h_*)^4={2g_2(h_*)^4+[g_2(h_*)^2+g_Y(h_*)^2]^2\over16}.
\]

We solve this boundary and run the one-loop SM couplings downward for
three **supplied** h_* values. The verifier instead uses analytic gauge
solutions and fixed-step RK4 for the Yukawa/quartic subsystem. These are
scheme-specific running-coupling calculations; neither gives pole masses.
Precision thresholds, higher loops, gravity and possible higher Higgs
operators are not calculated. Near Planckian h_* the last issue is material.

If the high minimum's energy gap itself supplies the observed residual,
the required quartic-level mismatch is about 4rho_Lambda/h_*^4. At
h_*=2.435e18 GeV it is 2.87e-120. This is an exposed sensitivity, not
precision achieved by the one-loop computation. QV1 has not selected either
the critical scale or the top/Higgs parameters to that accuracy.

## 4. A closed model that forms mu dynamically

This is a separate KV-inspired phenomenological extension, not the
nondissipative action with a silently changed integration constant.
The original constant-mu equation is replaced by an effective source law;
a microscopic three-form current or quantum completion is not supplied.
Take K=M² constant and radiation as the receiver. In flat FRW,

\[
3M^2H^2=\rho_V+\rho_r,\quad
\dot\rho_V=-S,\quad\dot\rho_r+4H\rho_r=S,\quad
S={\gamma\over36}|H|R^2,\quad R=4\rho_V/M^2.
\]

The curvature-production ansatz and its leading graviton coefficient
gamma=1/(32pi²) are credited to KV 1601.00601 with that paper's adiabatic
state and correction-term qualifications. This calculation does not derive
the ansatz from a microscopic action.

One covariant continuum specification uses normalized radiation velocity u,
Theta=nabla_mu u^mu, T_V^{mu nu}=-rho_V g^{mu nu}, Einstein's equation,
the transport law u^mu nabla_mu rho_V=-(gamma/36)|Theta/3|R², and
nabla_mu T_r^{mu nu}=nabla^nu rho_V. The latter includes momentum transfer
when gradients exist. Setting the entire transfer vector to S u^nu while
allowing arbitrary vacuum gradients would be inconsistent. The homogeneous
system is solved here; full perturbative stability and quantum noise are not.

For illustration choose the published-type convex function
epsilon(q)=C+q²+q^-2 at positive q, in reduced Planck units. Then

\[
\mu=2q-2q^{-3},\quad r_V=C-q^2+3q^{-2},\quad
\epsilon''=2+6q^{-4}>0,\quad\dot\mu=S/q.
\]

Changing C changes the zero-energy target mu, while the source drives mu
toward that target on the positive expanding branch. This demonstrates
formation rather than simply imposing its final value. No retuning of mu
is needed between the three tested *constant-C* cases. A time-dependent
threshold must include its energy source and may leave this basin.

For r=rho/M^4, h=H/M, tau=tM and N=ln a, define b=4gamma/9. The exact result is

\[
{dr_V\over dN}=-br_V^2,\qquad
r_V(N)={r_{Vi}\over1+br_{Vi}N},\qquad
r_r(N)=r_{ri}e^{-4N}+b\int_0^N e^{-4(N-z)}r_V(z)^2dz.
\]

The receiver integral is checked independently by quadrature. There is no
positive fixed point in this branch: the attractor is zero. From r_Vi=1,
reaching the reference 7.17e-121 requires about 9.91e122 e-folds at the
graviton coefficient. A vacuum-dominated early-start asymptote gives

\[
r_V(\tau)\simeq[2\gamma\tau/(3\sqrt3)]^{-2/3}.
\]

At the reference age 13.8 Gyr it leaves approximately 2.24e35 GeV^4 under
those assumptions. These are finite-history failures of this natural-rate
example, not a no-go for every dissipative theory or arbitrarily old history.
The negative-vacuum control expands briefly and turns around; the same
positive source cannot lift a negative rho_V to zero.

## 5. Explicit positive metastable version QV1-M

Keep the field content of section 1, set c_q=0 by matching the eventual
true vacuum, and supply the new potential/preparation choice

\[
V_M(s)={\lambda_s\over4}(s^2-v_s^2)^2+
{\delta\over4}[2-3s/v_s+(s/v_s)^3],\quad\delta>0.
\]

Both s=+-v_s are stationary, V(-v_s)=delta and V(+v_s)=0. They are local
minima if delta<4lambda_s v_s^4/3. At tree level,

\[
m_-^2=2\lambda_s v_s^2-{3\delta\over2v_s^2},\quad
m_+^2=2\lambda_s v_s^2+{3\delta\over2v_s^2},\quad
\delta={v_s^2\over3}(m_+^2-m_-^2).
\]

The false phase has an exact late homogeneous solution
H²=delta/(3M²), F-1=4delta/(3M²m_q²). Finite trajectories for an exactly
false preparation and a small perturbation remain near that minimum and
accelerate; the preparation above the barrier instead rolls to the true
phase and decelerates in the tested window. This demonstrates a positive
residual with consistent gravitational response. It has not selected delta.

At delta=0 the potential has exact s->-s symmetry. If this is a symmetry
of all couplings and the renormalization procedure, the two phases remain
degenerate when it is unbroken. This supports technical protection of a
small symmetry-breaking gap; it does not compute its value, perform loop
matching or protect against unspecified symmetry-breaking UV operators.
The displayed mass/gap relation is a tree-level quartic relation, not an
all-order anomaly-protected identity.

The q matching still refers to the final fully renormalized vacuum. It is
a supplied state rule, and false-phase preparation is another input. No
quantum tunneling rate, false-vacuum population or lifetime is established.
The two new choices are explicit in INPUTS.json; QV1 and its failures are
preserved beside them.
