# Formation prerequisites from the unchanged scalar potential

The new calculations use the classical zero-portal scalar sector; the
specified tiny portal supplies the previously computed production and
decay inputs. Its backreaction on these stationary coefficients is
negligible at the order being kept. This is not a quantum or thermal
completion. Natural units throughout.

## 1. Population and energy accounting

For X=f exp(-i omega t)/sqrt(2), D=0,

\[
U=\frac\lambda4(\phi^2-v^2)^2+
\frac{g^2}2\phi^2 f^2+\frac{\lambda_X}4 f^4,
\quad n_Q=\omega f^2.
\]

Take g=.1, lambda=.01, lambda_X=.001, v=1 GeV and the previously
solved E_Q=3639.242286 GeV at Q=39944. With the supplied final
B_DM=rho_DM/s_vis=4.4e-10 GeV, an all-object final population needs

\[
Y_{obj}=B_{DM}/E_Q,\qquad
Y_Q=QB_{DM}/E_Q.
\]

At T_vis=.020 GeV, s_vis=2pi^2*10.75*T_vis^3/45. Its n_Q is an
input calibrated to the desired end state, not a computed genesis.
No charge-violating interaction has been introduced. A global net
charge cannot arise from zero through the specified portal. Opposite
charges can be pair-produced; this observation does not forbid neutral
ensembles of objects and antiobjects.

For an initially empty +v vacuum and the parent's cold-source history,
all produced hidden energy is already bounded by

\[
B_{new}\leq\int\frac{\mathcal Q_{portal}}{s_{vis}H}dN.
\]

Even retaining all of it as stable matter gives at most 3.07e-8 of
B_DM in the reference enlarged response. No rearrangement of this
same energy can yield all of B_DM. A pre-existing false-vacuum
reservoir, charge population or extra producer changes the initial
state/problem and is not included in this empty-state conclusion.

If supplied free carriers instead assemble into these objects,

\[
B_{release}=Y_Qm_X-B_{DM}
=B_{DM}\left(\frac{Qm_X}{E_Q}-1\right),\quad m_X=gv.
\]

This is positive and must be conserved in the complete evolution.
For the dilute stationary condensate below, its energy per charge
differs from m_X only at relative order 1e-12, so using the free-carrier
release is adequate here. If fraction f_rest of B_release remains as
s rest mass after emission/redshifting, Y_s=f_rest B_release/m_s.
On the parent's radiation-era decay-temperature estimate,
rho_s/rho_rad approximately equals 4 f_rest B_release/(3T_decay).
This is a conditional energy diagnostic. An unknown formation spectrum
or branching ratio has not been fixed; neither a new nuclear limit
nor an actual injection fraction is derived from it.

## 2. Perturbative two-body channel and gas degeneracy

The original sXX* vertex is 2g^2v. Dividing the single-exchange
amplitude by the nonrelativistic external-state factor 4m_X^2 and
Fourier transforming gives

\[
V_{XX}(r)=-\alpha\frac{e^{-m_sr}}r,\qquad
\alpha=\frac{g^4v^2}{4\pi m_X^2}=\frac{g^2}{4\pi},
\quad m_s=\sqrt{2\lambda}\,v.
\]

The positive lambda_X contact interaction is repulsive. For the
attractive comparison Hamiltonian -d^2/(2mu dr^2)+V, the l=0
zero-energy Birman-Schwinger kernel is

\[
K(r,r')=2\mu\sqrt{|V(r)|}\min(r,r')\sqrt{|V(r')|}.
\]

It is positive with trace

\[
I=2\mu\int_0^\infty r|V(r)|dr
=\frac{2\mu\alpha}{m_s}=5.62698\times10^{-4},\quad\mu=m_X/2.
\]

Its largest eigenvalue is no larger than its trace. I<1 therefore
precludes the unit eigenvalue needed for a bound state, and negative
energies have still smaller free Green functions. Higher partial waves
add centrifugal repulsion; an additional positive contact interaction
cannot create binding. This proves absence for this leading
nonrelativistic comparison potential, not for all relativistic,
multibody or coherent-field configurations. The mediator is not
light compared with m_X, so a complete relativistic few-body analysis
is not silently claimed. Weak exchange and the large margin make
this a useful perturbative seed diagnostic.

For an ideal nonrelativistic Bose gas dominated by the supplied charge
sign, with negligible antiparticle density, the degeneracy threshold is

\[
T_{BEC}=\frac{2\pi}{m_X}\left[\frac{n_Q}{\zeta(3/2)}\right]^{2/3}.
\]

At the supplied density it is 1.06454e-7 GeV. This is a hidden gas
temperature, not a constraint on the visible bath. In adiabatic
decoupled nonrelativistic expansion n_Q proportional to a_scale^-3
and T_X proportional to a_scale^-2, so n_Q/T_X^(3/2) stays constant.
The attractive interacting gas may have more complex phase behavior;
this ideal-gas estimate is not a demonstrated condensation process.

## 3. Supplied rotating condensate and exact linear dispersion

A homogeneous stationary charged background on the positive-phi
branch has

\[
\phi_0^2=v^2-\frac{g^2}{\lambda}f_0^2,
\quad \omega^2=m_X^2-a f_0^2,
\quad a=\frac{g^4}{\lambda}-\lambda_X>0,
\quad n_Q=\omega f_0^2.
\]

Eliminating phi locally gives
U_eff(f)=m_X^2 f^2/2-a f^4/4 on this small-amplitude branch. Do not
extend its negative quartic past phi_0^2=0: the original two-field
potential remains bounded and must be retained during collapse.

Let eta=delta f, u=delta phi, and psi=f_0 delta(theta), for
X=(f_0+eta) exp[i(omega t+delta(theta))]/sqrt(2). The sign of the
conserved charge is immaterial here. For a Fourier momentum k,

\[
\ddot\eta+(k^2+A)\eta-2\omega\dot\psi+Cu=0,
\]
\[
\ddot\psi+k^2\psi+2\omega\dot\eta=0,
\qquad
\ddot u+(k^2+B)u+C\eta=0,
\]

where A=2lambda_X f_0^2, B=2lambda phi_0^2,
C=2g^2 phi_0 f_0. Hence the exact local growth determinant for
exp(Gamma t) is

\[
\left[(k^2+A+\Gamma^2)(k^2+\Gamma^2)
+4\omega^2\Gamma^2\right](k^2+B+\Gamma^2)
-C^2(k^2+\Gamma^2)=0.
\]

The code solves directly for positive Gamma^2 in scaled variables,
avoiding subtraction of nearly equal massive frequencies. The exact
positive band edge solves (k^2+A)(k^2+B)=C^2; its stable algebraic
expression is

\[
k_c^2=\frac{2(C^2-AB)}{\sqrt{(A-B)^2+4C^2}+A+B}.
\]

For k, Gamma much smaller than the mediator/massive-mode scales,

\[
\Gamma^2\simeq\frac{k^2(2af_0^2-k^2)}{4\omega^2},\quad
k_{peak}^2\simeq af_0^2,\quad
\Gamma_{peak}\simeq\frac{af_0^2}{2\omega}.
\]

An independent thermodynamic check uses the stationary chemical
potential mu_Q=omega:

\[
n_Q(\omega)=\frac{\omega(m_X^2-\omega^2)}a,
\quad c_s^2=\frac{n_Q}{\omega\,dn_Q/d\omega}
=-\frac{af_0^2}{2m_X^2-3af_0^2}<0.
\]

The full dispersion approaches -Gamma^2/k^2=-|c_s|^2 as k->0;
finite-k differences include the expected quantum-pressure term.
The reference mediator and nonrelativistic truncations are controlled
at relative order about 1e-12 during this onset. They are not
controlled through the eventual high-density collapse without solving
the full equations.

Use H=sqrt(pi^2 g_star/90)T_vis^2/M, M=2.435e18 GeV and
g_star=10.75. With Gamma_peak/H=4.5965e8, a seed needing S growth
e-folds has H Delta t approximately S H/Gamma_peak while S is modest.
The recorded shot-noise diagnostic chooses a coherent-state cell
with Q_cell=n_Q(2pi/k_peak)^3 and amplitude noise of order
Q_cell^-1/2. It gives S=.5 ln Q_cell and H Delta t=2.59e-8.
This is an illustrative seed estimate, not a derived primordial
fluctuation spectrum or a stochastic field simulation.

The wavelength charge Q_cell is a length-scale diagnostic only.
Even if it were a fragment charge, it is approximately 2.15e10,
whereas the target is 39944. The nonlinear efficiency, losses,
mergers and final distribution are necessary missing outputs.

For context, minimizing 2U/f^2 gives the bulk phi=0 phase with
f_bulk^2=v^2 sqrt(lambda/lambda_X) and
omega_floor=(lambda lambda_X)^(1/4)v. This is consistent with the
older charged-phase analysis, not an independently new mechanism or
a selection of finite Q.

## 4. Neutral wall tension and observational scope

At X=D=0 the neutral static wall is

\[
\phi(z)=v\tanh(z/\delta),\quad
\delta=\frac{\sqrt2}{\sqrt\lambda v}.
\]

Its gradient and potential energies are equal point by point. Therefore

\[
\sigma=\int dz\left[\frac12(\partial_z\phi)^2+
\frac\lambda4(\phi^2-v^2)^2\right]
=\frac{2\sqrt{2\lambda}}3v^3.
\]

The X and D potentials are nonnegative for neutral static variations;
zero fields minimize them for this kink. Occupied charged walls and
their network dynamics are a different question. The small even
portal leaves the two sign vacua degenerate. Adiabatically eliminating
the Higgs gives lambda_eff=lambda-kappa^2/(4lambda_H), a relative
shift of about 7.5e-17 at the reference kappa with
lambda_H=m_h^2/(2v_H^2); Higgs-gradient effects are also higher order.
Neither effect provides an odd bias term.

The comparison uses Caloni et al. 2602.20050v1, whose conservative
stable-network analysis gives sigma^(1/3)<0.81 MeV at 95% credibility.
It is expressed in tension, not merely the name of a transition scale.
For our fixed lambda it corresponds to v<1.77969 MeV under those
network assumptions. This rejects the specified ordinary stable
network outcome at GeV scales. It does not establish formation,
percolation, memory erasure, charged-wall dynamics or survival, or
apply the same likelihood to an early-decaying network. We do not
compute a new CMB significance from a tension ratio.
