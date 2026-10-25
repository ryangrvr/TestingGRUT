# Full release dictionary

All unadorned masses and fields in the solver are in H_i units. Set
a=H_i eta, with a=1 at release and a'=1. For chi_i=a z_i/H_i,
the physical zero-velocity preparation is chi_i'=chi_i at a=1.

The conformal potential is

\[
U=\lambda_X|\chi_X|^4+\lambda_Y|\chi_Y|^4
+\lambda_\Sigma(|\chi_\Sigma|^2-a^2\hat h^2/2)^2
+a^2\hat M^2|\chi_Y|^2
+2\operatorname{Re}[\kappa\chi_Y^*\chi_X^3
+a\hat\mu\chi_\Sigma^*\chi_Y^2].
\]

The full equations are chi_i''=-partial U/partial chi_i*.
No mediator elimination is used. This formulation agrees with the
review's conformal equations; the independent control evolves the
physical complex fields in cosmic time with 3H friction instead.

## Conserved gauge charge and energy work

\[
N_i=-2\operatorname{Im}(\chi_i^*\chi_i'),\qquad
N_X+3N_Y+6N_\Sigma=0.
\]

The physical energy and conformal canonical Hamiltonian differ:

\[
a^4\rho/H_i^4=\sum_i|\chi_i'-\chi_i/a|^2+U,\qquad
\mathcal H_c=\sum_i|\chi_i'|^2+U.
\]

The diagnostic independently integrates

\[
\frac{d\mathcal H_c}{da}=\partial_a U
=2a\hat M^2|\chi_Y|^2
-2\lambda_\Sigma a\hat h^2(|\chi_\Sigma|^2-a^2\hat h^2/2)
+2\hat\mu\operatorname{Re}(\chi_\Sigma^*\chi_Y^2).
\]

The reported work residual is normalized to the initial canonical
energy. It is not a claim that expanding matter energy is constant.
Gauss error is normalized to the initial X circular-orbit charge scale,
so an exactly aligned zero-charge solution has a meaningful zero error.

## Charge, yield and the compulsory companion

For X alone, E_X=|chi_X'-chi_X/a|^2+lambda_X|chi_X|^4. Minimizing quartic
energy at fixed charge gives a circular orbit with

\[
N_{\rm circ}(E_X)=2\sqrt{2\lambda_X}
\left(\frac{E_X}{3\lambda_X}\right)^{3/4},\qquad
\epsilon_X=N_X/N_{\rm circ}(E_X).
\]

Window statistics include mean(epsilon_X) and mean(N_X)/mean(N_circ),
which are different operations. This fraction excludes Y, Sigma and
vector particles; it is therefore an optimistic leftover diagnostic.

N_eff=N_X+3N_Y=-6N_Sigma is also reported. The mu term changes N_eff;
late Y/anti-Y mixing can leave N_eff oscillating even when N_X has
settled. Their means and standard deviations must remain visible.

\[
T_R=\left(\frac{90M_{\rm red}^2H_i^2}{\pi^2g_*}\right)^{1/4},
\quad s_i=\frac{2\pi^2}{45}g_{*s}T_R^3,\quad
Y_X=N_XH_i^3/s_i.
\]

These formulas assume instant visible reheating and conserved visible
entropy. They do not provide dilution or a late charge-decay history.

In the quadratic Y approximation at fixed Sigma, with
D=M_Y^4-4|mu Sigma|^2>0,

\[
\alpha_6=\mu\kappa^2\Sigma^*/D,\qquad
\beta_6=-M_Y^2|\kappa|^2/D,\qquad
\frac{|\beta_6|}{2|\alpha_6|}=\frac{M_Y^2}{2|\mu\Sigma|}>1.
\]

At the late vacuum this ratio is 14.1421. An angular X^6 term without
its negative radial companion is a different model. Even small Y
quartic corrections do not by themselves justify freezing Sigma at
its vacuum value during the release.

## Late Y excitation diagnostic

Let A=a^2 Mhat^2, B=2a muhat chi_Sigma and J=kappa chi_X^3. Where
A>|B| and Y self-interactions are small, the instantaneous forced
quadratic solution satisfies A y_p+B y_p*=-J:

\[
y_p=(-AJ+BJ^*)/(A^2-|B|^2).
\]

Its derivative follows by differentiating that same linear system,
including A', B' and J'. For d=chi_Y-y_p, the positive quadratic
excitation energy is

\[
E_{Y,\rm exc}=|d'-d/a|^2+A|d|^2+
\operatorname{Re}(B^*d^2).
\]

The approximate comoving number J_Y=E_Y,exc/(a Mhat) neglects the
small Majorana mass splitting and the time variation of its eigenvectors.
It is a diagnostic, not a decay or thermalization simulation. If its
mean is conserved to Y decay, the mean Y charge vanishes, each Y produces
three X quanta and no number changing occurs afterward, the extra pairs
per absolute net X charge are approximately 3 J_Y/(2 |N_X|).
These are explicitly conditional assumptions. Excitation energy is not
added twice to the total scalar ledger.

## Gauge direction and vector source

For canonical scalar coordinates the gauge-orbit norm is

\[
F_g^2/H_i^2=2\sum_i q_i^2|\chi_i|^2/a^2,\qquad
m_A/H_i=eF_g/H_i.
\]

Its values and history supply coefficients for a finite-k coupled
perturbation problem. A homogeneous solution with zero Gauss charge
does not supply a gauge-invariant longitudinal occupation number.
In particular, substituting the late vector mass into a constant-mass
inflationary relic formula does not evolve this changing Higgs background.

VECTOR_INITIALIZATION.json also supplies a declared vacuum initialization.
For stationary inflationary amplitudes, a light vector and the assumed
Bunch-Davies state, the superhorizon Goldstone limit has
sqrt(P_G)/H_i=1/(2 pi), sqrt(P_alpha)=H_i/(2 pi F_g).
The corresponding Proca convention is
P_AL/H_i^2=(k/H_i)^2/[4 pi^2 (m_A/H_i)^2]. These are initial-mode
normalizations, conditional on that state and limit, not an evolved
occupation or gauge-independent stochastic angle. Their coupled transfer
through release has not been calculated.

## CP and spatial scope

Complex conjugation maps each solution to one with all N_i reversed.
For the real-coupling action and a CP-even preparation, the ensemble
mean vanishes exactly at every H_i. An imposed common nonzero delta
is a CP-odd initial condition. Its survival is a necessary test of a
possible supplied preparation, not a preparation mechanism.

Sensitivity of homogeneous histories is evidence about that transfer
map. A domain size, spatial charge power spectrum or volume-weighted
cosmic mean requires an inflationary fluctuation spectrum and finite-k
evolution; neither is inferred solely from a collection of k=0 runs.
