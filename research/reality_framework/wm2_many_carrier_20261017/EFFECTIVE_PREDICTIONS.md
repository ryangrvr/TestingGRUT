# Defined domain, units, and predictions that can be rejected

**Calculated predictions, not measured results.** All-size D1 review remains
pending. This document supplies a narrow empirical interface for the same
candidate, not a claim that the full law describes nature.

## 1. What effective physics has actually been obtained?

On a retained recorded interior D1 gives a one-dimensional, finite-range,
interacting hard-core quantum medium with two conserved species and internal
states. The dimensional class, quantum formalism, external time and relational
generator rule are inputs. The common-space error bound controls its finite-
time use within the moving-geometry law. A many-body continuum or hydrodynamic
limit has not been established.

In the one-carrier channel, with the other species absent, H_bond is exactly
zero. A fixed internal level 0 gives the exact free path Hamiltonian

    H_s=-j_s sum_(x>=0)(|x+1><x|+|x><x+1|).                   (E1)

For bulk plane waves on an infinite chain the physical dispersion is

    E_s(k)=-2E_0 j_s cos(k a),
    E_s(k)-E_min=E_0 j_s(ka)^2 + remainder,
    |remainder|<=E_0 j_s(ka)^4/12.                            (E2)

The corresponding long-wavelength effective mass and maximum free group
velocity are hbar^2/(2E_0 j_s a^2) and 2E_0 j_s a/hbar. Thus the supplied
j values imply v_g,A/v_g,B=3/2 and m_B*/m_A*=3/2 in this channel. These
are consequences of supplied hopping ratios, not derived particle masses
or a universal speed of light. The dispersion is quadratic at low k, not
relativistic. No preferred-time structure or path dimension disappears here.

## 2. Observable dictionary and scales

| Mathematical object | Operational meaning | Calibration / qualification |
|---|---|---|
| Recorded event x | Label at rank x of an accepted word | Recording is a supplied intervention |
| d=number of occupied links | Logical relational separation | Not automatically SI spatial distance |
| n_sx | Species-resolved occupation probability | Calibrated preparation/readout |
| X_sx,Y_sx,Z_sx | Occupation-tagged internal measurements and rotations | Defined only on the declared slot levels |
| chi | Derivative of a specified pulse/readout experiment | Nonzero source variance must be checked |
| tau | Dimensionless evolution time | tau=omega_0 t_SI, omega_0=E_0/hbar |
| E_0 | Common energy unit | Supplied physical matching/calibration, not a derived fundamental scale |
| a | Optional physical length per logical cell | Independently calibrated mapping; absent from primitive label law |
| E_ij | Model gauge register | Not identified with a measured electromagnetic field |
| Receiver mode index | Supplied auxiliary mode chain | Not an extra emergent spatial dimension |

In particular the holdouts below are **population probabilities after localized
injection**, not phase-pulse responses of a number-eigenstate source. The
finite-density phase protocol P12 instead uses an adjacent coherent source
with variance 1/4.

## 3. Exact fixed-time holdouts

Inject one carrier at the root, fixed internal level 0, with its compensator;
the other species is absent. No remote charged receiver is prepared. The
half-line sine transform of E1 gives

    A_s(r,tau)=(2/pi) integral_0^pi sin((r+1)k)sin(k)
                                      exp(2i j_s tau cos(k)) dk
              =i^r[J_r(2j_s tau)+J_(r+2)(2j_s tau)],
    p_s(r,tau)=|A_s(r,tau)|^2.                                (E3)

The time tau=32 and channels were frozen before calculation in
PREDICTIVITY_CHARTER.md. The predictions, rounded for display, are:

| Recorded outcome | Species A | Species B |
|---|---:|---:|
| Position 0 | 0.774578072 | 0.893901755 |
| Position 1 | 0.211244974 | 0.103142370 |
| Position 2 | 0.013778093 | 0.002919175 |
| Position >=3 | 0.000398861 | 0.000036701 |

PREDICTIONS.json encloses each target with outward decimal endpoints. This
particular numerical calculation uses **exact rational alternating-series
enclosures**, not an assumed double-precision error budget. At z=1 or 2/3,
the series

    J_r(z)=sum_(k>=0)(-1)^k (z/2)^(r+2k)/[k!(r+k)!]

has decreasing positive term magnitudes. The first omitted term bounds the
remainder. Twenty-five terms, rational arithmetic, interval addition and
squaring give the enclosures; endpoints are rounded outward to 50 decimals.
Independent spectral quadrature and 64-site matrix propagation agree at
ordinary double precision. They are implementation controls, not experiments.

For a reflecting chain of ell sites, a common-space boundary Duhamel bound
on the root-injected state is (1/2)tau_ell(2j_s|tau|). Thus the probability
error relative to E3 is at most tau_ell(2j_s|tau|). At ell=64,tau=32 it is
at most (z^64/64!)/(1-z/65): below 8.01*10^-90 for A and 4.28*10^-101
for B. These are mathematical boundary errors, not device-noise estimates.

## 4. Two different empirical claims

**Effective-channel experiment.** An engineered 64-cell path implementing
E1 can test these predictions. Its physical logical chain and control apparatus
are supplied. Agreement tests the effective quantum-transport model; it does
not validate D1's receiver, geometry formation, extensive preparation or any
claim that the universe has D1's microscopic structure.

**Full D1 conditional transfer.** If the actual D1 recorded/interior protocol
is implemented, its population prediction also pays 2Delta(32), from P7, in
addition to the finite-prefix tail. Its actual conditional receiver is retained.
The huge receiver resources, qualifying-record probability and successful
matter injection cannot be replaced by the easier effective-channel experiment.
For instance N=2048,ell=512,b=1,m>=1024 gives a very small evaluated Delta
at tau=32, but no formation guarantee at that N follows from S9. The current
S9 sufficient formation threshold remains N>=2^28, pending review.

Neither domain has experimental data in this packet. No apparatus or dataset
has been selected to make the predictions look favorable.

## 5. Calibration, uncertainty and a fixed future rejection rule

Use a separate early-time A calibration dataset to determine

    beta_A=lim_(t_SI->0) p_A(1,t_SI)/t_SI^2=j_A^2 omega_0^2,
    omega_hat=sqrt(beta_hat_A)/j_A.

Control extrapolation/calibration uncertainty before the holdout. This fixes
one time scale. Do not fit j_B/j_A, the populations at tau=32, or a new E_0
to the holdout. Independently characterized preparation/readout and Hamiltonian
errors are nuisance **bounds**, not new fitting coefficients.

For the effective-channel test require before data collection:

    |omega_0/omega_hat-1|<=epsilon_omega=0.001,
    ||H_device/E_0-H_target||<=epsilon_H=10^-5,
    preparation trace distance<=epsilon_prep=0.001,
    readout effect-operator error<=epsilon_read=0.001.

These are proposed eligibility requirements, not measured capabilities.
Run exactly 10,000 independent preparation/measurement trials per species
at physical time t_*=32/omega_hat. Record the four outcome categories without
discarding unfavorable ones. If the pretest bounds are not established, the
empirical test is not ready; do not widen them after seeing outcomes.

The probability allowance, before sampling noise, is

    epsilon_sys,s <=128j_s epsilon_omega+66epsilon_H
                      +epsilon_prep+epsilon_read+epsilon_boundary. (E4)

Indeed |dp/dtau|<=4j_s, while Duhamel gives at most
2tau_actual epsilon_H; tau_actual<=33 under the calibration bound. The
coefficient 66 safely includes their product, rather than dropping it.
E4 is 0.00466 for A and 0.003993333334 for B, plus the negligible boundary
tail. A full-D1 test must additionally include its comparison/acceptance terms.

For eight held-out channel probabilities and family error alpha=0.01,
Hoeffding plus a union bound gives

    epsilon_stat=sqrt(log(16/alpha)/(2*10000))=0.019206456.

**Reject the claimed calibrated model/device identification** if any observed
frequency lies farther from its certified target interval than
epsilon_sys,s+epsilon_stat. The family false-rejection probability is <=1%
only under the stated calibration and independent-trial assumptions. No
post-test refitting, time reselection, extra outcome category, or success-only
renormalization is allowed by this protocol. Invalid calibration leaves an
unresolved implementation question, not evidence of new physics.

These predictions are also produced by conventional tight-binding quantum
mechanics with the same Hamiltonian and resources. They are **not distinctive
GRUT predictions about nature**. They make a narrow candidate operationally
testable; a broader theory still needs independent validation and empirical
comparisons outside data used to choose its inputs.
