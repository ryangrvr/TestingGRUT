# Unchanged-law diagnostic identities

Author derivations for the pinned finite B2 matter/fuel law. Independent
mathematical review remains pending. No full growing-geometry theorem is
obtained by studying these fixed records.

## 1. Sector and spatial continuity

q=N_+=N_-=2-F_L-F_R counts opposite-charge matter quanta. It does not assert
that every constituent is co-located in a composite. Let P_q project onto
q=0,1,2 and psi_q=P_q psi. Only conversion connects these sectors. With

    z_01=<psi_1,H psi_0>, z_12=<psi_2,H psi_1>,
    I_01=2 Im z_01, I_12=2 Im z_12,

the exact equations are p0'=-I01, p1'=I01-I12, p2'=I12.
Complex z, including its real part, is retained in evaluation; populations
are not propagated as a classical Markov chain. Hopping and all virtual
charged separations remain inside each psi_q.

Fuel patterns are mutually exclusive. P(F_L=F_R=1)=p0, whereas

    <F_L>=p0+P(F_L=1,F_R=0),
    <F_R>=p0+P(F_L=0,F_R=1).

The old 70.46% was EACH fuel's marginal, not the joint both-full event.
The accepted R is nested in q=2 and is not multiplied by these probabilities.

For an endpoint e define M_e=(n_+,e+n_-,e)/2 and L_e=M_e+F_e. Conversion,
density, field, mass, fuel and counterterm operators commute with L_e.
Only hopping across that endpoint's bond changes it. Consequently

    d<L_e>/dt=2 Im <L_e psi,H_hop psi> = -J_e,out,
    Delta F_e = Delta L_e - Delta M_e.                    (1)

Over theta in [2,6], Delta L_e is NET resource returned across the bond.
Delta M_e accounts for changed matter already at the source. Equation (1)
is an exact quantum continuity identity, not a particle-history assignment
or additive decomposition of hypothetical independent mechanisms.
It determines whether local conversion alone could explain fuel refill.

As a diagnostic baseline, delete hopping mathematically but keep Omega,
lambda and the original A counterterm fixed. End cells then factor exactly
and the one-cell matter probability is

    p_iso(theta)=576/745 sin^2(sqrt(745)*theta/24).

The two-cell sectors are binomial for this REFERENCE law only. This is not
a physical change or replacement operating protocol. Agreement at one
endpoint would not establish the same dynamics throughout the interval.

## 2. Low spectral cluster of the complete two-resource matrix

Let P span both fuels, one fuel plus one co-pair anywhere, and two
nonadjacent co-pairs. On eight events its dimension is 1+16+21=38.
Every P diagonal energy is <=4 and ||P H_convert P||<=2J, so
sup PHP<=4+2J. Every other bare configuration has either adjacent co-pairs
or a nonzero electric string. For this eight-event law its diagonal is
at least 9/2-A. The complete hopping norm is <=tau2=5/48 and conversion
norm <=2J. Thus

    QHQ >= h_star=9/2-A-tau2-2J >4+2J.

Min-max gives exactly 38 low eigenvalues below this gap, and every remaining
eigenvalue is >=h_star. This is a fixed-record MATTER spectral cluster,
not a claim of a storage-consistency gap in the full receiver theory.

Initial double fuel has mean E0=2Omega and variance 2J^2, including both
conversion amplitudes. Its exact spectral measure therefore satisfies

    w_high <=2J^2/(h_star-E0)^2=32/37027225 <8.65e-7.       (2)

This weight is conserved at every time. Virtual charged configurations
can still be present in a low eigenstate: a high ENERGY BAND weight is
not a projector onto bare separated-charge states. No such state is
discarded. The normalized low-band evolution differs from the full state
by norm <=sqrt(2 w_high), so any event changes by <=2sqrt(2 w_high).
Large endpoint failures cannot be attributed solely to populating the
high-energy spectral bands. The actual low-band dynamics remains to diagnose.

## 3. Spectral interference and recurrence

In finite volume all spectra are discrete. Write psi(t)=sum_alpha
c_alpha exp(-iE_alpha t)v_alpha, keeping every mode. The initial both-fuel
state is unique, so its return amplitude is exactly

    a_ff(t)=sum_alpha |c_alpha|^2 exp(-iE_alpha t).          (3)

For any projector X, p_X(t) includes both spectral diagonal and interference
terms. The diagnostic computes eigenvalue intervals from the pinned
residual/orthogonality account. Adjacent overlapping intervals are grouped;
within each group phases remain coherent. At the SAME endpoint, compare
the true p_X with sum_g ||X psi_g(t)||^2. Their difference is resolved
inter-group interference, not a replacement physical dephasing mechanism.
A finite-time average is separately computed with the exact sinc kernel
exp[-i(Ea-Eb)T/2] sinc[(Ea-Eb)T/(2pi)] over [0,T]. It is not substituted
for the failed endpoint or conflated with the spectral-block comparison.

Dominant weights, source and accepted-output support are reported for the
prepared state's modes. Source-heavy finite modes are not automatically
infinite-volume trapped states. Wave-packet motion is diagnosed by spatial
currents. A proof of permanent trapping or escape would need a declared
infinite-volume limit; this finite diagnosis does not supply one.

## 4. Interior fuel: channel and recoverability

Use the same three-event interior-fuel fixture from B2, initially co-pairs
on events 0 and 2 with empty middle fuel. It is a conditional information
fixture, NOT production from the original fuel seed. In the declared
configuration/color encoding, H=A tensor I+B tensor S, S=SWAP. Define
P_T=(I+S)/2, P_S=(I-S)/2, H_T=A+B, H_S=A-B and

    phi_T=exp(-itH_T)psi_spatial,
    phi_S=exp(-itH_S)psi_spatial,
    gamma=<phi_S,phi_T>.

If the spatial/fuel/field environment is discarded after collecting the
two ordered color carriers, the exact channel is

    E(rho)=P_T rho P_T+P_S rho P_S
              +gamma P_T rho P_S+gamma* P_S rho P_T.       (4)

This is an idealized stronger color access than merely reading two outgoing
pairs; it can include a color carried by fuel. Its collective algebra
B(C^3_triplet) direct-sum C_singlet is preserved exactly. One logical qubit
or a qutrit encoded wholly in the triplet remains intact. Arbitrary two
qubits also need coherence between the sectors.

If |gamma|=1, a known collective exchange-phase unitary corrects (4).
If |gamma|<1, a deterministic color-only recovery cannot restore all inputs:
the originally orthogonal |01>,|10> have output trace distance |gamma|,
which a subsequent channel cannot increase back to 1. Information is
correlated with the discarded registers, not destroyed globally. Correcting
arg gamma gives reference-entangled-input fidelity (5+3|gamma|)/8; this
number is an explicit correction's fidelity, not a claim of optimal recovery.

For return to the original spatial configuration r (two outer co-pairs,
empty fuel), the actual successful Kraus operator is

    M_r=a_r P_T+b_r P_S,
    a_r=<r|phi_T>, b_r=<r|phi_S>.

Success probabilities are |a_r|^2 on triplet inputs and |b_r|^2 on singlet
inputs. They must not be replaced by a qubit-independent arrival number.
Known joint phase correction removes arg a_r,arg b_r but leaves a filter.
Conditional fidelity for a maximally entangled 4-dimensional input after
that correction is

    (3|a_r|+|b_r|)^2 / [4(3|a_r|^2+|b_r|^2)].            (5)

The complete environment-access map is an isometry and has a perfect
inverse on its image; physically reversing the entire dynamics or an
equivalent coherent decoder is an additional resource, not the original
receiver algebra. Measuring the full spatial configuration r and applying
a configuration-dependent exchange phase plus inverse filter gives a
HERALDED exact identity subchannel with total original-trial probability

    p_herald=sum_r min(|<r|phi_T>|^2,|<r|phi_S>|^2).       (6)

This is an explicit possible recovery construction, not an implemented
autonomous protocol, an optimized decoder, or a rescue of B2's failed
outgoing event. A coarse fuel-occupation record generally leaves further
unobserved position/field correlations, which are separately accounted.
