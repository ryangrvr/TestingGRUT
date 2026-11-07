# GRUT physical obligations and a photon-sector bridge

6 November 2026 · Companion to the [integrated theory](GRUT_Integrated_Theory.md)

The integrated candidate supplies a probability law once the physical algebra, state, Hamiltonian, and access are given. A physical extension has to constrain those inputs. The construction below realizes part of its continuum matching condition (14) and tests whether that match fixes a microscopic law. It uses the free electromagnetic field of the low-energy target on a fixed flat background.

| Obligation | Current answer |
|---|---|
| Select access and preparation | The integrated theory supplies both. The stationary-data obstruction excludes a deterministic, basis-equivariant selection of concrete proper factors from commuting data alone. In this sector, lattice geometry and gauge invariance specify the detector dictionary, and an explicit ground-state rule specifies the vacuum within the chosen flux sector. The geometry, detector permission, and use of that vacuum are inputs. |
| Recover the adopted continuum | An explicit family of finite-dimensional regulated models converges to the free Maxwell detector histories specified below. The interacting chiral Standard Model, a dynamical metric, its renormalized source, and the physical clock dictionary still require a match. |
| Impose an extra physical restriction | The present GRUT accepts both Hamiltonians below. A **plaquette-only magnetic energy** at a nonzero physical spacing would exclude one of them, even though both have the same gauge symmetry, access permissions, ground-state rule, and continuum limit. The restriction is conditional here. |
| Distinguish it with matched calibration | All axis-aligned vacuum probe histories coincide for the two models. An off-axis polarization frequency splitting distinguishes them at finite spacing. The splitting tends to zero with the regulator, so identifying a physical spacing is essential to the test. |

The proof of the stationary obstruction and the operational meaning of the detector probabilities are in the existing theory. This sector uses their declared scope. [Stationary obstruction](../grut_foundational_access_20261106/REPORT.md) · [Integrated theory, Sections 3 and 7](GRUT_Integrated_Theory.md)

## Physical algebra, preparation, and Hamiltonians

Use continuous time and a cubic spatial three-torus of side $L$, with $N^3$ sites, odd $N$, and spacing $a=L/N$. The torus is an infrared regulator. Work in units with $\hbar=1$ and keep the calibrated wave speed $c$. The noncompact link potential $A_i(x)$ and its canonical momentum $E_i(x)$ have the kinematic bracket $[A_i(x),E_j(y)]=i\delta_{ij}\delta_{xy}/a^3$ before reducing Gauss’s law. Define

$$
D_i^+ f(x)=\frac{f(x+a e_i)-f(x)}{a},\quad
D_i^- f(x)=\frac{f(x)-f(x-a e_i)}{a},\quad
B_i(x)=\epsilon_{ijk}D_j^+ A_k(x).
$$

The components of $B$ are on oriented plaquettes. The derivative along $i$ of $B_i$ is centered on a cube. The gauge transformation is $A_i\mapsto A_i+D_i^+\chi$ and physical modes obey $\sum_iD_i^- E_i=0$. The curl is gauge invariant because the forward derivatives commute. Work on the nonzero transverse modes, taking the zero flux sector and omitting the harmonic potential and electric zero modes. This removes the non-oscillatory degrees of freedom from the declared sector. The same choice is made in continuum Maxwell on the torus.

The admitted Hamiltonian family is

$$
H_{a,\eta}=\frac{a^3}{2}\sum_x\left[
 \sum_i E_i^2+c^2\sum_i B_i^2
 +\eta c^2 a^2\sum_i(D_i^+ B_i)^2
 \right],\qquad \eta\ge0 .
\tag{P1}
$$

Both magnetic terms are gauge invariant, invariant under cubic rotations, and local on the spatial regulator. The last term can be positive even though the discrete Bianchi identity gives $\sum_iD_i^+B_i=0$. The positive quadratic Hamiltonian preserves the Gauss constraint and has one normalized Gaussian vacuum on the nonzero physical modes. Select this vacuum for each Hamiltonian by the same minimum-energy rule. This is a preparation for the free-field comparison. It specifies neither the state of the cosmological target nor the gravitational effect of its vacuum energy.

The access dictionary uses smeared gauge-invariant $E$ and $B$ and the instruments below. The allowed microscopic generator terms include both magnetic terms of (P1). Hence $\eta=0$ would be a restriction on the actual dynamics within the same access table. For a choice of real smooth spatial smearing $f=(f_E,f_B)$, define $F(f)=\int(f_E\cdot E+f_B\cdot B)\,d^3x$ in the continuum and use the same Fourier components on the lattice, with link and plaquette midpoint phases. Couplings in $f$ make $F(f)$ dimensionless. The detector uses the Weyl unitary $W(f)=e^{iF(f)}$ and the two Kraus operators

$$
K_\pm(f,\varphi)=\frac{1\pm e^{i\varphi}W(f)}2 .
\tag{P2}
$$

They give a normalized instrument since $W(f)$ is unitary. The experiment family consists of a specified finite collection of ordered histories of these instruments at supplied inertial times in $[0,T]$. This tests the joint probability law, including disturbance by an earlier detector. It fixes the state class to the prepared vacuum and the displayed detector protocols. Gauge-invariant regional access is used without assuming independent spatial tensor factors on the constrained space. The lattice and inertial clock are given. The fixed flat metric places this calculation at the quadratic photon order of the adopted matter target, with charged matter and metric response omitted. [Gauge algebras and constraints](https://arxiv.org/abs/1312.1183) · [Low-energy gravity EFT](https://arxiv.org/abs/gr-qc/9512024)

## Observer and detector

The perspective of an observer is an operational specification. In the continuum target let a detector follow a timelike worldline $x_o(\tau)$ with proper time $\tau$, four-velocity $u_o^\mu$, and carried spatial axes $e_{o,i}^\mu$. In units with $c=1$ and metric signature $(-,+,+,+)$, an electric readout along one axis uses the gauge-invariant field component

$$
E_{o,i}=u_o^\mu F_{\mu\nu}e_{o,i}^\nu .
$$

The electric sign convention agrees with the canonical variable in the rest frame. A finite detector uses a smearing of this component, and its instrument and clock identify the corresponding probabilities in the integrated theory. Different velocities, orientations, and couplings can give different outcomes from one field state. In a common coordinate representation these are distinct physical protocols. A mere change of representation carries the state, dynamics, and instrument maps together. If $\Gamma$ is the corresponding algebra isomorphism, the transported state is $\omega\circ\Gamma^{-1}$ and each transported map is $\Gamma\circ\mathcal J_x\circ\Gamma^{-1}$ or $\Gamma\circ\alpha_t\circ\Gamma^{-1}$. The $\Gamma$ factors cancel in the composition in (1) of the integrated theory. The same physical experiment therefore has the same probabilities in either representation.

For a light wave in geometric optics with phase $\phi$, let $k_\mu=\nabla_\mu\phi$. A detector samples phase along its own worldline, so its measured angular frequency is

$$
\omega_o=-\frac{d\phi(x_o(\tau))}{d\tau}=-k_\mu u_o^\mu .
\tag{P0}
$$

For a wave moving in the positive $x$ direction in flat spacetime, an observer receding along $x$ at speed $v$ measures $\omega'=\gamma(1-v/c)\omega$ relative to an observer at rest in the original frame. This difference follows from contracting the same wave covector with two different four-velocities. A light ray has a null tangent and accumulates zero proper time. A rest frame and a proper-time clock require a timelike trajectory. The ray can be traced by an affine parameter, and its emission and arrival can be timed by detectors. The null interval specifies propagation and supplies no first-person clock or observer state for the photon. [Carroll, Lecture Notes on General Relativity](https://arxiv.org/abs/gr-qc/9712019) · [Einstein Online, Doppler effect](https://www.einstein-online.info/en/spotlight/doppler/)

The matching below can compare finitely many fixed inertial detector trajectories. Sample each detector at given proper times and express its smoothed electric components as linear combinations of the lab-frame smeared $E,B$ at the corresponding lab times. These are among the declared probes. The worldlines, axes and couplings determine the readout given the background. They supply no rule selecting the background or its state. The finite-spacing test fixes the detector and geometry used for the calibration and the comparison, so an ordinary Doppler shift is common to the two models on the calibrated sector.

## A finite family and its detector limit

With midpoint Fourier conventions the derivative multipliers have magnitude $\hat k_i=2\sin(k_i a/2)/a$, with the sign of $k_i$. Write $C_{\hat k}v=i\hat{\mathbf k}\times v$. A mode is physical on $\hat{\mathbf k}^\perp$. On that plane its frequency-squared matrix is

$$
Q_{a,\eta}(\mathbf k)
=c^2 C_{\hat k}^\dagger
 \big[1+\eta a^2\operatorname{diag}(\hat k_1^2,\hat k_2^2,\hat k_3^2)\big]C_{\hat k}.
\tag{P3}
$$

The first term is $c^2|\hat{\mathbf k}|^2$ times the identity on the two transverse polarizations. The Fourier symbol of Gauss's law removes the longitudinal oscillator. For each fixed nonzero $\mathbf k$, (P3) tends to $c^2|\mathbf k|^2$ on the continuum transverse plane as $a\to0$ at fixed $\eta$.

The spatial lattice has finitely many independent real oscillators. For a finite-dimensional candidate, additionally keep occupations $0,\dots,K$ for each of its diagonal oscillators. The projection $P_K$ commutes with the Hamiltonian and keeps its vacuum. Its range has dimension $(K+1)^{2(N^3-1)}$. Compress each Kraus operator to $P_KK_\pm P_K$. The sum of the compressed effects is at most the identity on this range. Add the positive square root of its deficit as an overflow Kraus operator, so the finite instrument is normalized. Use the full matrix algebra on the range and declare the compressed instruments there. This auxiliary number cutoff can be expressed in the eigenmodes of each Hamiltonian. Their uncut physical detector dictionary is common, and the cutoff agrees exactly on their common axial oscillators. The cutoff makes a finite probability model. Spatial commutation and tensor factorization of its compressed detectors are separate conditions.

**Number cutoff.** At fixed $a$, a word of finitely many free evolutions and Weyl unitaries sends the vacuum to a coherent state up to a phase. Its mean total occupation is bounded by the square of the sum of the displacement amplitudes in the word. The squared norm removed by $1-P_K$ is at most this mean divided by $K+1$, by the number distribution tail bound. Expand an uncut history branch into such words using (P2). Inserting the cutoff projections changes the branch vector by $O((K+1)^{-1/2})$, by telescoping the contractions. Squared norms give the probabilities. There are finitely many outcomes and the limiting no-overflow probabilities sum to one, so the overflow probability vanishes at the same order. Extend the uncut distributions by zero for overflow outcomes.

**Spatial cutoff.** On the Brillouin zone, $(2/\pi)|\mathbf k|\le|\hat{\mathbf k}|\le|\mathbf k|$ and $\eta a^2\hat k_i^2\le4\eta$. The positive frequencies are between $(2/\pi)c|\mathbf k|$ and $c\sqrt{1+4\eta}|\mathbf k|$. The error $\hat{\mathbf k}-\mathbf k$ is $O(a^2|\mathbf k|^3)$. Identifying the transverse planes by a rotation of size $O(a^2|\mathbf k|^2)$ shows that each frequency differs from $c|\mathbf k|$ by $O(a^2|\mathbf k|^3)$. This follows by bounding the stiffness-matrix difference by $O(a^2|\mathbf k|^4)$ and using the lower frequency bound to take square roots.

For a mode of frequency $\omega$, the electric covariance weight is $\omega/2$ and the magnetic weight has two curl factors and the position covariance $1/(2\omega)$. All field-strength covariance weights are therefore bounded by a constant times $|\mathbf k|$. The differences of their amplitudes are $O(a^2|\mathbf k|^3)$ after the transverse identification. The factors $e^{\pm i\omega t}$ contribute an additional error of order $a^2 T|\mathbf k|^4$ to the covariance weights. Products of smooth smearing coefficients make these bounds summable, and their high-momentum tails outside the lattice zone are $O(a^2)$ or smaller. Both real covariance and commutator pairings thus converge with error $O(a^2)$ for the declared probes on this fixed torus.

Products of Weyl operators reduce to the exponential of a sum of linear fields times the central commutator phase, since $[F(f),F(g)]$ is scalar. Their vacuum expectations are Gaussian and depend only on those covariance and commutator pairings. The joint instrument probabilities are finite sums of such products. Consequently for every fixed finite family of the displayed protocols,

$$
\operatorname{TV}\!\left(p_{\pi,\eta}^{(a,K)},p_\pi^{\rm Maxwell}\right)
\le \frac{C_{\pi,a,\eta}}{\sqrt{K+1}}+ C_{\pi,\eta,L,T}a^2 .
\tag{P4}
$$

The constants depend on the chosen smearing, couplings and history length. This gives ordered limits, occupation cutoff removed at each spatial spacing and then $a\to0$. The error covers joint probabilities at those detector resolutions. The convergence supplies the free flat photon part of the matching condition (14). It does not recover dynamical geometry, interacting matter, or an internal clock for the continuum. The usual use of a regulator is compatible with both members of the family, including arbitrary fixed positive $\eta$.

## A same-calibration discriminator at finite spacing

**Calibration.** Choose the probes supported on modes with $\mathbf k$ along a coordinate axis. For such a mode its magnetic field has zero component along that axis, and the other two derivatives along their own axes are zero. The extra term of (P1) vanishes for every axial polarization. The Hamiltonian and vacuum of the calibrated modes, hence their entire displayed detector histories, are identical for $\eta=0$ and any $\eta>0$. In particular, they have the same axial dispersion

$$
\omega_{\rm axis}^2(k)=\frac{4c^2}{a^2}\sin^2(ak/2) .
$$

Any calibration of $c$ and $a$ using this curve is shared. This statement treats the curve as calibration data of the finite model. No measurement fixing a physical spacing is supplied.

**Withheld direction.** Choose $\mathbf k=(q,q,0)$ with $q\ne0$ in the zone. Put $u=4\sin^2(aq/2)/a^2$. Compare the transverse polarization $e_z$ and the in-plane polarization $(1,-1,0)/\sqrt2$. The in-plane polarization produces only a $B_z$ component, so its extra energy vanishes because $\hat k_z=0$. The $e_z$ polarization produces $B_x,B_y$, each of squared magnitude $u|A_z|^2$. Equation (P3) gives

$$
\omega_{\parallel}^2=2c^2 u,\qquad
\omega_z^2=2c^2u+2\eta c^2 a^2u^2,\qquad
\boxed{\frac{\omega_z^2-\omega_{\parallel}^2}{\omega_{\parallel}^2}
 = 4\eta\sin^2(aq/2)} .
\tag{P5}
$$

Here the label $\parallel$ means in the plane and transverse to propagation. The conditional plaquette-only law $\eta=0$ predicts a vanishing split. The allowed comparator $\eta=1$ predicts a positive split with the same calibration. Both polarizations are read by the same inertial detector at the same lattice-frame wavevector. Changing coordinates leaves the comparison of its two readings invariant. A physical finite lattice would supply a rest frame and axes beyond the observers themselves. Covariance under a change of representation leaves this physical background in the comparison. Exact boost symmetry holds only in the local continuum Maxwell target here. It is stable because its stiffness is positive on the physical plane. This use of polarization dispersion is within established gauge-field kinematics. [Kostelecký and Mewes](https://arxiv.org/abs/0905.0031)

The frequencies are accessible to the declared detectors. For a normalized real oscillator mode, its electric momentum has symmetrized vacuum correlation $(\omega/2)\cos(\omega t)$. A two-read sign correlation with (P2), $\varphi=\pi/2$ and probe strengths $\varepsilon$ has leading term $\varepsilon^2\langle\{E(t),E(0)\}\rangle/2$. This follows by expanding the signed instrument $\sum_s s K_s^\dagger X K_s=-\varepsilon\{E,X\}/2+O(\varepsilon^2)$ and the signed second effect $-\varepsilon E(t)+O(\varepsilon^3)$. Normalizing by the zero-delay coefficient gives the cosine, so the splitting is a detector consequence rather than a gauge choice.

In the continuum limit the right side of (P5) is $O(a^2q^2)$ and vanishes. Thus even this detector-level continuum match, the same access and state rules, and an exact calibration on the axial sector leave $\eta$ undetermined. The plaquette-only law would need a physical selection principle and an independently identified finite scale and orientation for (P5) to become a prediction about nature. Including interactions would also require matching the coefficient in their presence. Both Hamiltonians are examples of ordinary quantum dynamics. The integrated GRUT currently imposes no law choosing one of them.
