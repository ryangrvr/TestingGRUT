# Mathematical foundations of the integrated framework

The following is a common specification for the established-physics part of GRUT. A model supplies objects of this form with its own validity domain. The type of specification is common even where the underlying Hilbert spaces, actions, constraints and protocols differ. No tensor product or interaction across the archived cases is defined by putting them in the same notation.

## 1. State, algebra and dynamics

Let $\mathfrak A$ be a unital physical observable algebra with state $\omega$ satisfying $\omega(I)=1$ and $\omega(A^\dagger A)\ge0$. A density matrix gives $\omega(A)=\operatorname{tr}(\rho A)$ in the finite case. A subsystem is represented by a specified inclusion $\mathfrak A_X\subseteq\mathfrak A$, or by an operationally declared instrument interface. The restricted state is $\omega_X=\omega|_{\mathfrak A_X}$. For finite tensor factors this is represented by a partial trace. Local continuum algebras and gauge constraints may prevent a canonical factorization. A geometric causal complement and an operator commutant can be identified only when the relevant duality and representation conditions hold.

One chooses the dynamical case appropriate to the model:

- **Ordinary closed dynamics.** Supply a self-adjoint $H$ or another admissible specified reversible dynamics. With a calculational time parameter $t$, $\rho(t)=e^{-itH/\hbar}\rho_0e^{itH/\hbar}$. The formal parameter alone is not a physical clock observable. Sectors, symmetries and conserved charges follow from the actual generator and operator domains.
- **Constrained dynamics.** Supply kinematic variables, constraints $\widehat C_a$, a representation, a physical inner product and an admissible observable algebra. For the finite cyclic F5 case, the physical space is the image of a group-averaging projector. General constraints can require a rigging or reduction construction. A first-class classical algebra, $[\widehat C_a,\widehat C_b]$ in a quantum representation, physical positivity and control of observables are separate issues.
- **Effective dynamics.** Specify which variables are retained, the joint state and preparation class, approximations, scale and error estimates. A reduced channel may be CPTP on arbitrary allowed inputs when an appropriate preparation assignment such as a common product reservoir exists. Restriction of a single global evolution always defines local expectations on that global state, while arbitrary continuations from a local marginal can depend on unavailable correlations.

For example, an exact time-independent projection identity for a finite linear evolution $\dot \rho={\cal L}\rho$ with projectors $P,Q=I-P$ gives
\[
\frac{d}{dt}P\rho(t)=P{\cal L}P\rho(t)
  +\int_0^t ds\, P{\cal L}Q\,e^{(t-s)Q{\cal L}Q}Q{\cal L}P\rho(s)
  +P{\cal L}Q\,e^{tQ{\cal L}Q}Q\rho(0).
\]
Here the exponentials act on the $Q$ space. This follows by solving the $Q$ equation and substituting it into the $P$ equation. The memory kernel, initial-correlation term, spectral assumptions and approximations used to replace them by local response are properties of that model. A stationary Gaussian complement yields different simplifications from a driven or reused one. The identity does not choose the operator $P$, the generator or the preparation. [Closed-whole formal core](../grut_closed_whole_theory_20261019/FORMAL_CORE.md) · [I1 reduction](../grut_i1_integrated_20261125/EFFECTIVE_THEORY.md)

## 2. Interventions, clocks and records

A quantum instrument consists of completely positive outcome maps with trace-preserving sum in the Schrödinger picture. Its Heisenberg adjoint maps the identity to the effects $E_x$, with $\Pr(x)=\omega(E_x)$. To predict a sequence include the appropriate intervening dynamics and the adjoints of the actual instruments. Composing maps preserves the order and distinguishes conditioning on a record from discarding it. The measurement statistics postulate is adopted.

For one clock-and-record procedure let $E_{\tau,r}$ be the joint physical effect, and $D_\tau=\sum_r E_{\tau,r}$. The conditional record law is
\[
p(r\mid \tau)=\omega(E_{\tau,r})/\omega(D_\tau),\qquad \omega(D_\tau)>0.
\]
A product of readout effects is justified when it describes the joint instrument. Gauge-invariant compression of a joint effect and multiplication of separately compressed effects can be different operations. In well-defined complete clock charts, $U$ transports the state and all effects, so the numerator and filter probabilities are invariant. Conditioning on unmatched readouts, tracing out a reference or changing the calibration is a change of observable data. Approximate matching has a finite operator-norm bound that deteriorates when the filter is rare. [F5 foundations and bound](../grut_f5_closure_consistency_20261127/DERIVATIONS.md)

In a finite model a record can be written by a reversible interaction acting on a ready memory $M$. Its joint states, accessible effects, conditional decay and readout fidelity specify what can be learned. The record need not contain the whole preceding state. A reservoir occupation can record a defect count without retaining the original arrangement in that one observable. A restricted record may lose the phase or inter-subsystem correlations needed for a prediction. A pointer effect and a causal use of its result must be established separately. For observers with limited access, a minimal linear summary for a given collection of future probabilities is described by the span of their backward-evolved effects in a declared state class. Changing the future experiment changes that span. [Foundational access](../grut_foundational_access_20261106/REPORT.md) · [D6 causal controls](../grut_d6_adaptive_20261125/CONTROLS.md)

## 3. From microscopic coupling to response and EFT

A probe with perturbation $\delta H(t)=j(t)B$ has, to first order and with that sign convention, response kernel
\[
\frac{\delta \langle O(t)\rangle}{\delta j(t')}
= -\frac{i}{\hbar}\theta(t-t')\,\omega([O(t),B(t')]).
\]
The response depends on the supplied representation, dynamics, state and physical probe. Calibrated correlation functions can constrain dispersion and geometry, but access restrictions can identify only an equivalence class. Ward identities, positivity, anomaly cancellation and symmetries restrict admissible responses in their stated domains. They do not provide a single action by themselves. [F1 and F2 geometry](FOUNDATIONAL_OBSTRUCTIONS.md) · [F4 consistency audit](../grut_f4_charge_selection_20261126/ESTABLISHED_CONSTRAINTS.md)

A representative covariant low-energy action used by the matter and gravity lanes has the form
\[
S_{\rm EFT}=\int d^4x\,\sqrt{-g}\left[
\frac{M_P^2}{2}R-\rho_\Lambda+\mathcal L_{\rm SM}
+\mathcal L_{\rm added}+\sum_i c_i\,\mathcal O_i/\Lambda^{d_i-4}
\right]+S_{\rm boundary}.
\]
The field content, representation, gravitational coupling, finite coefficients, counterterm scheme, preparation and boundary conditions require physical matching. Familiar gravitating field dynamics, including $\nabla_\mu T^{\mu\nu}=0$ on shell in the applicable covariant class, follow within the chosen EFT. F3 shows that even the full matched homogeneous stress history and ordinary classical constraint closure leave effective matter characteristics undetermined. Low-energy mode stability and a common causal slicing are additional checks. They do not amount to a nonlinear or nonperturbative quantum-gravity construction. [F3 model and constraints](../grut_f3_coupled_consistency_20261126/MODEL_AND_CONSTRAINTS.md) · [Candidate 1](CANDIDATE_1_STATUS.md)

A bridge from a finite model to an EFT needs a common physical preparation, a matching map for operators and couplings, and controlled predictions for the same observables. A fitted frequency, trace, or entropy by itself can leave many microscopic theories consistent. In the pinned archive, I1 supplies a specific controlled open-system limit. Detector matching, relativistic locality, a full physical state dictionary for cosmological sources and a continuum gravitational quantum limit have separate obligations.

## 4. Thermodynamics and organization

For a closed unitary interaction starting in $\rho_A\otimes\tau_B$ with $\tau_B=e^{-\beta H_B}/Z$, the standard finite-reservoir identity is
\[
\Delta S_A-\beta Q_A
=I(A:B)_{\rm f}+D(\rho_{B,\rm f}\Vert\tau_B)\ \ge0,\qquad Q_A=-\Delta E_B .
\]
Use dimensionless entropy and $\beta=1/(k_BT)$. Global entropy conservation gives $\Delta S_A+\Delta S_B=I(A:B)_{\rm f}$. The Gibbs formula gives $D=\beta\Delta E_B-\Delta S_B$. Interaction, control and work terms require their own consistent endpoint accounting. A non-Gibbs or initially correlated complement requires a different accounting. This is not a consequence of closure without preparation hypotheses. Entropy and work of sensing, correlations, fair preparation, resetting and exhausting a reservoir must be compared jointly. [Closed-whole thermodynamic bounds](../grut_closed_whole_theory_20261019/FORMAL_CORE.md) · [I2 work accounting](../grut_i2_persistence_20261125/THERMODYNAMICS.md)

An organizational test specifies a physical order observable or a set of acceptable outputs and a scoring history. An operational historical label may be an initial codeword, a domain-wall count or a matching rule. A persistence test can concern a first exit, a finite readout window or a final configured boundary. Repeated measurement and unchecked evolution define different tasks. Information acquired through one interaction can influence a later action if the physical ability to use it is provided. It cannot label histories that leave the complete available state identical beyond prior inference. The resulting organization and predictive advantages are effective properties of the dynamics, access and resources. [Developmental definitions and mechanisms](DEVELOPMENTAL_EFFECTIVE_THEORY.md)

## 5. Bridges still missing

| Boundary | Needed for a stronger synthesis |
| --- | --- |
| Closure to laws or an initial state | An independently warranted selector or admissibility restriction demonstrated before conditioning on the target data. F5 establishes freedom in its supplied cyclic class. |
| Nominal parts to physical access | A physical dictionary for the relevant observable algebras and clocks, not a preferred split inferred from a symmetry-invariant stationary state alone. |
| Relational networks to spacetime | A local observable net and controlled propagation and continuum limits in a physically supported mechanism. Inverse recovery of effective connectivity needs independent calibrations. |
| Matter EFT to cosmology | A physically prepared continuous history, matched stress and remaining work, and a causal transfer and abundance prediction with calibrated errors. |
| Information to a new physical law or consciousness | An independently testable mechanism or mapping. The archive establishes ordinary records, control and selective access only in its physical instruments. |
