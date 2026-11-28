# Complete mathematical framework

This appendix specifies the working interfaces and gives the reductions on which the synthesis relies. Algebraic quantum mechanics in the declared representations, instruments, linear response, Gibbs thermodynamics, Einstein theory and EFT enter as conventional physics. The task-specific Hamiltonians, counterexamples, operational comparisons and validated scalar certificates are GRUT constructions. Each proof uses its own preparations and validity regime. [Manuscript](GRUT_WORKING_THEORY_FINAL.md) · [Reproduction index](REPRODUCTION_INDEX.md)

## 1. Physical datum and reduction

In a quantum example the datum is a positive normalized $\omega$ on a physical observable algebra $\mathfrak A$, its dynamics or constraints, a family of accessible algebras $\{\mathfrak A_X\}$ and the preparations and instruments. In a finite representation a valid instrument $\{\mathcal I_r\}$ has completely positive maps, $\sum_r\mathcal I_r$ trace preserving, and effect $E_r=\mathcal I_r^\dagger(I)$. After a channel $\mathcal E$, $\Pr(r)=\operatorname{tr}[\rho\mathcal E^\dagger(E_r)]$. A history of accessible controls and readouts similarly pulls its joint event back to an effect on the initial state. The real span $\mathcal V$ of those effects, with identity included, gives exact predictive equivalence for that task. Two states agree on all such histories precisely when their difference pairs to zero with every element of $\mathcal V$. Discarding a joint or phase-sensitive effect can merge preparations of distinct predictive consequence. [Conventional mathematical foundations](../grut_s1_integrated_synthesis_20261127/MATHEMATICAL_FOUNDATIONS.md) · [Finite access result](../grut_foundational_access_20261106/REPORT.md)

With a bounded time-independent full generator $\mathcal L$ and a projection $P$ onto retained operators or states, let $Q_P=1-P$. Solve the complementary block equation by variation of constants, then substitute it in the retained block. This gives
$$
\begin{aligned}
\partial_t(P\rho)(t)
={}&P\mathcal L P\rho(t)
+P\mathcal L Q_P e^{tQ_P\mathcal L Q_P}Q_P\rho(0)\\
&+\int_0^t P\mathcal L Q_P e^{(t-s)Q_P\mathcal L Q_P}
         Q_P\mathcal L P\rho(s)\,ds .
\end{aligned}
$$
Energy and information can flow between retained and eliminated variables of a closed whole. The memory kernel and initial-correlation term depend on the split and preparation. For $H\mapsto H+j(t) B$, first-order retarded response in the unperturbed dynamics is
$$
\frac{\delta\langle O(t)\rangle}{\delta j(t')}=
-\frac{i}{\hbar}\Theta(t-t')\omega([O(t),B(t')]).
$$
That is the conventional response whose access sensitivity is tested by F1 and F2. For any observation map $R:\mathcal X\to\mathcal D$, a target $L$ is determined by the data iff it is constant on every fiber of $R$: necessity follows by composition, sufficiency by defining $f(d)=L(x)$ for any preimage of $d$ so $L=f\circ R$ on the image. Calibration and covariant transport are part of $R$. [Whole-system formal core](../grut_closed_whole_theory_20261019/FORMAL_CORE.md) · [F1 proof](../grut_f1_operational_reconstruction_20261126/DERIVATIONS.md)

## 2. Relational states, time and events

For F5 let $a,b\in\mathbb Z_N$ with $N$ odd prime, $s\in\{0,1\}$, $\nu_s=1+s$, $r=b-\nu_s a$. Given unitary prefixes $V_{s,r,a}$ with $V_{s,r,0}=I$, the shift maps $(a,b)\to(a+1,b+\nu_s)$ and acts on matter by $V_{s,r,a+1}V_{s,r,a}^\dagger$. Thus the wrap compensates the product of the preceding gates, and $T^N=I$. Its invariant-sector projector and seed isometry are
$$
P_0=N^{-1}\sum_{j=0}^{N-1}T^j,\qquad
J|s,r,\psi\rangle=N^{-1/2}\sum_a|a,r+\nu_sa,s\rangle V_{s,r,a}|\psi\rangle.
$$
Orthogonality at fixed $a$ gives $J^\dagger J=I$. Invariance determines the whole vector from the $a=0$ components, so $JJ^\dagger=P_0$. Reductions $\sqrt N\,\langle a|_A J$ and $\sqrt N\,\langle b|_B J$ are unitary onto the complete remaining spaces because each clock label determines the other given $(s,r)$. The latter uses $a_*=(b-r)/\nu_s$. A kinematic effect $F$ induces $F_{\mathcal L}=J^\dagger FJ$, equal to the compression of the gauge average of $F$. In chart $X$ use $U_X\rho U_X^\dagger$ and $U_X F_{\mathcal L} U_X^\dagger$. Use this rule for the joint numerator and the filter denominator. Ratios use only filters of nonzero weight. The frozen example has joint weight $1/20$ and filter weight $3/40$, with conditional $2/3$ in both charts. An independent tick comparison yields $1$ versus $1/2$ for distinct effects. [F5 construction](../grut_f5_closure_consistency_20261127/MATHEMATICAL_FOUNDATIONS.md) · [F5 derivations](../grut_f5_closure_consistency_20261127/DERIVATIONS.md)

Any finite unitary sequence fits into a sufficiently long cycle in $r=0$. In another relative sector a single $-I$ prefix gives a residual sign submatrix with rows $(-1,1)$ and $(1,1)$, of determinant $-2$. Thus nonseparable clock dependence and an unrestricted supplied history can coexist with exact probability transport. A laboratory calibration to proper time, dynamical locality restrictions and the mechanism preparing the seed remain separate physical premises. [F5 obstruction](FOUNDATIONAL_OBSTRUCTIONS.md#f5--matched-clock-events-do-not-select-dynamics)

## 3. Records, organization and repair

D1 uses $\mathrm{CNOT}_{S\to M}$, then $\mathrm{CNOT}_{M\to A}$, then $\mathrm{CNOT}_{A\to S}$. On basis input $(s,0,0)$ these steps yield $(s,s,0)$, $(s,s,s)$ and $(0,s,s)$. By linearity an arbitrary input moves into correlations of $M,A$, giving one ebit for a balanced pure qubit. A density matrix restricted to one of those registers discards the relative phase. This is reversible with the joint registers and uses a prepared set of controls. [D1 evidence](../grut_i1_integrated_20261125/EVIDENCE_MAP.md)

D2 constructs an autonomous, local finite spin model. For physical labels $s_i$, the orientation is $r=s_0$ and wall numbers come from $u_i=s_i\oplus s_{i+1}$. The inverse recursively determines every $s_i$, so the basis change is unitary. Adjacent wall exchange is $gX_{i+1}(1-Z_iZ_{i+2})/2$. Including the reservoir, the wall-plus-reservoir open path has exact fermion form
$$
H=\Delta Q+g\sum_{j=0}^{m+L-2}(c^\dagger_{j+1}c_j+c^\dagger_jc_{j+1}),\quad Q=\sum_j c_j^\dagger c_j .
$$
The path hopping operator has single-particle norm at most $2g$. Thus $H\ge(\Delta-2g)Q$ with the supplied $\Delta=4g$. In the selected sector identity of the orientation persists, the empty-wall state is the ordered state, and the conserved count can be retained in the prepared reservoir. Finite windows follow from exact transport and the leakage bound in the proof. For the joint order and count claim the diagnostic reference must be carried in the compared input states. [D2 proof](../grut_d2_development_20261123/PROOF.md) · [Corrected control](../grut_d2_development_20261123/CONTROL_CORRECTION_01.md)

At a repair interface, a subsequent label-identifying measurement $\{E_b\}$ induces a POVM $\{\mathcal E^\dagger(E_b)\}$ on the available inputs. Therefore an average recovery score for distinguishable target labels is bounded by the optimal available-state discrimination. In the equal-prior $d$-label case with states on a $D$-dimensional space, $\rho_b\le I$ gives $\sum_b\operatorname{tr}(F_b\rho_b)/d\le D/d$. Identical available states erase the distinction at that cut. D3 realizes flagged erasure of the substrate with a retained record and conditional reconstruction. D4 uses a substrate and two encoded copies. One wrong copy has a recoverable majority in its separated fault window. Two flips can give the opposite majority and three identical flips produce the other codeword, indistinguishable to a decoder with only those inputs. Bounds for repeated separated faults rely on drain capacity and time margins. [D3 proof](../grut_d3_record_repair_20261124/PROOF.md) · [D4 proof](../grut_d4_internal_recognition_20261124/PROOF.md)

## 4. Development, acquisition and quantum function

D5 writes an initially blank tape from an acquired bit $M$ and a generative bit $G$. At the exact completion edge its intended symbol is $A_i=M\oplus [G(i\bmod2)]$ for valid inputs. The history register is a physical head. With work gates $U_j$ define $V_0=I$, $V_{j+1}=U_jV_j$ and $\mathcal D=\sum_j|j\rangle\langle j|\otimes V_j$. The frozen coupling $\omega\sqrt{(j+1)(L-j)}/2$ dresses the spin-$L/2$ chain so that $H=\mathcal D(\omega J_x\otimes I)\mathcal D^\dagger$. The propagated seed is $\sum_j a_j|j\rangle V_j|s\rangle$ with
$$
a_j(t)=(-i)^j\sqrt{\binom Lj}\cos^{L-j}(\omega t/2)\sin^j(\omega t/2).
$$
At $T=\pi/\omega$ this is exact completion. The full one-head dynamics recurs up to common phase at $2T$. Single-pulse bounds use concentration on an idle band where $V_j$ is constant. If the band has weight $q$, a pulse on the full tape differs in state-vector norm from the band-constant comparison by at most $2\sqrt{1-q}$. The recurrence, calibration and bypass all limit interpretation. [D5 theorem](../grut_d5_prospective_20261124/PROOF.md) · [D5 limits](../grut_d5_prospective_20261124/LIMITS.md)

D6 probes two equiprobable programmed interaction signs with output states $\rho_0,\rho_1$. Maximizing the binary decision score over $0\le E\le I$ chooses the positive spectral projector of $\rho_0-\rho_1$. Therefore
$$
p_{\rm success}=\tfrac12+\tfrac14\|\rho_0-\rho_1\|_1 .
$$
For the registered probe overlap of magnitude $\cos2\theta$ on the calibration interval, the norm is $2\sin2\theta$. The sample preparation and post-acquisition access determine which score can be used at the later actuator. [D6](../grut_d6_adaptive_20261125/PROOF.md)

For D7 and I1, let $b$ have weights $\pi_b$, conditional accessible programme $\rho_b$, and qubit coherence factor $\alpha_b$ with common population flips $u=\Pr(1|0)$ and $v=\Pr(0|1)$. With $\lambda=1-u-v\ge0$, $\bar\rho=\sum_b\pi_b\rho_b$ and $C=\sum_b\pi_b\alpha_b\rho_b$, a Bell input separates into diagonal reference blocks and an off-diagonal programme block $C/2$. For Stinespring blocks $V_i=\langle i|_O V|i\rangle_Q$ the Bell fidelity is
$$
F=\tfrac14[u+v+\lambda\operatorname{tr}\bar\rho(V_0^\dagger V_0+V_1^\dagger V_1)
        +2\operatorname{Re}\operatorname{tr}(V_1^\dagger V_0 C)] .
$$
The blocks are contractions. The polar unitary of $C$ can be implemented conditionally while leaving the computational populations unchanged, giving $F_{\rm opt}=(1+\lambda+2\|C\|_1)/4$. Discarding the programme makes the coherence term $|\operatorname{tr} C|$. D7 has $\lambda=1$ and $\alpha_b=e^{-2i(-1)^b\phi}$. I1 adds prepared thermal contacts and distinguishes the theoretical optimum from a restricted actuator. Its zero usable gain at registered jitter is compatible with strictly positive ideal gain. [D7](../grut_d7_functional_recovery_20261125/PROOF.md) · [I1 channel and energy derivations](../grut_i1_integrated_20261125/DERIVATIONS.md)

## 5. Thermodynamics, memory value and reduced observability

For an initially uncorrelated Gibbs bath, total entropy conservation gives $\Delta S_A+\Delta S_B=I(A:B)_f$. With $D(\rho'_B\|\tau_B)=\beta\Delta E_B-\Delta S_B$ and $Q_A=-\Delta E_B$, these yield the manuscript's entropy balance. A work or reset budget still requires an explicit bath, Hamiltonian and allowed policy. I2's first-exit costs apply to its frozen menu, not to arbitrary storage designs. [I2 ledger and mathematics](../grut_i2_persistence_20261125/DERIVATIONS.md)

In I3 a symmetric sign flips at rate $q$. Solving $\dot p=q(1-2p)$ for retention gives $(1+e^{-2q\tau})/2$. A measurement with error $\epsilon$, memory flip rate $\mu$ during its unshielded time $\tau-\delta$, and transfer error $\eta$ gives
$$
p_A=\tfrac12+\tfrac12(1-2\epsilon)e^{-2q\tau-2\mu(\tau-\delta)}(1-2\eta).
$$
If the common action produces $s_h$ when initially matched and $s_l$ when mismatched, incident prediction must have $p_A\ge(r_\star-s_l)/(s_h-s_l)$ for target $r_\star$ when $s_h>s_l$. For the screened null a no-hit event has no new material observation until binding. A rate constraint and convexity of $e^{-x}$ bound relaxed material-feedback success by $1-e^{-\pi\overline c T_r}$, where $\pi$ here is the birth fraction in the bond dynamics. This is an interface theorem. A cue before binding can invalidate that null argument. The archived scores compare an acquired record with screened policies, a weak cue and a perfect present cue using the declared cost ledger. [I3 analytical gate](../grut_i3_predictive_memory_20261126/FEASIBILITY.md) · [I3 results](../grut_i3_predictive_memory_20261126/RESULTS.json)

## 6. Geometry and dynamics from partial physical access

For unit-mass springs with $K=I+L$ and calibrated port embedding $B$, expansion of $R(z)=B^\mathsf T(K+zI)^{-1}B$ at large $z$ determines all moments $B^\mathsf TK^jB$. The Gram matrix of $K^iB e_\alpha$ is composed of these moments. An equality of responses defines an isometry of the two Krylov spans that intertwines their stiffness matrices and ports. With full visibility it is a full orthogonal transformation. In the common-pin connected class the constant vector is the unique pin eigenvector with its sign fixed at the port. The hidden block therefore fixes its constant vector, and at most two hidden sites admit only permutations. Three hidden sites allow a weighted cycle and a pendant triangle with all root responses equal. The exact hidden block, fixing site 0, is
$$
\frac13\begin{pmatrix}2&-1&2\\-1&2&2\\2&2&-1\end{pmatrix}.
$$
Apply the same map to coordinates and canonical momenta. Observable calibration of the hidden sites can invalidate that equivalence. The continuum conformal obstruction instead follows from $P_{\Omega^2g}(\Omega^{-1}\varphi)=\Omega^{-3}P_g\varphi$ for $P_g=\Box_g-R_g/6$ in four dimensions. With the conformal factor one where the causal source/readout is calibrated, selected retarded probe responses have exact blind directions. [F1 proof and witness](../grut_f1_operational_reconstruction_20261126/DERIVATIONS.md)

For F2 set $G=(K+zI)^{-1}$ at fixed $z>0$. Connected positive spring networks give $u_i=G_{0i}>0$. A singleton loading $d>0$ causes root loss $r_\varnothing-r_i(d)=d u_i^2/(1+dG_{ii})$. Two such values determine $u_i$ and $G_{ii}$. Write $\alpha=G_{ii}+1/d_i$, $\beta=G_{jj}+1/d_j$. The pair loss as a function of $g=G_{ij}$ is
$$
F_{ij}(g)=\frac{\beta u_i^2+\alpha u_j^2-2gu_i u_j}{\alpha\beta-g^2}.
$$
The inverses of principal spring submatrices are entrywise nonnegative. Schur elimination then gives $0\le g\le\min(G_{ii}u_j/u_i,G_{jj}u_i/u_j)$. On that interval
$$
F'_{ij}(g)=-\frac{2(\beta u_i-gu_j)(\alpha u_j-gu_i)}{(\alpha\beta-g^2)^2}<0 .
$$
Thus a pair setting recovers one cross entry, and the complete data recover $K=G^{-1}-zI$ using $1+2(n-1)+\binom{n-1}{2}$ settings. Tags and load strengths are inputs. Poor conditioning under weak access is compatible with injectivity. Singleton all-frequency controls and nonconformal causal tests have the separate scopes in the source. [F2 spring theorem](../grut_f2_complementary_interactions_20261126/SPRING_DERIVATIONS.md) · [Relational and cosmological evidence](FINAL_EVIDENCE_MAP.md)

For F3 take the adopted Einstein action with matter $P(\phi^I,X^{IJ})$, where $X^{IJ}=-g^{\mu\nu}\partial_\mu\phi^I\partial_\nu\phi^J/2$. For scalar normal velocities $y^I$ and spatial metric $h$, the regular Legendre branch has $\pi_I=\sqrt h\,P_{IJ}y^J$. Matter normal and tangent constraint densities are $\mathcal H_m=\pi_I y^I-\sqrt h P$ and $\mathcal D_{m,i}=\pi_I\partial_i\phi^I$. At fixed momentum the derivatives of the Legendre transform give $y^I$ and $\sqrt h P_{IJ}h^{ij}\partial_j\phi^J$ as the momentum and gradient variations. Antisymmetrization of two smeared normal variations gives the usual hypersurface bracket $D[h^{ij}(N\partial_jL-L\partial_jN)]$. The mixed normal brackets with gravity cancel since the matter metric dependence is algebraic. This common closure leaves the principal symbol
$$
\mathcal P_{IK}(k)=P_{IK}g^{\mu\nu}k_\mu k_\nu-P_{IJ,KL}u^J(k)u^L(k),
\quad u^J(k)=\nabla^\mu\phi^J k_\mu
$$
sensitive to the matter kinetic Hessian. In the example $P_1=X+X^2/\Lambda^4$ plus a constant background spectator, the active speed squared is $(1+2x)/(1+6x)$ at $x=X/\Lambda^4$. The matched canonical field obeys $\dot\chi=\sqrt{2XP_{1,X}}$ and $U(\chi(t))=XP_{1,X}-P_1$. Its energy and pressure match on the whole background and its cone differs, giving $1$ versus $17/19$ at $x=1/32$. Extra equality of the full kernel on each metric-null covector over a regular two-gradient branch forces the kinetic Hessian to vanish by a rank-two projection and a quadratic polynomial identity. This is a conditional kinetic restriction. [F3 equations and brackets](../grut_f3_coupled_consistency_20261126/MODEL_AND_CONSTRAINTS.md) · [Full-cone hypothesis](../grut_f3_coupled_consistency_20261126/CONE_OBSTRUCTION.md)

## 7. Residual physical freedom in matter laws

F4 uses left-handed weights $(q,u,d,\ell,e,n,h)$ in units where $h=3$ and the physical hypercharge is $Y=y/6$. Dirac Yukawa invariance gives $(u,d,e,n)=(-q-h,-q+h,-\ell+h,-\ell-h)$. The weak mixed anomaly gives $3q+\ell=0$, which also kills the cubic trace $-6h^2(3q+\ell)$. The quotient congruence $y+2t_3+3t_2=0\bmod 6$ imposes $q=1+6k$ and thus $Y=Y_{\rm SM}+3k(B-L)$ for integer $k$. Spin global-anomaly cancellation is imported for the chosen quotient. Ordinary connected phases can be broken for each $k$ by $QQQL$ and a gauge-neutral dressing of $n^c n^c$ with powers of $H^Ti\sigma_2D_\mu H$ or its conjugate. Their independent $(B,L)$ violations span both ordinary directions. Local operator-dimension bounds depend on $k$. With the charged-vector unit fixed, the neutron-to-electron charge ratio is $3k/|1+3k|$. Changing $k$ changes observable ratios and is physically distinct. Neutral neutron evidence selects the canonical member within the stated family. [F4 proof and measurement](../grut_f4_charge_selection_20261126/DERIVATIONS.md) · [Obstruction inventory](FOUNDATIONAL_OBSTRUCTIONS.md)

## 8. Matter EFT, continuum spectrum and physical readout

The covariant action and original matching, the rotating fixed-charge reduction, the relative one-loop determinant and the continuum linear certificate are developed together in the [Candidate 1 final report](CANDIDATE_1_FINAL_REPORT.md). C1 uses the stationary Legendre transform and reports strict one-loop differences with numerical allowances. C2 uses an exact stationary-profile enclosure and inertia of the half-line radial operator. For a positive normalized charge vector $e$ with $L_+$ invertible,
$$
L_Q=L_++4\Omega^2|e\rangle\langle e|,\qquad
1+4\Omega^2\langle e,L_+^{-1}e\rangle=\frac{\bar Q'}{I}.
$$
The verified unconstrained Morse index is one and $\bar Q'<0$, so the fixed-charge correction removes that direction. Angular ordering and symmetry zero modes complete the leading continuum result. The reporter and gravitational calculations then use that prepared state and their own readout assumptions. C3's elastic cross section carries a calibrated hidden-field wave but lacks a physical source/readout chain. C4's material quadrupole is an externally forced leading contribution with a separate mass-matched comparator. [C2 certificate method](../grut_c2_continuum_stability_20261127/CERTIFICATE_METHOD.md) · [C3](../grut_c3_empirical_feasibility_20261127/REPORT.md) · [C4](../grut_c4_gravitational_accessibility_20261127/ANALYTICAL_GATE.md)
