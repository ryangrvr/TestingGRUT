# GRUT-W1: a working finite operational formulation

9 November 2026. **GRUT now has an explicit, calculable effective formulation for a finite record-action system and a controlled predictor of the adopted massless continuum response.** The constructive correction is to retain the information that future readouts can use: parity coherence in the record model, and the actual frequency distribution in the continuum model. The resulting predictions follow from specified equations and parameters rather than a fit to the target outcomes.

This is a working mathematical/effective candidate within quantum mechanics. The physical inputs, apparatus correspondence and probability law are stated below. W1 does not supply a new fundamental dynamical or state-selection law, an empirical confirmation, or a completed lattice-to-continuum detector dictionary.

## 1. The theory's primitive and prediction law

Specify one operational parent

    G=(A,rho0,H,I_access).

A is the physical observable algebra, rho0 a positive trace-one preparation, H the Hamiltonian, and I_access the actual site algebras, permitted couplings and instruments. These are physical inputs. Mathematical tensor factors alone do not establish independently accessible pieces of an apparatus. The finite theory adopts the Born rule and unitary dynamics, with hbar=1:

    U_t(rho)=exp(-iHt)rho exp(iHt),
    J_x(rho)=sum_b K_xb rho K_xb^dagger,
    sum_xb K_xb^dagger K_xb=I,
    p(x1,...,xn)=Tr[J_xn U_dt_n ... J_x1 U_t1(rho0)].

This law fixes both outcome probabilities and backaction. A selected branch is divided by its positive trace. Ignoring an outcome retains the nonselective channel; omitting the measurement is another physical protocol. Under its closed unitary evolution, a time-independent parent conserves its total energy distribution. A driven reduced apparatus may exchange energy with controls that its model must identify.

The word responsive has an operational meaning: changing an allowed intervention changes the probability of a declared later readout. A record is a correlation that remains reliable under a stated class of disturbances; an actuator uses that correlation through its specified coupling. These definitions do not introduce a consciousness variable, a purpose, or an extra physical force.

## 2. The effective-description principle that makes this formulation work

Fix the preparations and future experiments before compressing the state. Retain enough expectations to determine every backward-propagated readout effect, or bound the discarded contribution. This is the operational sufficiency condition of the existing candidate, made explicit here.

In finite dimension, suppose retained variables are Tr(rho X_j), including X_0=I, and all states are admitted. Every relevant probability has a linear reconstruction from these variables precisely when its effect belongs to span_R{X_j}. If an effect lies outside that span, its orthogonal component supplies a traceless Hermitian Y. For sufficiently small epsilon, rho+/-=I/d+/-epsilon Y are valid states with the same retained variables but different predictions. Thus no reconstruction from those variables alone can be universally correct for that experiment.

For a restricted preparation family the condition can be weaker, but the family and error criterion must be declared. This explains the two corrections in W1. Population-only parity memory discards transverse components that later pulses can convert to action. The R1 three-source first-moment approximation discards the distribution of field frequencies that controls later response. Neither defect can be repaired merely by declaring the retained variables sufficient.

## 3. A complete finite record-action model

The four work qubits are signal S, memory carriers M1,M2, and actuator A. Prepare |+>_S |Phi+>_M1,M2 |0>_A and write P=Z_M1 Z_M2. The work sequence is: write S into parity using CNOT(S,M2); apply the historical within-parity hold; rotate M2 by theta=1 about X; apply the declared parity channel; apply exp(-i beta P/2); rotate M2 by theta2 about X; and copy parity to A using CNOT(M1,A) and CNOT(M2,A).

The coherence parameter eta in [0,1] defines

    Lambda_eta(rho)=(1+eta)rho/2+(1-eta)P rho P/2.

It preserves the two parity sectors and multiplies coherence between them by eta. Eta is distinct from the old hold-angle symbol lambda. Eta=1 is the identity channel; eta=0 is nonselective parity dephasing. Measuring the two memory bits individually is a different instrument. The usual projector measurement rule is the adopted quantum rule, not a new GRUT postulate; see [Luders' original rule, translated](https://arxiv.org/html/quant-ph/0403007v2).

The memory basis (m1,m2) can be rewritten as (g,p)=(m1,m1 xor m2). The Bell blank then becomes |+>_g|0>_p. The writer acts on S,p; the disturbances and phase act on p; the actuator couples p,A; the historical hold acts only on g. Therefore the clock-conditioned work predictions for Z_S,P,Z_A have an exact three-logical-qubit description on S,p,A. Measurements of the discarded complementary memory phase and arbitrary clock coherences require the larger description.

The logical Bloch vector after writing is (0,0,s), where s is the signal sign. The first X rotation gives (0,-s sin(theta=1),s cos(theta=1)). The parity channel attenuates the transverse components by eta; the parity phase rotates them by beta; the second X rotation converts the surviving transverse component into parity. Consequently

    p_agree(theta=1,theta2,beta,eta)
      = [1+cos(theta=1)cos(theta2)
           -eta sin(theta=1)sin(theta2)cos(beta)]/2.

This is the working finite theory's quantitative prediction. A disabled disturbance sets the corresponding theta to zero. At the historical angles theta=1=theta2=pi/3, beta=2pi/3,

    p_agree=5/8+3eta/16.

Thus eta=1 predicts 13/16 and eta=0 predicts 5/8. Independent population flips predict the eta=0 value for every beta. That is a restricted comparator, not every classical explanation. The general formula makes predictions away from the familiar example; it is not an interpolation fitted to the two displayed endpoints.

The finite pulse model in FINITE_MODEL.txt specifies positive durations and exact bounded Hamiltonians for every operation. Its illustrative 660 ns schedule has no instantaneous-pulse approximation. The duration is a model specification, not an achieved hardware speed or fidelity. A finite random I/P operation realizes eta=0 without an infinite dissipative rate; a finite-rate Markov dephasing semigroup is a different implementation and reaches eta=0 only asymptotically. Unknown detuning and other noncommuting pulse errors require an explicit remainder.

## 4. The same finite work model has an explicit closed parent

Introduce an environment qubit E prepared in sqrt((1+eta)/2)|0>+sqrt((1-eta)/2)|1>, and couple it by controlled-P. Tracing E gives exactly Lambda_eta. Add a nine-position propagation clock and two conserved selector qubits D1,D2. The new parent's dimension is 9*2^7=1152.

Let U_k be the eight specified work links, including controlled-P, and define

    H=sum_(k=0)^7 sqrt((k+1)(8-k))/2
      [|k+1><k| tensor U_k + |k><k+1| tensor U_k^dagger].

Every link is present throughout the autonomous evolution. With V0=I, V_(k+1)=U_k V_k and T=sum_k |k><k| tensor V_k, direct multiplication gives

    H=T(J_x tensor I)T^dagger,

where J_x is spin 4. Starting at tick 0, every supported clock tick has the ordered work state V_k rho_work V_k^dagger. The clock distribution is binomial with 8 trials and p=sin^2(t/2), and transfers to tick 8 at t=pi. This uses the standard engineered spin-chain construction ([Christandl et al.](https://arxiv.org/abs/quant-ph/0309131)); its work-dependent block identity is verified directly here.

This is a new dilation. The historical 512-dimensional parent retains its own specification and verdict. A real four-qubit scheduled experiment would test the work protocol, not automatically either autonomous clock Hamiltonian. The new parent's access table must physically permit its clock-work couplings; their engineering is not proved by the algebra. Its finite clock can reverse and undo the environment coupling, so irreversible noise or a cosmological arrow of time is not derived.

## 5. Calibration turns the formula into a physical prediction

The zero/one-disturbance agreement values do not determine beta, eta or the sign of an angle. Finite countermodels reproduce the same calibration values (1,3/4,3/4) while giving two-pulse predictions from 1/4 to 1. Phase-sensitive calibration is therefore necessary before the target prediction is opened.

Use independent quadratures for the X rotation angles; conditional X/Y probes of the parity phase for beta; cross-parity coherence for eta; and within-parity Bell coherence to qualify the instrument. The last control matters: a local Z_M2 channel gives the same agreement curve and immediate populations here but damages a Bell state that the intended parity channel preserves. Agreement statistics alone do not identify the physical channel.

For independently established angle/eta radii, one conservative global event bound is

    |Delta p| <= (delta_theta=1+delta_theta2+delta_beta+delta_eta)/2
                 + D_prep + delta_readout + sum(delta_channel).

The angle radii are in radians. D_prep is half the trace-norm preparation error, delta_readout an effect-operator error, and delta_channel a justified half-diamond channel bound on the relevant input class. Tighter interval-aware bounds and finite pulse remainders are derived and checked in FINITE_MODEL.txt. These budgets are conditional on independently justified apparatus assumptions; fitting the target outcome cannot supply them.

In the registered theoretical uncertainty example, the actual deterministic corner discrepancy is 0.025968264775, below its derived allowance 0.0581217286339. This is a check of the inequality, not a measured calibration. The separate 960,000-attempt structured-calibration example assumes additional certificates and is not a demonstrated budget-qualified EX1 calibration.

PHYSICAL_ACCEPTANCE.txt retains the M1 requirement: freeze a named apparatus, a supported calibration confidence set, raw prediction bands and an error ledger before held-out measurement. Retain erasures, count each error once, and require both registered accuracy and joint prediction-budget consistency. The adapter is still incomplete; no physical data or hardware job is included in W1.

## 6. A corrected continuum response predictor

The continuum instance is the adopted massless Schwinger scalar response at L=5, g=1, with Dirichlet modes e_n(z)=sqrt(2/L)sin(n pi z/L), omega_n=sqrt((n pi/L)^2+1/pi), and Phi=phi/sqrt(pi). The N22 native cell-source profiles and their apparatus boundary constants are unchanged from D1/R1. This is a second specified instance of the operational framework; the shared framework does not prove the lattice and continuum descriptions are physically equivalent.

For centered vacuum probes, define b_in=int f_i(z)e_n(z)dz and a_n=b_n/sqrt(2 pi omega_n). Then

    M(t)=sum_(n>=1) a_n a_n^T exp(-i omega_n t).

W1 retains the actual frequencies and positive spectral weights:

    M_K(t)=sum_(n<=K) a_n a_n^T exp(-i omega_n t),
    C_K=M_K(0),  R_K(t)=C_K^(-1/2)M_K(t)C_K^(-1/2).

K=1024 was fixed before the new comparisons; 16,64,256 were diagnostic levels only. This changes the approximation's retained dynamics, not the source dictionary. The response is normalized at t=0 and contractive at every time in the exact mathematical model.

For the omitted covariance Delta_K>=0, ||M(t)-M_K(t)||<=||Delta_K|| uniformly in t. If each zero-extended step source has total variation B_i, then

    ||Delta_K|| <= E_K=sum_i B_i^2 L^2/(2 pi^4 K^2).

With C_K eigenvalues bounded below by m and above by q,

    ||R_infinity(t)-R_K(t)|| <= E_K/m + E_K*q/m^2.

This uses variance tails and an inverse-square-root bound, not a finite H-squared moment or a short-time series. For the primary three native sine probes, a priori Poincare, edge and variation estimates give exact rational bounds m>=1323/163328, q<=1/16, E_K<=7/(9K^2). Therefore the exact 1024-mode predictor satisfies

    sup over all real t ||R_infinity(t)-R_1024(t)||
      <=3678389/4608866304
      approximately 0.0007981114567822361 <0.001.

This is the new mathematical result: a specified continuum-response approximation with a uniform modal-truncation bound below 0.1%. It fixes the R1 three-source predictor's loss of spectral information. It does not repair the separate R1 lattice-to-continuum mismatch, validate selective detector histories, or identify the detector map physically.

The theorem concerns exact mathematical coefficients and phases. Floating frequency errors accumulate phase with time, so it is not an arbitrary-time machine-precision guarantee. The executable was checked at the six frozen times 0,.17,.61,1.09,2.03,3.17 against fixed references 16384/32768. Beta_left was retained as a separate scalar witness outside the primary three-source theorem.

## 7. Verification results and their scope

Source/protocol freeze: [6461c320b87bea70c97916b625e2b14511a3a1a5](https://github.com/ryangrvr/TestingGRUT/commit/6461c320b87bea70c97916b625e2b14511a3a1a5). The 500-point finite-model regression was inspected before this freeze; the autonomous-parent evolution and W1 continuum response comparisons were not. All three registered calculations completed in 3.656 seconds within the fixed resource caps. No angles, times, thresholds, source weights or primary mode count were adjusted after evaluation.

| Test | Result | What it establishes |
|---|---:|---|
| 500 four-qubit density-matrix cases | Maximum formula error 3.33e-16 | Parametric finite-model identity and implementation |
| 64 autonomous-parent records | Maximum conditional TV 3.44e-15 | New closed dilation matches its ordered work predictions |
| Primary 1024-mode predictor, six new times | Maximum normalized difference 7.80926299243e-07 | Finite-time floating implementation check |
| Exact primary uniform truncation theorem | 0.0007981114567822361 | Mathematical all-time truncation accuracy under the specified continuum model |
| Reference 16384-to32768 change | 1.97242833215e-09 | Numerical reference stability at the registered times |
| Independent frozen R1 implementation | 8.67e-16 | Source and response implementation agreement |

All 10 spectral numerical-quality gates passed. Three separate audits passed, checking the 500-point formula/budget arithmetic, the 64 parent records and clock statistics, and all spectral levels/times/tail decisions. Each verified the relevant frozen source provenance; there are 20 hashed source files. The audits do not claim a second full density/evolution/source computation: the parent conditional TVs are reaggregated from archived scalars, and the spectral audit reconstructs matrices and comparisons from archived mode results. Source audits and mathematical proofs remain distinguishable from numerical checks.

The diagnostic spectral levels were:

| Retained modes | Maximum measured normalized error | Floating all-time tail estimate |
|---|---:|---:|
| 16 | 0.00176158871796 | 0.48480928101 |
| 64 | 0.000199322781846 | 0.0302426159219 |
| 256 | 1.28462305982e-05 | 0.0018897462459 |
| 1024 | 7.80926299243e-07 | 0.000118107157924 |

The primary choice remains 1024. Its exact a priori certificate is looser than its data-dependent floating tail estimate and much looser than the six sampled errors; these are three distinct quantities.

## 8. The boundary of the working theory

For finite-dimensional Hermitian H,G, Tr[G,[H,G]]=0. A positive continuum current cI therefore cannot be approximated with relative error below 100% uniformly over every finite-dimensional state. The dictionary must name a physically justified state/source domain, include any claimed reached detector branches, and carry explicit error bounds. DICTIONARY_BOUNDARY.txt gives that proof and the state/preparation/pulse requirements.

The previous guarded |kappa|<=0.7057 condition remains a restricted current result on archived R1 matrices. It is not a rigorously certified physical pulse window, a full continuum dictionary, or evidence for the failed selected-state current claims. N22 retains its completed finite massless pass. B2's N20 component failure, D1/R1 mismatches, withheld massive transfer and P1 interacting-parent obligations remain separate.

What can now be called working is the finite operational theory conditional on its explicit inputs, its exact record-action formula, its closed mathematical dilation, and the adopted-continuum response approximation with the stated theorem. What remains open is physical realization with calibrated held-out data, a justified full lattice/continuum instrument correspondence, and any additional empirical content beyond the adopted quantum framework. The current candidate supplies no law selecting the universe's algebra, state, couplings or constants, and no derived gravitational, consciousness or cosmological completion.

The immediate physical use is concrete: measure the specified apparatus parameters independently, compute the frozen raw-outcome predictions and budget, and confront the held-out record-action data. A failure must change or reject the stated model at its declared scope. The scope does not need to expand to consciousness or cosmology for that test to be meaningful.

## Reproduction and source map

The archive contains all frozen study sources, imported source files, raw outputs, resource records and independent audit programs/results. SHA256.json covers every payload except itself. The production runner requires a clean Git checkout of the freeze above, Python 3.12, numpy 2.3.5 and scipy 1.17.0 on Windows, and a new external output directory. The copied sources support inspection; their folder is not itself a Git checkout.

CORE_THEORY.txt is the earlier source audit of the existing integrated candidate; its initial no-production statement refers to that audit, not this completed W1 suite. FINITE_MODEL.txt and AUTONOMOUS_DERIVATION.md define the new working instances. SPECTRAL_DERIVATION.txt proves the new response bound. PHYSICAL_ACCEPTANCE.txt and device_adapter.template.json retain the unfinished physical qualification explicitly.

The adopted core and original Construction 1 are pinned to the [integrated candidate](https://github.com/ryangrvr/TestingGRUT/blob/4300c1836c4967a268fe3b028fa190e0f9886153/research/reality_framework/grut_integrated_candidate_20261106/GRUT_Integrated_Theory.md) and [Construction 1 specification](https://github.com/ryangrvr/TestingGRUT/blob/4300c1836c4967a268fe3b028fa190e0f9886153/research/reality_framework/grut_construction_1_20261106/SPECIFICATION.md). No historical report is rewritten, and no hardware job, massive test, merge or release occurred.

Offline audit commands and their limits are in AUDIT_REPRODUCTION.txt.
