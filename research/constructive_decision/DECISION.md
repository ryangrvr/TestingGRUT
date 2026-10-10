# Constructive-theory decision: a calibrated effective model, no new law claim

2026-10-10. Authorized development branch. Starting revision:
`b1746edec0021b35925f7b8168d6eccb8de44e9e`.

## Decision and scientific limit

**Choose Approach B for the next executable milestone.** Use an independently
calibrated spin-dependent force to forecast qubit coherence and oscillator
energy in held-out control protocols and independently prepared temperatures.
Implement the transfer, including uncertainty, rather than commissioning
another reconstruction countermodel.

**This is a conventional effective model, not a distinctive GRUT theory.** Its
predictions follow from ordinary quantum mechanics once the Hamiltonian, state
and controls below are supplied. A future unmeasured outcome can be predicted
without being a new universal law. This packet establishes the calculation and
the prospective test; it contains no apparatus calibration or experimental
outcome. It does not satisfy the stronger requirement of excluding a model that
the same full standard-physics baseline permits by an additional GRUT principle.

No independently supported new interaction principle was identified by this
synthesis. Inventing a shared constant or coupling between sectors to obtain a
desired forecast would supply the missing hypothesis, not derive it. Approach A
therefore has no proposed law here. **Card 2 stays unopened; two slots remain.**
The killed Card 1 is neither repaired nor repackaged. External evaluator/CR-5
clearance, independent scientific review and canonical integration remain
separate and unestablished by this work.

## What the accepted record licenses

| Earned result | Consequence for this decision |
|---|---|
| R1 identifies response modulo a separately certified common carrier | Keep carrier and readout calibration explicit; no fundamental-law credit |
| Level-0/S5 does not select a temporal generator | Supply a measured physical generator for an effective theory |
| Quantum-lift results do not derive the physical algebra | Supply oscillator and qubit algebras; price that premise |
| Recon02–04 preserve real freedom under different environmental controls | Specify the complete applied control history; do not demand uniqueness across different experiments |
| Interaction-Origin 01 leaves higher interactions invisible to lower-sector channels | Do not infer an interaction from observations insensitive to it; calibrate the force itself |
| GR2 leaves spin, reach, coupling class and causal structure supplied | No lab-to-gravity transfer without an independently established bridge |
| P3/P4 shows identical two-point kernels need not determine non-Gaussian influence | Declare the full Gaussian state/model assumption; do not infer it from a covariance fit |

These are limitations on reconstruction from the stated information, not
evidence that a universal selector exists. They allow a finite, conditional
effective model with measured inputs. They do not select it as nature's ontology.

## The supplied model

Take one qubit with Pauli operator `Z` and `M` independently addressable harmonic
modes. The Hilbert space is

\[
\mathcal H=\mathbb C^2\otimes\bigotimes_{j=1}^M\ell^2(\mathbb N_0),\qquad
[a_j,a_k^\dagger]=\delta_{jk},\quad Z^2=I.
\]

`M` is finite and frozen before the holdout. Physical mode frequencies
`omega_j > 0` are measured in radians per second. In the frame rotating with the
known qubit and oscillator free Hamiltonians, apply the specified force

\[
H_I(t)=\hbar Z\sum_{j=1}^M
\left[f_j(t)e^{i\Omega_jt}a_j^\dagger+
f_j(t)^*e^{-i\Omega_jt}a_j\right].                 \tag{1}
\]

`Omega_j` is the signed force detuning, not the physical oscillator frequency.
`f_j` has units `s^-1`; it is a bounded, piecewise constant calibrated complex
force envelope. The command protocol selects calibrated envelopes and their
signs. Finite switching errors must be enclosed by the force/timing error
budget, not dismissed as free pulses. All times use a supplied laboratory clock.
In a laboratory oscillator frame the force coefficient is
`f_j(t) exp(i(Omega_j-omega_j)t)`; the positive physical free Hamiltonian remains
`sum hbar omega_j a_j^dagger a_j`. Completing its square shows a finite lower
bound at every finite applied force. The explicit displacement propagator below
defines the evolution of this finite-mode model on the oscillator Fock space;
there is no Markov limit or divergent ancilla stream.

At time zero supply a qubit superposition and the product state

\[
\rho_E=\bigotimes_j\rho(n_j),\qquad
\rho(n)=\frac1{n+1}\sum_{r=0}^{\infty}
\left(\frac n{n+1}\right)^r|r\rangle\langle r|,
\quad 0\le n<\infty.                              \tag{2}
\]

The qubit/environment preparation is initially a product. Each mode state is
centered, isotropic Gaussian with measured occupation `n_j`. If it is genuinely
thermal, `n_j=(exp(hbar omega_j/(k_B T_j))-1)^-1`. **Never replace `omega_j` by
the drive detuning in this temperature relation.** Squeezing, intermode
correlations, non-Gaussian populations, heating during the sequence, force
nonlinearity and unmodelled modes require a different predeclared model or a
quantitative error certificate. Occupation alone does not certify Gaussianity.

For trapped-ion implementations, the rotating-wave and Lamb-Dicke reductions
are additional effective-model premises with device-specific error obligations.
No available paper certifies those errors for an unspecified apparatus here.

## Exact parameter transfer

Define dimensionless displacements

\[
\alpha_j(t)=-i\int_0^t f_j(s)e^{i\Omega_js}\,ds.    \tag{3}
\]

The commutator of (1) at two times is a scalar times `Z^2=I`. Thus the Magnus
series terminates after its scalar second term. For either eigenvalue `z=+1,-1`,

\[
U_z(t)=e^{i\Phi(t)}\prod_jD_j(z\alpha_j(t)),\qquad
D(\alpha)=e^{\alpha a^\dagger-\alpha^*a}.          \tag{4}
\]

The phase `Phi` is common to the two branches and cancels in their relative
coherence. Taking the trace of `U_+ rho_E U_-^dagger` gives the product of
characteristic functions at `2 alpha_j`. The state (2) has
`Tr[rho(n) D(beta)]=exp[-(n+1/2)|beta|^2]`. Therefore

\[
V(t)=\exp[-\Gamma(t)],\qquad
\Gamma(t)=2\sum_j(2n_j+1)|\alpha_j(t)|^2.          \tag{5}
\]

This is exact in (1)–(2), including finite memory and revivals. A monotonic
positive-rate master equation has not been assumed. In a branch, displacement
changes `a_j` by `z alpha_j`; the zero initial mean gives

\[
\Delta E_j=\hbar\omega_j|\alpha_j|^2,
\qquad
-\log V=2\sum_j(2n_j+1)\frac{\Delta E_j}{\hbar\omega_j}. \tag{6}
\]

`Delta E_j` is the change in the physical oscillator energy. Compare free-mode
energies with the force switched off at the initial and final measurements;
switching work and its errors are part of the external control, not omitted
resources. The externally driven force supplies the mode energy increase.
This is not a calculation of total apparatus dissipation. It is not energy generated by decoherence or
an equilibrium heat current. Equations (5)–(6) jointly forecast coherence and
mode energy from the independently calibrated forces, frequencies and states.
The right-hand side of (6) can alternatively be tested using independent
energy measurements, with their own errors. No coherence fit is needed.

For a single constant real force `g` and nonzero `Omega`, a duration `t` gives

\[
\Gamma_{\rm constant}=4(g/\Omega)^2(2n+1)
(1-\cos\Omega t).
\]

Changing the force sign halfway through that duration gives

\[
\Gamma_{\rm reversal}=32(g/\Omega)^2(2n+1)
\sin^4(\Omega t/4).
\]

At a complete detuning period the first sequence has `V=1`; the sign-reversed
sequence has `V=exp[-32(g/Omega)^2(2n+1)]`. This is a useful held-out contrast,
not a GRUT discovery. At `Omega=0`, evaluate the integral (3) directly; there is
no division by zero. Temperature transfer changes the independently measured
`n`, retaining the independently calibrated generator only if its stability has
also been certified.

## Calibration before any held-out coherence measurement

| Quantity/premise | Independent source | Holdout restriction |
|---|---|---|
| Physical frequency and force detuning | Mode spectroscopy and clock/drive frequency calibration | Fixed values and intervals; detuning not used as a thermal energy gap |
| Complex force for each allowed drive setting | Population-prepared `Z=+1,-1` runs and calibrated motional quadrature response | No superposition coherence data used to estimate the force |
| `n_j` and state class | Independent thermometry plus preparation validation | No fitting occupation, squeezing or model order to holdout coherence |
| Readout probabilities for `+X,-X` | Separate spin readout calibration in the aligned Ramsey frame | Fixed before holdout; phase drift enclosed in model error |
| Applied waveform and timing | Recorded commands plus measured force and clock bounds | Whole-segment error bounds, including switches |
| Residual model discrepancy | Independent control experiments or analytical apparatus certificate | Bound declared before holdout; no fitted background decay |

For centered states, the population-conditioned displacement difference is
`<a>_+ - <a>_- = 2 alpha`. With `X=a+a^dagger`, `P=i(a^dagger-a)`, a separate
quadrature measurement gives
`alpha=((<X>_+-<X>_-)+i(<P>_+-<P>_-))/4`. Short pulse slopes calibrate the force;
the subsequent holdout uses different durations/sign sequences and fresh
preparations. The readout may involve further physical operations; their
calibration and resource costs remain inputs.

Freeze model inventory, calibration, protocol, uncertainty bounds and analysis
before collecting the held-out outcomes. The CLI accepts these files and
produces hashes, predictions and binomial acceptance regions. It has **no**
argument for fitting observed coherence. A separately recorded holdout is
compared against the already frozen prediction. A failure rejects the joint
model/calibration/error certificate; it does not by itself refute quantum
mechanics or establish GRUT. Any revision requires fresh holdout data.

## Whole-domain uncertainty and falsification

Suppose nominal segment endpoints are `a_l,b_l`, duration `d_l`, duration
uncertainty `e_l<d_l`, force `f_l` with disk radius `r_l`, and detuning uncertainty
`Delta Omega`. Let `E_a=sum_{k<l} e_k`, `E_b=E_a+e_l`. Actual endpoints can move
by at most these quantities. For every admitted force, detuning and duration,

\[
|\alpha-\alpha_0|\le R=\sum_l\left[
r_l(d_l+e_l)+|f_l|\left(E_a+E_b+
\min\{2d_l,\Delta\Omega(b_l^2-a_l^2)/2\}\right)\right]. \tag{7}
\]

Proof: split the integral difference into force error on the actual interval,
endpoint motion at fixed detuning, and detuning change on the nominal interval.
Use unit modulus, interval symmetric-difference length at most `E_a+E_b`, and
`|exp(i delta s)-1| <= min(2,|delta|s)`. The displayed minimum bounds the whole
integral even if it is loose. This is a whole-domain analytic enclosure, not
agreement at sampled parameter points.

Set `A_- = max(0,|alpha_0|-R)`, `A_+ = |alpha_0|+R`,
`n_-=max(0,n_0-r_n)`, `n_+=n_0+r_n`. Then

\[
2\sum_j(2n_{j,-}+1)A_{j,-}^2\le\Gamma\le
2\sum_j(2n_{j,+}+1)A_{j,+}^2.                     \tag{8}
\]

Exponentiation reverses these bounds for visibility. Positive physical
frequency intervals similarly enclose (6)'s per-mode energy prediction.
The software inflates elementary floating-point results by a disclosed
roundoff allowance; this is **not** a formally certified interval-arithmetic
implementation. Physically claimed calibration coverage and model error are
also not certified merely by entering them into JSON.

Let `p_+,p_-` be the calibrated probabilities for bright detection on `+X,-X`.
The expected raw bright probability in the phase-aligned holdout is

\[
p(V)=\tfrac12[p_+(1+V)+p_-(1-V)].                 \tag{9}
\]

Its extrema over the readout and visibility intervals occur at their eight
corners. Widen by the independently justified absolute probability model
error, then clip to `[0,1]`. For `m` frozen settings, `N_s` independent shots
per setting and total sampling failure probability `alpha_s`, use

\[
\eta_s=\sqrt{\log(2m/\alpha_s)/(2N_s)}.            \tag{10}
\]

A setting fails if the observed bright fraction lies outside its predicted
interval widened by `eta_s`. Hoeffding plus a union bound gives joint sampling
coverage at least `1-alpha_s`. Adding the claimed simultaneous calibration
failure probability `alpha_c` gives at most `alpha_s+alpha_c` false rejection
under the model and certificates. Drift or dependent shots invalidates this
claim unless separately bounded. Energy measurement is a second holdout, not
covered by this spin-shot acceptance rule.

## Conventional comparison and novelty audit

The strongest conventional comparator is **the identical model (1)–(2)**.
It predicts all of (3)–(10). The assembled model is already realized in
spin-dependent-force and oscillator characteristic-function frameworks.
Haljan et al. [S1] supply the controlled-displacement construction; Fluhmann
and Home [S2] give direct characteristic-function readout and independent
calibration comparisons. Sawyer et al. [S3] experimentally used spin dephasing
to probe motional temperature and state distributions. Even the intended
quantum/thermal connection therefore has an explicit conventional precedent.

Consequently a successful holdout would validate a calibrated apparatus model
and its parameter transfer. It would not remove a point from the full
standard-physics joint image by a new principle. A broader conventional model
can use squeezing, non-Gaussian states, extra modes or additional interactions;
the premises and independent certificates, not a new law, exclude these here.
Recon04 and Interaction-Origin 01 already explain why omitted dynamics cannot
be reconstructed from insensitive observations. This packet does not add a new
countermodel campaign.

GRUT's potential contribution here is disciplined assembly, premise tracking,
common calibration and prospective validation. No numerical agreement, new
notation or unification of the bookkeeping is claimed as new physics.

## Assumption and information inventory

**POSTULATED/SUPPLIED:** quantum algebra and Born measurement rule; qubit/mode
partition; finite mode inventory; positive physical frequencies; fixed spin
axis and linear force operator; drive envelopes and detunings; laboratory
clock; Gaussian product state and initial system/environment independence;
preparation/readout; waveform, drift and omitted-dynamics certificates; energy
readout; resource ceilings. The implemented schema carries `M` physical
frequencies, `M` detunings, `M` occupations, a finite library of complex force
amplitudes, their uncertainty bounds, and all protocol durations/signs. These
are measured inputs, not derived constants.

**DERIVED within those premises:** displacement dynamics, coherence,
physical-mode energy, relation (6), and the analytic uncertainty enclosure.

**NOT DERIVED:** persistent subsystem existence; a uniquely selected readout or
interface; quantum mechanics; gravitational couplings; a universal relation
between lab couplings and gravity; initial thermalization; universal observer
identity; a new GRUT restriction. This is not an Appendix-B law submission and
receives no candidate score or purported `L_0` price. The explicit premise
inventory must not be mistaken for eliminating those information costs.

In particular the shared parameter for a lab force is **not** identified with
gravitational coupling or a cosmological memory time. Existing GR2, USL tabletop
and Gamma_T gates do not license that transfer. No gravity forecast is emitted.

## Completion status

The selected effective model, mathematical transfer, calibration contract and
forecast implementation are concrete. **32 author software checks passed**;
the largest coherence/occupation discrepancy against independent finite-Fock
Hamiltonian evolution was about `7.8e-15` at cutoffs 24 and 40. Reproduction tests compare the transfer
against independent finite-Fock Hamiltonian evolution and check the uncertainty
and input contract. They are author-level software controls, not new theory
evidence. `results.json` reports them.

No real calibration bundle or held-out measurements are present. Thus the
experimental stage is **AWAITING INDEPENDENT CALIBRATION**, not validated.
The implementation can use a supplied independently calibrated bundle without
refitting; its receipt is not needed to publish this constructive decision.

**Scientific verdict:** no new fundamental or distinctive effective GRUT law
established; no candidate selected for Card 2. Approach B is selected as a
testable conventional effective-model milestone. A future distinctive physical
commitment needs its own explicit authorization and pricing, rather than being
hidden in a forecast correction. This completes the requested synthesis and
constructive decision; it does not answer the theory-of-reality question.

## Sources and source boundary

- [S1] Haljan et al., *Spin-dependent forces on trapped ions for phase-stable
  quantum gates and motional Schrodinger-cat states*,
  https://arxiv.org/abs/quant-ph/0411068 (especially controlled displacements
  and the effective interaction Hamiltonian).
- [S2] Fluhmann and Home, *Direct characteristic-function tomography of quantum
  states of the trapped-ion motional oscillator*,
  https://arxiv.org/abs/1907.06478 (characteristic-function definition, readout
  and comparison to independent calibrations).
- [S3] Sawyer, Britton and Bollinger, *Spin dephasing as a probe of mode
  temperature, motional state distributions, and heating rates in a
  two-dimensional ion crystal*, official author/institution record:
  https://www.nist.gov/publications/spin-dephasing-probe-mode-temperature-motional-state-distributions-and-heating-rates
  (the experimental precedent, not a certificate for an apparatus here).
- Repository boundaries: `GRUT_WORKING_THEORY_SYNTHESIS_01.md`,
  `GRUT_PROGRAM_STATE_SYNTHESIS_01.md`, `GRUT_MODEL_FRAMEWORK.md`,
  `SIGNATURE_AUDIT.md`, `GRUT_PREDICTION_GATE_GAMMA_T.md`,
  `program/gates/Q1_USL_TABLETOP_GATE.md`, `P3_P4_OWNER_RULING_01.md`,
  `GR2_CAMPAIGN_SYNTHESIS_01.md`, and the existing `research/` packets.
  The manifest pins the local source bytes; their historical governance status
  is not promoted by citing them.
