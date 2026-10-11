# CP09: exact recovery of ideal response by clock conditioning

2026-10-11. Analytical and executable clock-reference construction. CP09 law unchanged, no CP10 or Stage3 card. This is a distinct reference construction from the previous backreacting two-qubit apparatus. It is not a claim that that apparatus now implements ideal response.

## Exact construction

For finite CP09 H|n>=E_n|n>, supply a d-dimensional clock with H_c|n>_c=-E_n|n>_c. Given a system preparation |psi0>=sum a_n|n>, construct the stationary correlated state

|Psi>=sum a_n |n>_c|n>, (H_c+H)|Psi>=0.

Define normalized clock test states |t>_c=d^(-1/2)sum exp(+iE_n t)|n>_c. Conditioning gives

_c<t|Psi>=d^(-1/2) exp(-iHt)|psi0>.

Thus the normalized conditioned system state is exactly the known evolved state. The global correlated state is stationary; prediction uses a clock reading rather than a global evolution parameter. The clock spectrum, energy matching, constraint and entangled preparation are SUPPLIED. This removes external evolution time from this conditional representation, not from the explanation of physical dynamics. For finite d a constant can shift H_c to nonnegative energies and shift the total constraint eigenvalue accordingly; this changes no correlations. It does not supply an arbitrary-N physical clock or fix absolute units.

## Detector response actually recovered

Use frozen N=2 CP09 and the disjoint minimal/maximal-role observables already defined in 8bb18f4. Prepare two clock-system histories whose conditional initial system states are exp(+/-i epsilon B)|ground>, epsilon=1e-5. Differentiate the conditional A expectations. Results:

| clock reading | conditioned impulse response | ideal Kubo response |
|---|---:|---:|
| .1 | .004121694452 | .004121694452 |
| .4 | .052588513534 | .052588513539 |
| 1 | .029091102673 | .029091102668 |

Conditioned vectors agree with direct evolution within 1.2e-15. Total constraint residual is 1.1e-15. This is exact correspondence earned by algebra plus independent matrix calculation, not a newly selected law. The construction applies generally to any finite Hermitian H: successful response recovery therefore cannot distinguish CP09 as fundamental physics.

## Crucial limitations

Clock test states at different readings are not orthogonal: squared overlap between .1 and .4 is .85736. They are legitimate separately specified rank-one binary measurement effects, but do NOT constitute exclusive sharp-time outcomes or automatically form a normalized continuous time POVM. The reported 1/d is the probability of each individually tested binary effect; it is not a probability density to integrate over all t.

The imposed stationary entangled preparation already encodes the entire supplied H. It is not obtained by running the previous autonomous recorder, and preparing each impulsed history is not the same as physically performing an intervention at a clock reading inside one history. Multi-time interventions require a specified instrument and apparatus, beyond these conditional one-time expectations.

The coordinate check t'=2t uses t=t'/2 and recovers identical correlations. This is a coordinate reparameterization, not proof of equivalence between physically different clocks or reference frames. A finite clock remains recurrent and spectrally tuned. No monotonic physical time, local reference access, light cone, continuum phase or spin2 sector has been generated.

## Verdict

The requested ideal detector response is recovered EXACTLY in an explicitly supplied stationary clock-conditioned reference model. External evolution time is replaced by a clock correlation in that representation. A dynamically prepared, localized and calibrated clock-apparatus realization remains unresolved. These achievements and limitations should not be conflated.

Run python research/reality_framework/cp09_clock_conditioned_response/check.py (NumPy/SciPy). Expected CP09_CLOCK_CONDITIONED_IDEAL_RESPONSE_PASS. Results retained. Author mathematical control, not empirical validation. The construction is conventional relational quantum-clock mathematics; no novelty claim.
