# R1-01: a bounded working sector, and a sharper physical obstruction

Completed 9 November 2026. **R1 identifies a useful finite lattice predictor and a weak unitary-kick current regime. The stronger continuum detector correspondence remains unresolved.** All 15 numerical quality gates and the independent audit pass. Five registered points completed in 247.469 seconds, with peak committed memory 1.138 GiB under the 4 GiB cap.

The protocol and code were frozen at [18fc45c13b981992e8a1c6845f9b4182e6b4d5a2](https://github.com/ryangrvr/TestingGRUT/commit/18fc45c13b981992e8a1c6845f9b4182e6b4d5a2) before production. N=14,18,20,22 and refined N22 used the hashed historical ground states. No new larger lattice, detector refit, tolerance change or physical retry occurred. A postproduction algebraic supplement is explicitly separate from the registered decisions.

## What works

The three-source lattice predictor was built only from static covariance C and first energy moment K. At new held-out times 0,.29,.83,1.57, its maximum whitened response error is 0.000292773845542, or 0.029277%. Including the frozen 1e-4 guard gives 0.039277%, within the 1% target. This is a completed finite response-prediction task on the registered source span and times; it does not validate nonlinear detector histories or a real apparatus.

The zero-delay one/two-instrument moments also pass: the maximum error 0.000916443190124, plus 1e-4 reserve, is below .002. The worst is the beta_left second moment in the coincident `--` branch. The static first energy-moment matrix passes its1% criterion. The prospectively tested unitary kick strengths kappa=.5 and 1 pass both the disturbance-from-vacuum and absolute continuum-current criteria.

These are distinct successful tests. Combining them does not establish that a finite-strength nonlinear trajectory is predicted to 1%, because that combination was not tested.

## What the bounded diagnostic separates

For fixed source operators G_i, let Xi_i=(G_i-<G_i>)|Omega>. Compute

    C=Xi^dagger Xi,
    K=Xi^dagger(H-E0)Xi,
    M(t)=Xi^dagger exp[-i(H-E0)t]Xi.

M_ii(t)/C_ii removes scalar gain and offset. Generalized spectra of (K,C) and (Hermitian M(t),C) remove invertible source mixing within the registered three-sine span. Both exact signed amplitude/shape decompositions are archived. They are algebraic decompositions, not unique causal percentages.

At N22 the largest normalized-response discrepancy is 0.0490637839606, for sine3 at t=1.57. The centroid discrepancy is 0.0343830815787; the cosine-spectrum discrepancy is 0.0255944718207. All fail their .005 criteria. An overall source gain/offset, or invertible constant mixing of these three sources, cannot repair these discrepancies. Changes outside the source span, prepared state and boundary correspondence remain possible, so this is not a unique attribution to H alone.

The exact lattice source-current rule supplies the zero-evolution diagnostic:

    S_gh(psi)=<[G_g,[H,G_h]]>_psi
             =-sum_j v_gj v_hj <H_hop,j>_psi.

The matched continuum rule is S_c=(v v^T)/(pi a), independent of the reached state. The finite lattice instead changes its bond expectations after a kick, with the independently checked identity h_j(kappa)=cos(kappa w_j)h_j(0) for this real initial state.

At zero free-evolution time, the worst coincident selective-branch current mismatch is 0.0574416128539 (about 5.74%), and the twice-unitarily-kicked state differs by 0.0704938396835 (about 7.05%). These discrepancies already exist without accumulated propagation. The zero-delay moment pass does not establish the stronger source-current dictionary across reached states. Likewise it can coexist with the 4.26% relative covariance miss, because the moment and covariance criteria have different absolute/relative normalizations.

The earlier D1 delayed-moment error .00397959128 remains a miss. R1's passing coincident moments and failed reached-state current rules show why assigning all of D1's error to either static calibration or free evolution would be unjustified.

## Registered N22 decisions

Each row uses the maximum across standard and refined N22. Values below already include the 1e-4 guard/reserve. `pass` in assessment.json means numerical eligibility only, never an aggregate theory pass.

| Test | Guarded metric | Limit | Decision |
|---|---:|---:|---|
| Static covariance, relative primary matrix | 0.04265537577 | 0.01 | miss |
| Static first energy moment, relative primary matrix | 0.007025817639 | 0.01 | pass |
| Normalized response, all four probes | 0.04916378396 | 0.005 | miss |
| Generalized centroid spectrum, relative | 0.03448308158 | 0.005 | miss |
| Generalized cosine spectrum, absolute | 0.02569447182 | 0.005 | miss |
| Lattice three-source predictor, whitened norm | 0.0003927738455 | 0.01 | pass |
| Continuum three-source predictor, whitened norm | 0.05963223905 | 0.01 | miss |
| Zero-delay moments and insertions, absolute | 0.00101644319 | 0.002 | pass |
| Current after one operation, worst reached state | 0.03863321122 | 0.01 | miss |
| Current after coincident selective operations | 0.05754161285 | 0.01 | miss |
| Current after two coincident unitary kicks | 0.07059383968 | 0.01 | miss |
| Current after delayed selective operations | 0.024671983 | 0.01 | miss |

The continuum three-source predictor fails its own full-response comparison by about 5.95% before guard. The corresponding numerical tail estimates are far smaller. This is a failure of that finite source-span predictor for the continuum reference; it is not a failure of the full continuum calculation. Additional spectral content outside the three-source span is a plausible mechanism. Assigning the observed magnitude specifically to ultraviolet first-moment weighting would require a separate prospective decomposition.

The lattice predictor's success therefore supports a finite numerical model. It does not demonstrate continuum equivalence. A next derivation should retain the physical source prescription and account for omitted field modes and their memory effects, or derive physically justified Hamiltonian/source/boundary corrections before a new freeze. Retuning gains on the failed response data cannot solve the invariant discrepancies.

## Unitary-kick operating range

The table gives the worst guarded metric over both N22 runs. These are instantaneous unitary states U_kappa|Omega>. They are not the selected K+/- states or the nonselective detector channel obtained by ignoring the recorded outcome.

| kappa | Current disturbance from vacuum | Current mismatch to continuum | Both 1% criteria |
|---|---:|---:|---|
| 0.5 | 0.001204034358 | 0.006833360239 | pass |
| 1 | 0.004513406038 | 0.006337215259 | pass |
| 2 | 0.01770997881 | 0.01867627278 | miss |
| 4 | 0.06984522452 | 0.07059383968 | miss |

A postproduction sufficient bound follows from the exact finite-matrix structure. When each initial bond hopping h_j<=0, define A_j=-h_j v_j v_j^T>=0. Then

    0 <= S0-S_kappa <= delta S0,
    delta=max_j[1-cos(kappa w_j)] <= kappa^2 max_j|w_j|^2/2.

If e0=||Sc^(-1/2)(S0-Sc)Sc^(-1/2)||, the reference mismatch is bounded by e0+(1+e0)delta. All required signs and positive current matrices hold for the stored N22 inputs, and the registered cosine identities pass. With the unchanged 1% allowance and 1e-4 guard, the conservative sufficient radius across both runs is 0.7057128233124. Thus **|kappa|<=.7057** satisfies this current-only sufficient condition in the stored-matrix calculation. The actual registered kappa=1 point passes beyond this conservative sufficient interval.

This supplement evaluates algebra on archived results; it runs no extra physical trial and changes no gate. Its numerical evaluation is not an interval-certified exact-ground-state or continuum theorem. It does not bound selective branches, delayed responses, general states or experimental errors. The complete proof and reproducible calculation are in CURRENT_BOUND.txt and current_bound.py/json.

## Resolution and numerical controls

| N | Static covariance relative | Normalized response discrepancy | Lattice predictor residual | Continuum predictor residual |
|---|---:|---:|---:|---:|
| 14 | 0.120056727 | 0.113199926 | 0.00154390343 | 0.157994349 |
| 18 | 0.0670908465 | 0.0711461965 | 0.00061810254 | 0.0922368221 |
| 20 | 0.0527802205 | 0.0578660279 | 0.000418242177 | 0.074030353 |
| 22 | 0.0425553758 | 0.049063784 | 0.000292773846 | 0.0595322391 |

These are descriptive trends, not an extrapolated continuum bound. The same-size N22 refinement changes the archived comparison payload by at most 1.03472785895e-13. The largest mode-doubling change is 4.13385812248e-08; the largest conservative reference-tail estimate is 1.34674794368e-06. Step sources have finite first energy form; this study assumes no finite H^2 moment or strong residual bound.

The independent audit verified 21 source hashes, 5 checkpoint hashes, historical detector weights/boundary, branch contractions, source-current identities, diagnostics, guards and resource records. Its maximum assessment difference was 9.29811783124e-16; all tampering controls were rejected. This audit rechecks archived arithmetic and integrity, not a second implementation of the full evolution or Hamiltonian. The production run itself rebuilt H and checked ground residuals.

## The route to a physical working model: EX1-01

EXPERIMENT_PROTOCOL.md is a prepared hardware-agnostic experiment, not a Bell-style device-independent certification. It tests a finite Construction 1 operational prediction on four work qubits with external scheduling. Its independently verified ideal predictions are signal-action agreement 13/16 with coherence and 5/8 after parity dephasing: contrast 3/16. The immediate registered readout distribution is unchanged. The ideal verifier passes 22 checks and the raw-count scorer passes 12 synthetic controls. None of these are experimental observations.

The protocol fixes preparation, operations, parity-preserving dephasing, raw readout including erasures, control experiments, separate calibration/validation freezes, uncertainty requirements and rejection rules. Acquisition is capped at 1,000,000 calibration attempts plus 140,000 held-out validation attempts. It requires an independently justified apparatus confidence set and sampling model; a best fit or successful simulator replay cannot qualify the device. The hardware identity, compiled schedule, calibration model and adapter remain unfilled because no device was selected. No paid job or hardware call occurred.

A successful calibrated, held-out experiment would establish a finite physical realization of this operational model under standard quantum mechanics and test the restricted population-only comparator. It would not establish the autonomous clock model, Schwinger detector correspondence, continuum limit, gravity or a new GRUT law. The repository's integrated candidate adopts standard quantum dynamics and measurement rules; a distinctive physical restriction or selection principle still requires its own fixed prediction and independent empirical test.

The evidence therefore supports a working finite numerical predictor and a concrete physical test protocol. It does not yet support calling GRUT an empirically established fundamental theory.

## Provenance, preservation and reproduction

N22 remains its completed finite massless numerical pass under its adopted dictionary. Its component headroom 5.9324e-7 and combined headroom 2.5910e-4 are unchanged. B2's N20 component miss, D1's delayed-moment miss, withheld massive transfer and P1's interacting-parent obligations remain separate.

The ZIP contains the exact frozen source payloads, all raw results/audit/resource records, five historical binary inputs, EX1 protocol and software, and the postrun bound supplement. SHA256.json covers every payload except itself. For production reproduction use a clean checkout of the freeze above, Python 3.12, numpy 2.3.5, scipy 1.17.0 on Windows, and the supplied inputs; follow PROTOCOL.md. Source copies in the ZIP permit inspection and offline arithmetic audit but do not themselves contain Git metadata required by the production runner.

For an offline audit, preserve the delivered archive, make a copy of raw omitting its existing audit.json, and run:

    python source/research/reality_framework/grut_response_separation_20261109/independent_audit.py --root source --inputs inputs --output RAW_COPY

For the algebra-only bound calculation, use current_bound.py --raw raw --output NEW_BOUND_JSON. For the hardware ideal/scorer checks, follow EXPERIMENT_PROTOCOL.md. All result writers refuse overwrite. No massive test, larger-size extension, merge or release was performed.
