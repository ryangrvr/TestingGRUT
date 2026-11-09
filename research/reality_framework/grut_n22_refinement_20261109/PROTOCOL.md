# N22-01: one frozen same-rule massless refinement

This follow-up responds to the [research review](https://github.com/ryangrvr/TestingGRUT/pull/14#issuecomment-6733231617). Base: B2 result commit 1fa9f643e01f2b1a012e82cb8061592cba2a40ef. B2's combined measured budget passed (0.000887384268 < 0.001); its separate lattice component failed. Both statements and the withheld N20 massive transfer remain unchanged whatever N22 finds.

The authoritative parameters are in protocol.json. This directory imports B2's model.py and point.run without modifying either. Recompute weights at a=5/23 with the identical quadratic stencil, boundary values and normalized Gaussian. Do not reuse the N20 coefficient array. L=5, sigma=0.8, coupling=2, phase=0, m_R=0, m_lat=-a/8, preparation rule, Hamiltonian, ghost endpoints, and all seven timing sequences are unchanged. The previously named confirmation timings are inherited validation cases, not new blind timings for this successor.

## Registered execution and stop rule

Run exactly one N22 standard solve and one N22 refined solve, sequentially. The neutral-sector dimension is 705432, versus 184756 at N20. Standard eigensolver tolerance is 1e-11; refinement uses 1e-13, an independently specified starting-vector seed, and two subintervals per propagation interval. This changes the numerical solve, not the physical ground-state preparation. Continue to compute 8192- and 16384-mode references and conservative analytic Gaussian tails.

Limits: 900 seconds per point, 2000 seconds total, one BLAS thread, and 4 GiB committed memory per child process. Windows Job Objects enforce the memory limit; sample resource telemetry and emit a heartbeat every 30 seconds. A miss, exception, nonzero exit, or resource-cap failure is archived. Stop without automatic retry, N24 extension, changed weights, reassigned allowances or massive run. The runner is Windows-specific; the scientific model and analysis are ordinary NumPy/SciPy. Memory implementation follows [Microsoft's documented extended job limits](https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_extended_limit_information).

Before committing, run inherited algebra checks at N22, the synthetic oscillator/N6 checks, import/syntax checks, and toy child-process tests of normal exit, memory-allocation rejection and timeout termination. Test synthetic decision boundaries and a sign-reversed forecast. These prechecks do not inspect N22 production probabilities. Commit source, protocol and prechecks before production; run only from that clean commit and record source/input hashes. Preserve raw probabilities, ground states, peak memory, logs and failure output. Outputs are external and never overwritten.

## Nominal thresholds and numerical margin

Both N22 solves must satisfy the original component thresholds separately: geometry TV <=0.00025; lattice-to-matched TV <=0.000740; direct target TV <=0.001; matched+geometry+0.00001 <=0.001. Keep the original normalization, eigensolver, reference-tail, doubling, refinement, disturbance and finest-step thresholds. The new finest step is N20-to-N22, checked against both solves. Validate resolution metadata and every outcome distribution.

Report a nominal pass only if all these nominal gates pass. Additionally report a numerically separated component pass only if, for each history,

    max(standard matched TV, refined matched TV) + 0.00001 <= 0.000740.

This uses the full existing numerical reserve as a conservative operational margin guard on the component; it does not enlarge or redistribute a component allowance. If nominal gates pass but this guard fails, label the outcome `nominal_pass_margin_unresolved`. A miss of the nominal threshold remains a miss even if other components have spare allowance. The reserve and refinement are numerical diagnostics, not a rigorously proved probability-error bound. Report both low eigenvalues and their computed gap; do not treat the gap or a bare residual as a certified ground-state error bound.

## Signed forecast diagnostic

Preserve every delta_N(h,s)=p_lat,N(h,s)-p_matched,N(h,s) and delta_N/a_N^2 for N10,14,18,20 and both N22 solves. Freeze the forecast delta22_pred=delta20*(21/23)^2, without subtraction from any measured result. For each history and each N22 solve, calculate

    R_h = half_L1(delta22/a22^2 - delta20/a20^2) / half_L1(delta20/a20^2).

The forecast diagnostic passes only if all fourteen comparisons have R_h<=0.10; for a zero denominator require numerator<=1e-10. This 10% tolerance is a prospective diagnostic choice, not derived physics, and is separate from the finite probability gates. Report every history, not just the maximum. No fitted or extrapolated correction enters any detector or probability.

From inspected N20 data, the largest predicted matched TV is 0.000729507006 at [0,2]. Its nominal headroom is 0.000010492994; after the 0.00001 guard the forecast leaves only about 0.000000492994. This narrow prediction is not a bound or a pass.

## Interpretation and transfer restriction

An N22 pass supports only the registered finite massless experiment under the adopted lattice-to-cell-average dictionary. It does not derive that dictionary, establish an asymptotic rate, certify the ground state rigorously, or show a distinctive GRUT law. A forecast pass is separate evidence about the observed signed trend. A failure is preserved without describing it as a total-budget violation when the measured total passes.

No massive run is authorized by N22-01, even if all numerical gates pass. B2's massive rule was specifically N20; N22 cannot retroactively clear it. Any later massive study must state its own resolution and dictionary. P1's interacting-parent cutoff and overflow obligations remain separate. No merge, release or Zenodo deposit.

## Reproduction

Use Python 3.12 with numpy==2.3.5 and scipy==1.17.0. From a clean source freeze on Windows:

    python research/reality_framework/grut_n22_refinement_20261109/run_n22.py --output NEW_EXTERNAL_DIRECTORY

No existing output directory is accepted. execution.json distinguishes process completion from the scientific verdict in assessment.json.
