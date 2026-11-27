# C1 fixed-charge quantum channel comparisons

**Verdict A within the frozen scalar one-loop calculation.** The original Candidate 1 benchmark has negative energy differences for emission of one elementary carrier and splitting into two identical half-charge solitons. Both intervals include the registered numerical allowances. The fixed-charge continuum saddle premise, exact continuum error bounds, higher-loop effects and external replication remain open at this stage.

The branch was frozen at [a3ff9f4](https://github.com/ryangrvr/TestingGRUT/commit/a3ff9f40a1588501a24afc87de996668ef237374) on [S1 draft PR #29](https://github.com/ryangrvr/TestingGRUT/pull/29). The original session ended during packaging. This packet retains the completed freeze-pinned replay recovered from the saved archive. The archived energies and production log were recovered byte-for-byte. The six profile splines were regenerated using only the stationary solves and checked against the recovered run. [Protocol](PROTOCOL.json) · [Reproduction](REPRODUCTION.md)

## Channel results

All energies below are divided by $g_Xv$. The strictly expanded elementary pole mass is $1.0001232694040532$.

| Q | Classical energy | Scalar one-loop correction at retained cutoffs | Numerical allowance on the sum |
| ---: | ---: | ---: | ---: |
| 39,944 | 36,392.42286136 | 2.78783262 | 0.02668966 |
| 39,943 | 36,391.62286052 | 2.78777941 | 0.02668966 |
| 19,972 | 19,729.91009896 | 1.65541564 | 0.02668972 |

The frozen convention is initial minus asymptotic final energy,
\[
\Delta_1=E(39944)-E(39943)-m_{X,\mathrm{pole}},
\qquad
\Delta_{1/2}=E(39944)-2E(19972).
\]

| Comparison | Center | Signed interval |
| --- | ---: | ---: |
| One-carrier emission on matched grids | −0.2000692023 | [−0.203320826, −0.196817579] |
| Integrated chemical potential | −0.2000692023 | [−0.201140826, −0.198997579] |
| One-carrier emission with separate absolute allowances added | −0.2000692023 | [−0.253448528, −0.146689876] |
| Equal fission into the specified daughter branch | −3067.92033523 | [−3068.00040433, −3067.84026613] |

The archived tree comparisons were −0.1999991521 and −3067.39733656. For emission the matched spectral difference is −0.000018492 and the restored-graph difference is +0.000071711. The direct and integrated chemical potential centers differ by $3.71\times10^{-11}$. The integrated check varies the stationary background at five frequencies and uses tangent boundary-value solves for $\mathrm dQ/\mathrm d\Omega$. Both routes share the physical determinant, subtraction and matching. [Raw results](RESULTS.json)

## Verification of the calculation

The original action and matching use $g_X=0.1$, $\lambda_\phi=0.01$, $\lambda_X=0.001$, $m_{D0}/(g_X v)=1$, $h_D=0.01$, $\lambda_D=0.001$, $v/M_P=10^{-5}$ and $\mu_{\overline{\rm MS}}/(g_X v)=1$. The calculation keeps the leading flat scalar saddle and new-scalar loops. The physical branch was continued at the three integer charges from the original hashed profile. Each branch passed the registered charge, virial, monotonicity, tangent and finite-grid reduced-spectral controls. The half-charge saddle has $\Omega=0.8850741952$ and $\mathrm dQ/\mathrm d\Omega=-124417.66$. [Original action](../grut_candidate_1_20261018/GRUT_CANDIDATE_1.md) · [Raw results](RESULTS.json)

At fixed charge the Routh reduction of the cyclic phase adds $4\Omega^2|e\rangle\langle e|$ to the radial amplitude stiffness. Eliminating the normalized zero coordinate in the Euclidean pencil yields $\det P_{\rm full}=\zeta^2\det P_{\rm reduced}$. The original complete vacuum sets the relative normalization. Phase and translation modes contribute the restored infrared logarithm with total multiplicity four. In the Legendre transform the frequency shift multiplies $K_0'(\Omega_0)+Q=0$, so the first-order energy retains only $K_1(\Omega_0)$ on the classical saddle. The carrier mass uses the strict first-order pole expansion at the same matching. [Routhian and Legendre derivation](../grut_candidate_1_energy_20261019/ENERGY_DERIVATION.md) · [Pole matching](../grut_candidate_1_fluctuations_20261019/POLE_THRESHOLDS.md)

For the reporter, the quadratic form of $-\partial_\rho^2+\ell(\ell+1)/\rho^2+1+0.01S^2$ is bounded below by the identity. The vacuum block has mass squared $1.01$ and the relative potential is $0.01(S^2-1)$. This relative potential enters the determinant and both Born terms. The restored graphs use the fixed $\overline{\rm MS}$ convention and laboratory frequency harmonics. The pole mass includes the reporter tadpole and vacuum shift. [Independent critique](INTERNAL_CRITICAL_REVIEW.md)

## What the allowances support

The absolute allowances combine analytical omitted-frequency and omitted-angular estimates on the supplied truncated profiles (about 0.014870 per energy) and 0.011820 in empirical allocations for radial, frequency, support, graph and classical convergence. The paired emission allowance combines the paired-tail estimate 0.0000516234, empirical allocations totaling 0.0032 and a pole allowance of $10^{-8}$. The half-charge errors enter the fission interval twice. The separate-allowance result adds $0.0266896579656+0.0266896579674+10^{-8}$ to the emission center computed from the absolute energies.

The analytic tail estimates use profile envelopes evaluated in floating-point arithmetic on numerical splines. The discretization and support allowances are empirical. Validated exact-profile and half-line error control would be a separate result. Higher-loop, gravitational and nonperturbative corrections fall outside the intervals. [Independent critique](INTERNAL_CRITICAL_REVIEW.md) · [Negative results and limitations](NEGATIVE_RESULTS.md)

The verdict applies to the specified one-carrier channel and equal split into the continued branch. Asymmetric fission, other partitions, a lifetime, and formation require separate physical results. [External RAI review packet](RAI_REVIEW_PACKET.md)
