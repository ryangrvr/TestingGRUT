# Independent finite-volume fixed-charge calculation

18 October 2026. Same Candidate 1 action, \(a=1,b=0.1\), and target \(\bar Q_*=399.4384624223012\). This is a separate numerical **method and implementation cross-check**, performed within the same project/team; it is not external independent validation or empirical evidence.

## Result and its scope

A cell-centered finite-volume energy minimization, using no stationary ODE/BVP solver, finds the same lower-frequency stationary configurations at three fixed charges. Two frozen starts at each of two cell counts produce matching energies. Relative to the refined BVP reference, the energy discrepancy shrinks by approximately a factor four when the cell spacing is halved. A subsequent, separately frozen one-step Newton check reduces the fine-grid stationarity residuals to near floating-point accuracy while preserving the original optimization records.

The numerical conclusion is agreement of two discretizations/methods on these **radial classical stationary states**. This does not establish that they globally minimize energy, exclude all possible fragments, prove a continuum error bound, show nonradial nonlinear stability, determine a quantum lifetime, or validate the particle sector against data.

## Frozen experiment and retained records

`independent_protocol.json` was written before the optimizations. It fixed \(L=40\), \(N=800,1600\), charges \(\bar Q_*,0.75\bar Q_*,0.5\bar Q_*\), both starts, analytic-gradient preflight, bounds, optimizer options and time budgets. No physical coupling was changed. All twelve requested runs completed in approximately 9.50 seconds; none was skipped. The objective and gradient code do not import the BVP or branch-scan implementation.

The two starts were:

- The existing target BVP shape, with its radius scaled by \((\bar Q/\bar Q_*)^{1/3}\), used only as initial fields.
- An independently specified broad Gaussian, \(F=1.25\exp[-r^2/(2\sigma^2)]\), \(S=1-0.8\exp[-r^2/(2(1.1\sigma)^2)]\), \(\sigma=4(\bar Q/\bar Q_*)^{1/3}\).

The maximum directional gradient-check discrepancy was \(3.66\times10^{-9}\) relative to the stated normalization. All original L-BFGS-B results returned a relative-energy `ftol` stop; **none** met the initially requested projected-gradient threshold \(3\times10^{-7}\). Their projected-gradient infinity norms ranged from \(4.06\times10^{-6}\) to \(7.60\times10^{-5}\), and the largest populated-cell unweighted Euler residual was about \(8.36\times10^{-3}\). These are explicitly retained stationarity limitations of the first pass, even though energies had converged closely between starts.

`independent_refinement_protocol.json` was then frozen before Newton refinement. It preserved the twelve original records and specified only the three \(N=1600\) Gaussian-origin outputs, at most six steps, the same bounds and physical parameters, sparse analytic Hessian plus rank-one solve, a Hessian-vector preflight and a 30-second limit. Every state accepted its first full Newton step; no backtracking or further iteration was required. The Hessian-vector preflight discrepancies were \(1.66\times10^{-10}\) and \(7.69\times10^{-10}\). The Newton loop took approximately 0.04 seconds, excluding Python startup.

## Discrete functional and analytic derivatives

With \(h=L/N\), cell centers \(r_i=(i+1/2)h\), and exact shell volumes

\[
V_i=\frac{4\pi}{3}\left[((i+1)h)^3-(ih)^3\right],
\]

the eliminated fixed-charge functional is

\[
\bar E_{\bar Q}=\frac{\bar Q^2}{2I}
+\frac12\sum_{i=0}^{N-2}c_i[(F_{i+1}-F_i)^2+(S_{i+1}-S_i)^2]
+\frac12c_{\rm out}[F_{N-1}^2+(S_{N-1}-1)^2]
+\sum_iV_i\left[\frac a4(S_i^2-1)^2+\frac12S_i^2F_i^2+\frac b4F_i^4\right],
\]

\[
I=\sum_iV_iF_i^2,\quad\Omega=\bar Q/I,\quad
c_i=4\pi[(i+1)h]^2/h,\quad
c_{\rm out}=4\pi L^2/(h/2).
\]

The outer values \(F(L)=0,S(L)=1\) are half a cell away. The origin has zero face area; there is no extra origin penalty. Charge is imposed through its exact eliminated kinetic energy, not obtained by solving at a guessed frequency.

Variables are \(u_i=\sqrt{V_i}F_i\), \(w_i=\sqrt{V_i}(S_i-1)\), with \(u_i\ge0\) and \(-\sqrt{V_i}\le w_i\le0\). Let \(L_V\) be the transformed gradient-energy matrix. Its interior off-diagonal entry is

\[
(L_V)_{i,i+1}=-c_i/(\sqrt{V_iV_{i+1}}),
\]

and its diagonal is the sum of adjacent face coefficients divided by \(V_i\), with the outer half-cell coefficient added at the last cell. Then

\[
g_u=L_Vu+\sqrt V(S^2F+bF^3)-\Omega^2u,
\quad
g_w=L_Vw+\sqrt V[a(S^2-1)S+SF^2].
\]

All products in this expression are componentwise. The exact Hessian is a sparse local block plus one rank-one term:

\[
H_0=\begin{pmatrix}
L_V+\operatorname{diag}(S^2+3bF^2-\Omega^2)&\operatorname{diag}(2FS)\\
\operatorname{diag}(2FS)&L_V+\operatorname{diag}(a(3S^2-1)+F^2)
\end{pmatrix},
\]

\[
H=H_0+\frac{4\Omega^2}{I}|(u,0)\rangle\langle(u,0)|.
\]

Sparse LU of \(H_0\), followed by Sherman–Morrison, gives the Newton update without forming a dense rank-one matrix. The local block need not itself be positive; omitting the rank-one term would solve the wrong fixed-charge problem. Linear solve residuals were below \(1.72\times10^{-18}\) in the recorded Newton steps. This stationarity check does not independently prove Hessian positivity.

## Quantitative comparison

The \(N=800\) values below are arithmetic means of both original starts; their energy differences are at most \(1.12\times10^{-10}\). The \(N=1600\) values are the separately refined Gaussian-origin states. No run was selected for best agreement with the BVP reference.

| Charge fraction | \(\bar E\), \(N=800\) | \(\bar E\), \(N=1600\), Newton | Refined BVP \(\bar E\) | \(\Omega\), \(N=1600\), Newton | BVP \(\Omega\) |
|---:|---:|---:|---:|---:|---:|
| 1 | 363.9174355292 | 363.9216078576 | 363.9229985516 | 0.799997913660 | 0.800000001972 |
| 0.75 | 282.6228841783 | 282.6263217720 | 282.6274675594 | 0.830505996789 | 0.830508938618 |
| 0.5 | 197.2954383785 | 197.2976750889 | 197.2984205537 | 0.885069028188 | 0.885074813117 |

The fine-grid energy discrepancies are approximately \(-0.001390694,-0.001145787,-0.000745465\). They are ordinary discretization discrepancies, not evidence for a distinct branch or lower continuum energy. Relative to the BVP values, the coarse-to-fine error ratios are about \(4.0002,4.0002,4.0004\).

Assuming a leading \(h^2\) error, \((4E_{1600}-E_{800})/3\) gives \(363.9229986337,282.6274676366,197.2984206591\), within approximately \(8.22\times10^{-8},7.72\times10^{-8},1.05\times10^{-7}\) of the BVP values. This is a Richardson consistency check with an assumed leading error form, **not** a rigorous error enclosure or evidence from a third independent discretization. Both minimization grids use \(L=40\); this calculation adds mesh refinement, not a new domain refinement.

The fine-grid stationarity checks are:

| Charge fraction | Full gradient infinity norm | Projected gradient infinity norm | Populated-cell \(F\) Euler residual | Populated-cell \(S\) Euler residual |
|---:|---:|---:|---:|---:|
| 1 | \(4.37\times10^{-12}\) | \(4.37\times10^{-12}\) | \(2.85\times10^{-12}\) | \(6.35\times10^{-13}\) |
| 0.75 | \(5.08\times10^{-12}\) | \(5.08\times10^{-12}\) | \(4.67\times10^{-12}\) | \(1.19\times10^{-12}\) |
| 0.5 | \(9.14\times10^{-12}\) | \(9.14\times10^{-12}\) | \(1.10\times10^{-11}\) | \(5.96\times10^{-13}\) |

The Newton energy changes from the original fine-grid optimizer outputs were only about \(10^{-10}\), showing that the original energy agreement was robust even though their stop criteria had not certified stationarity.

For the specific equal-halves branch competitor, the independent energy difference \(2E(\bar Q_*/2)-E(\bar Q_*)\) is \(30.67344123\) on the coarse grid and \(30.67374232\) on the fine grid. The BVP value is \(30.67384256\). This supports the positive gap for that particular same-branch channel; it does not exclude every split, every remnant, a disconnected lower branch, or quantum decay.

## Reproduction and boundaries of interpretation

Run from the repository root, with one BLAS thread if desired:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python research/reality_framework/grut_candidate_1_stability_20261018/independent_minimize.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python research/reality_framework/grut_candidate_1_stability_20261018/independent_refine.py
python research/reality_framework/grut_candidate_1_stability_20261018/independent_analyze.py
```

The first two scripts refuse to overwrite existing experiment outputs. To reproduce the experiment, use a fresh copy/output location while preserving the protocol and source inputs. Source hashes, implementation hashes, versions and stop records are saved with the outputs. `independent_comparison.json` describes the original optimizer pass and retains its stationarity caveat; `independent_refinement_results.json` separately records the successful Newton correction. Nothing in the later refinement erases the first-pass stop behavior.

The field class is radial with \(F\ge0,0\le S\le1\). Two starts do not exhaust that class. An excellent discrete root need not be a globally minimizing state, and solving a radial energy problem does not test arbitrary three-dimensional perturbations. The finite-volume/BVP agreement addresses method dependence of the reported stationary energies and frequencies. The separate evolution, branch/fission scan, quantum corrections and possible physical preparation remain separate obligations.
