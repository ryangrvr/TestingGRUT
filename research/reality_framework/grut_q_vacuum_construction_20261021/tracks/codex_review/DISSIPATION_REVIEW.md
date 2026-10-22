# Dissipation verification and handoff correction

22 October 2026 UTC. Codex, one author. This is a new targeted review
of the published `78352e2` packet, not a resumed Claude agent or a claim
that all partial verifier programs completed. Its decisive equations
also agree with the earlier separate `ca2c867` build.

## Dispositions

| Claim | Disposition | Deciding check or limit |
|---|---|---|
| A varying effective mu can be inserted while matter stays separately conserved | Incorrect without extra structure | Receiver term Q=q mu_eff_dot is required; at constant K, d(3H²-rho_V-rho_r)/dt=0 only with both opposite source terms |
| KV16a approaches Minkowski on the positive radiation-plus-vacuum branch | Confirmed for the specified source | Closed finite-N law; four independent cosmic-time integrations with two initial radiation densities |
| Radiation removes the attractor | Not the result of these equations | Vacuum at a fixed N is independent of initial radiation density; the receiver follows a separate retarded integral |
| Natural-rate relaxation prepares today's density in the stated thermal history | Fails | 9.91e122 e-folds from M^4, compared with a reference entropy span 72.5; after BBN self-drain is negligible |
| 9M^4/(4 gamma N) is a lower floor | Corrected: upper envelope | rho_Vi/[1+A rho_Vi N] < 1/(A N), A=4 gamma/(9M^4) for rho_Vi>0 |
| Positive attractor at reference age gives observed density | Fails for the assumed early-start attractor history | 2.2435e35 GeV^4, about 8.90e81 reference densities |
| Any negative branch is restored by the same law | Refuted by its sign | Q>=0 implies rho_V never increases; negative-vacuum control reaches H=0 |
| KSV16 selects Minkowski uniquely | Fails | At K=f/3, epsilon=-f²/2+f^4/6, static H²=epsilon'-epsilon/f and u=2epsilon/f-epsilon'; the source vanishes for each of these de Sitter states |
| K22 Table III can be reproduced | Confirmed in new variables | Evolve H,r_V,r_r with Raychaudhuri and both energy equations; reconstruct xi=2+sqrt(r_V) at mu=4, lambda=4 |
| K22's late-time branch enforces its original algebraic action equation | Does not follow | For general lambda at mu=4 the r_V=0 endpoint has epsilon'(xi)-mu=2sqrt(4-lambda), nonzero except lambda=4 |
| K22 H^5 stops when vacuum energy reaches zero | False with positive radiation | Q=g_tilde H^5>0 at rho_V=0 and rho_r>0 |
| The basin's coordinate-box fraction is a probability of our universe | Unsupported | No initial-state measure is derived |
| The fitted finite-range separatrix exponent proves 140 digits of tuning | Unsupported extrapolation | No asymptotic error control over the quoted extrapolation; measure-zero selection in that model is the supported structural point |
| The measured approximate (H/m_q)^6 deviation is universal | Unsupported outside the tested model/preparations | Keep the finite-run result; it is not a class-wide theorem |

## Independent calculation

`code/review.py` imports none of the original track programs or the
earlier Codex code. In M=1 units it evolves

\[
\dot H=-2\rho_r/3,\quad
\dot\rho_V=-4\gamma|H|\rho_V^2/9,\quad
\dot\rho_r=-4H\rho_r+4\gamma|H|\rho_V^2/9,\quad\dot N=H.
\]

For rho_Vi=0.1 or 1, rho_ri=0.2 or 10 and N_end=12, the closed vacuum
flow agrees within 5.6e-16 relatively. The independent receiver
quadrature agrees within 1.5e-9. The maximum sampled normalized
constraint is 1.1e-11. One tighter-tolerance run is included. These
figures describe these integrations, not a physical error budget.

For K22 with G=1, Q=g_tilde H^5, Raychaudhuri becomes
Hdot=-16 pi rho_r/3. At tau=1e6 the new integration gives
H=0.04301238, r_V=0.0002208344 and r_r=8.5569e-10,
reproducing the printed Table III. The original partial verifier's
large basin grid is not rerun or declared complete.

The BBN bound in `review_results.json` is for vacuum **self-drain**
on the radiation-only law, with initial rho_V no larger than the
reference radiation density. Plasma traces and later potential changes
require their own consistent history and can add signed offsets. No
universal probability or no-go claim for all reheating histories follows.
The extrapolated Planck-scale examples are diagnostics of this specified
effective law, not controlled quantum-gravity calculations.

The positive branch can relax with radiation present, and it can contain
a prepared standard radiation era. The reviewed failure is that this
law does not select the observed residual in the worked finite history.
That is the replacement for the handoff shorthand "fails with radiation".

## Reproduce

Python needs numpy, scipy and sympy. This review used an isolated
SymPy 1.14.0 installation in scratch; it did not modify the runtime's
global packages. Versions are in the output.

```sh
python tracks/codex_review/code/review.py
```

Run from the construction directory. This writes only
`tracks/codex_review/results/review_results.json`.
All 17 mathematical/numerical dispositions passed in the recorded run.
