# Verification actually performed

One author; no subagents or independent replication. No outside review of
this proposal is claimed. Claude's owner-relayed review concerns the
parent identity, not this calculation; his old 51 checks were not rerun.

## New mathematical checks

- Derive the radial static-source qualification directly from its field
  equation and the parent mass/depletion identity. Evaluate both source
  integrals on all three parent profiles. Maximum relative source-identity
  residual is 5.9e-13. This checks a defined static quantity, not a newly
  calculated quantum emission amplitude.
- Derive the proposed Higgs source from the portal operator and independently
  by completing the uniform-background potential square and differentiating
  the single-scale mass. The source integral and mass formula agree within
  1.6e-12 on the stored profiles. No new portal saddle was solved.
- Normalize the nonrelativistic amplitude with reduced mass, confirm its
  contact differential recoil cross section integrates to the total, and
  check the common mass-normalized cross section across all three charges.
- Derive the positive-source sufficient Born-kernel bound. At the declared
  kappa it is <=1.94e-4. The formal 1/10 GeV one-event examples remain below
  .002; the original high-scale extrapolation does not define a valid point.

## New numerical checks

- Hash all three frozen profile inputs. Check source normalization, F(0),
  small-q moment expansion, radial-grid decimation and a different
  PCHIP-plus-adaptive-quadrature evaluation at q/(g_X v)=1.
  Largest decimation difference is 1.15e-12 in F; largest interpolation/
  quadrature difference is 2.51e-8. This does not estimate omitted physics.
- Normalize the shifted/truncated Maxwellian speed distribution by direct
  quadrature. Check the analytic inverse-speed formula at seven thresholds
  against direct integration, including its escape edge.
- Check raw spectrum factorization, its integral, and recoil-grid decimation.
  Halving the reference recoil grid changes the 5-40 keV integral by
  4.55e-11 relative.
- Fresh isolated re-execution reproduces all 378 numerical JSON leaves,
  nonnumerical provenance and the raw 1401-by-6 reference CSV byte for byte
  in this environment. The check uses the same author, implementation and
  libraries; it is not an independent physical-method replication.

`code/verify.py --replay` passes 97 internal conditions before the manifest
exists, or 99 including its two content-integrity conditions afterwards.
Check counts describe the script, not a confidence level or validated theory.
Python syntax, JSON, local report links and git whitespace were also checked.
The manifest covers only this new packet, not prior scientific records or
the append-only coordination log.

Environment: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. No new dependency
installation. Calculator and verification need no network.

## Explicitly not performed

No experimental efficiency/resolution or isotope likelihood, backgrounds,
accepted-event forecast, exclusion region, held-out test, new discovery,
cosmological formation/abundance, portal heating, long-lived-scalar decay
network, BBN/CMB analysis, finite-portal saddle, new determinant, complete
renormalization-group evolution, full nuclear response or quantum lifetime.
No total vacuum matching was changed. No old result or canonical GRUT-RAI
file was edited, and no empirical target or new action was adopted.
