# RC05_CASE_06 — Chiral perturbation theory (ChPT) quantitative prediction — positive effective calibration

> REALITY_CHECK_05 calibration case; frozen protocol, unmodified.

## Sources

- **J. Gasser, H. Leutwyler, Ann. Phys. 158, 142 (1984); Nucl. Phys. B 250, 465 (1985)** — one-loop ChPT: pion mass relation m_pi^2 = M^2[1 - (M^2/F^2)(L-bar log M^2 + k_pi)] with L-bar = ln(mu^2/M^2) `CANONICAL_PRIMARY`
- **G. Colangelo, J. Gasser, H. Leutwyler, Nucl. Phys. B 603, 125 (2001)** — pi-pi scattering lengths predicted at two loops with error estimate; later confirmed by experiment `CANONICAL_PRIMARY`
- **S. R. Beane et al. (NPLQCD), Phys. Rev. D 77, 014505 (2008) / DIRAC experiment (Phys. Lett. B 619, 70, 2005)** — lattice/experimental determination consistent with ChPT pi-pi scattering prediction `CANONICAL_PRIMARY`

## Input ledger

- low-energy constants: F_pi, Sigma → EMPIRICAL_INPUT — fixed from OTHER observables (pi-pi scattering / lattice), not from the predicted observable
- chiral symmetry group SU(2)_L x SU(2)_R → ASSUMED (QCD symmetry, empirically grounded)
- expansion parameter (M^2/F^2), loop order → ASSUMED/DECLARED truncation with error estimate
- renormalization scale mu → DECLARED convention; physical result scale-independent at fixed order

## Model reconstruction

Effective Lagrangian built from Goldstone fields with chiral symmetry; perturbative expansion in M^2/(4 pi F)^2; observables computed as loops + counterterms with LECs fixed independently.

## Provenance map

prediction -> loop+counterterm computation -> LECs fixed from independent observables -> chiral symmetry (from QCD).

## Identifiability

Observable uniquely determined at declared order once LECs are fixed; truncation error estimated systematically.

## Prediction analysis

Satisfies all seven protocol prediction criteria: coefficients/structure fixed before comparison; inputs independently declared (LECs from other observables); no fitting to the target; observable specified; quantitatively testable; derivation traceable. PREDICTIVE_WITHIN_SCOPE = YES. pi-pi scattering lengths (Colangelo-Gasser-Leutwyler 2001) subsequently confirmed — a genuine effective-theory prediction.

## Effective status

EFFECTIVE — with declared regime (M^2/(4 pi F)^2 << 1) and systematic truncation uncertainty.

## Adversarial attack

Would the instrument call ChPT 'not predictive' because it is effective? No: effectiveness does not preclude prediction where the seven criteria are met — this is exactly the calibrated distinction. Would it reject because LECs are empirical? No: the LECs are fixed from observables OTHER than the predicted one. Both attack vectors failed.

## Final classification

```json
{
  "effectiveness": "EFFECTIVE",
  "derivation": "DERIVED (within declared order and declared inputs)",
  "prediction": "PREDICTIVE_WITHIN_SCOPE",
  "inputs": "EMPIRICAL_INPUT (LECs) + ASSUMED (symmetry, order)"
}
```

Confidence: HIGH | Instrument pass: True

Limitations: Sources canonical-primary, not re-fetched.
