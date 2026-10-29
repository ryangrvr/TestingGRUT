# Corrections and interpretation changes

1. **Future gas heat and the 0.841% number.** In the immediately preceding `grut_m1_open_stress_20261029/REPORT.md` and coordination log, we said the proposed 0.841% result was arithmetically incorrect. That is too broad. At the saved `a/a_0=10` checkpoint, using the same hypothetical cold-conversion coefficient `f=0.3368502039`, the cold excess is 0.8222566687%, the gas heat is 0.0165556570%, and the held remainder is 0.0556240848%:

   * Retain present heat: `0.8222566687 + 0.0165556570 + f(0.0556240848) = 0.8575493100%`.
   * Fully redshift present heat: `0.8222566687 + f(0.0556240848) = 0.8409936530%`.

   The earlier criticism is correct only if the heat is retained; 0.841% is arithmetically correct under the additional full-cooling assumption. Neither scenario is established by a future evolution. The parent file is preserved for provenance; this correction supersedes the overbroad statement.

2. **Inflationary Hubble value.** The actual candidate stores `134284273491.30815 GeV`. The handoff's `134.284 GeV` is off by nine orders of magnitude.
3. **20 MeV radiation clock.** The parent expansion code actually computes `H/m_X = 1.9301e-21` at 20 MeV. The handoff's `3e-14` is not the ordinary 20 MeV value for this benchmark.
4. **Refit scaling.** The parent's fixed-dimensionless-history yield is proportional to `H_i^(3/2)`, so the 0.857549% mass-budget shift corresponds to a −0.567645% scale shift, conditional on that scaling and fixed other factors. It is not −0.29%.
5. **Curvature operators.** In four dimensions, `R Σ*Y²` requires one inverse mass and `R Σ*X⁶` requires five when the Wilson coefficients are dimensionless. Signs depend on whether the operator is written in the potential or Lagrangian; [REPORT.md](REPORT.md) states the convention used.
6. **Charge fraction, mass fraction, speed.** The inherited 96-bin checkpoint has a 6.595% free-charge fraction, a 6.974% cold-plus-held mass fraction, and `v_rms=0.06838c`; these are distinct quantities and depend on the assumed kinetic history.
7. **A geometric first-hit estimate is not a bound on coherent wave feedback.** Under the supplied uniform ballistic hypothesis, the proposed crossing estimate is sound as an average geometric estimate, not a theorem for correlated waves.
