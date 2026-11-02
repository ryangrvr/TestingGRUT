# Brief: the decisive calculations for M1r (2 November 2026)

The owner asked us to "work it all out". Read [WORKING_THEORY.md](WORKING_THEORY.md) and [BRIEF.md](BRIEF.md) first;
the earlier rules apply. This pass computes the items of WORKING_THEORY.md §8:
1. The end-to-end hidden entropy per charge S/Q of the actual M1r release, which decides the branch: free X if
   S/Q > 3.47, Q-matter if S/Q < 3.47.
2. The Q-matter branch's microphysics and its observable: sticking, saturated vapour, ⟨Q⟩, vapour fraction and the
   transfer function against Lyman-α and satellites.
3. Higgs survival through inflation and reheating with ξ_H = 0.1.

## M1r (the repaired model)

**Light sector.**
- m_X = g v = 0.0924089637 GeV, with g = 0.1 and v = 10 m_X.
- λ_X = 10⁻³; λ_φ = 0.03.
- V_φ = λ_φ/4(φ² − v²)² − εv(φ − v)³, with ε = 0.002.
- Object family: E(Q), Ω(Q) from tracks/fx_formation/code/family.py (Q_min = 56709; binding above 71263); bulk
  Ω_b = 0.775865 at 449.31 m_X³; spinodal 233.2 m_X³; surface tension 98.8 m_X³.

**Heavy sector.**
- M_Y = 30 H_i, h = 40 H_i, μ = 7.159456 H_i, κ = 3.4586847×10⁻³.
- **λ_Y = λ_Σ = 0.5.**
- Hidden U(1) coupling e: scanned here. The declared value 10⁻⁶ is the free-X branch; the Q-matter branch needs
  e ≳ 2×10⁻⁶.

**Holder.**
- −(g_I e^{iα_I} φ_I Σ*Y² + h.c.), with g_I = 8.306742×10⁻⁷ (calibrated at λ = 1) and α_I = 2.192144074.
- Prepared state at λ = 0.5: (90.800, 6.888, 28.576) H_i, lightest mass² 10.19 H_i².
- Code: tracks/fx_holder/code/rptri.c (a copy of our rp_reheating engine with μ_eff = μ + Gφ_I/H_i). Known bug:
  forced-Y uses 2μ_eff where 2μ_eff* is needed; it affects early N_YL/N_YH only. Fix it in your copy.

**Inflaton and gravity.** (m_I, Γ_I) = (3, 0.5) H_i; ξ_X = −5/6; ξ_H = 0.1; zero SM portals.

**Induced-portal renormalisation condition.** Adopt "vanish at μ_R = M_Y" as the stated benchmark. Report the
sensitivity to "vanish at M_P" (Codex's skeptic found −17% in H_i).

**Release at λ = 0.5 with the trilinear holder on (3, 0.5).**
- Y/Y_Q = 0.990733; N_YL = 0.567, R_Y = 0.379.
- ε (charge/energy) ≈ 0.127.
- H_i calibration: 1.29×10¹¹ GeV on the free-X relic, 1.50–1.52×10¹¹ on Q-matter. Dimensionless histories are
  H_i-independent; physical timing (decays, onset) uses H_i.

## Known results to build on (do not redo)

- **Heat budget (tracks/fx_chain).** At the old spot release: total R = 4.9665 (band) / 5.7033 (broad) at e = 10⁻⁶,
  split as vectors 4.0538 / 4.7905 and mediators 0.9127 (Y_L 0.8056 + cascade vectors 0.1072). The skeptic found
  the cascade vectors scale as e^{−1/2} at fixed N_YH.
- **Vector production (tracks/hw_gates/code/gate1_vector.py and gate1b_pieceB.py).** Longitudinal Proca modes:
  piece A (inflationary, Graham–Mardon–Rajendran) and piece B (reheating era).
- **The S/Q machinery.** Relaxed Hartree plasma (tracks/fx_chain/code/heat_onset.py; skeptic hartree_indep.py;
  wt_critic/code/hartree_copy.py and chain_consistency.py, where c = N_c/E^{3/4} is universal); survival threshold
  tracks/fx_chain/skeptic/code/survival.py (3.4716 bulk, 2.93 at Q = 2×10⁵).
- **Kinetics and vapour (tracks/fx_formation).** Kinetics (Fuchs mergers, Wagner ripening, capture freeze-out),
  vapour thermodynamics (skeptic bessel_vapour.py) and the free-X fraction 1 − η ≈ 1/[Y_h(2.5 + 0.224/T_f)].
- **Transfer (tracks/fx_observables).** vapour_transfer.py, plus the skeptic's two_fluid_check.py; tight coupling
  fails for ζ ≳ 5×10⁻⁴.
- **X–core elastic scattering (tracks/m1_transport/skeptic/code/sk_xcore.py).** Three-channel radial solver:
  σ ≈ 332/m_X² at k = 0.068 m_X.
- **Higgs and SM running (research/reality_framework/grut_closure_forcing_20261027/tracks/cf_vacuum).** SM RGE and
  V_eff code; Λ_I and h_max = 2.7×10¹¹ GeV at M_t = 172.6.

## Rules

- Write only in `research/reality_framework/grut_working_theory_v2_20261101/tracks/<key>/` (NOTES.md, code/,
  results/*.json; nothing over 2 MB); scratch in /tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/<key>/.
- Git is read-only. No commits.
- At most 2 threads or worker processes (other agents share the 4 cores). Lattice runs are allowed only in tracks
  that say so, and no larger than 64³. Run in the foreground (Bash timeout up to 600000 ms per call; chunk with
  restart files); never end your turn with a background job running.
- Claim labels [D] [C] [Q] [E] [A] [R]. Choices are model choices, not consequences of closure. Never claim
  validation, exclusion of the real world or a new theory. Data comparisons are retrospective; cite with arXiv
  numbers. For current limits use WebSearch via ToolSearch "select:WebSearch".
- Return a digest of at most about 700 words.
