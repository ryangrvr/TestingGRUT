# Stage 2 brief: solve the construction, derive its consequences, audit

Read BRIEF.md first. All of its rules apply. This addendum gives the stage-1 results and defines the construction
that stage 2 solves.

## Stage-1 results (the records are tracks/literature, tracks/equations, tracks/ledger)

- **Equations.** The brief's equations are confirmed two ways: a blind derivation and the KV papers read at full
  level from TeX.
  - (q) eps'(q) − ½K'(q)R = mu, a single global constant.
  - (g) K G = T − (eps − mu q) g − (g Box − nabla nabla) K.
  - Off-shell identity: nabla·E = −q d mu. Matter is conserved exactly when mu is constant.
- **Scalar-tensor form.** For fixed mu the theory is Brans-Dicke with omega_BD = 0 (equivalently f(R)):
  - m_q^2 = (2/3) K eps''/K'^2;
  - G_eff(r) = G*(1 + e^{−m r}/3);
  - natural scales give m_q ≈ 2e18 GeV.
  - For m_q ≫ H the theory reduces to GR + V_m(t) − c*(mu) + O(H^2/m_q^2), where c*(mu) is the matter vacuum level
    at which that mu is in equilibrium.
- **Energy origin.** With mu = mu0[eps + c], rho_V,eq = 0 for every c, and dmu0/dc = 1/q0: the tuning is carried by
  the state.
  - With a fixed number mu, c gravitates one-for-one.
  - mu is conserved through transitions. The final vacuum is Minkowski only if mu = mu0 of the **asymptotic, final,
    fully renormalized** eps.
  - Pre-transition or initial-equilibrium readings leave −ΣΔ. The electroweak step alone is −4.7e54 rho_Lambda,
    giving an AdS crunch.
  - Re-imposing zero pressure after each transition is inconsistent with conserved mu.
- **KV themselves** choose mu = mu0 of the final full eps (1102.3152 §6) and concede that "the fine-tuning moves to
  mu".
- **The KV electroweak "kick"** rho_V ~ E_ew^8/E_P^4 (the same scaling as AHHKM astro-ph/0005111) is not derived from
  the action.
  - Its classical part, rho_V0 = (9/8) kappa^2 H^4 K'^2/eps'', is a transient proportional to R^2.
  - The freezing into a constant rests on KV's phenomenological relaxation law (0905.1919 eq (4.1)), which Klinkhamer
    later called hypothetical.
  - KV calibrate E_ew to the observed value (3–9 TeV, 4.7 TeV).
- **In the action**, a second-order term never survives as a constant:
  - the transition oscillation energy is (4/3)Δ²/(K m_q²)(x/sinh x)² with x = π m_q τ/2, and m_q/H ≈ 1e32 at the
    electroweak scale;
  - the matter-induced term is rho_V ≈ rho_m²/(12 K m_q²) ≈ 5e-169 GeV⁴ today.
- **GRUT contribution.**
  - "No exterior => zero vacuum pressure" is not derivable from FORMAL_CORE. Nothing is GRUT-native.
  - GRUT supplies wording, slot bookkeeping and one structural restatement: THM 2.15 (protect or form). If mu is
    conserved it must be prepared; if it relaxes, it is not protected.
- **Exposed data.** tracks/ledger/results/ledger.json: G-dot, BBN/CMB G, Cassini, Eöt-Wash, DESI DR2 w0wa, early dark
  energy, Omega_c h², FIRAS. The unscreened Yukawa coupling alpha = 1/3 requires m_q ≳ 4 meV.

## The construction to be solved (QV1) and its labelled variants

**QV1 (credited to Klinkhamer–Volovik).** It consists of:
- the action of BRIEF.md;
- matter = Standard Model plasma plus dust and radiation, minimally coupled;
- the state condition (P_q): mu = mu0[eps_inf], where eps_inf = eps(q) plus every q-independent vacuum energy of the
  final T → 0 vacuum, to all orders;
- generic eps(q), K(q) with natural scales (q0 ~ M², q0² eps'' ~ M⁴, q0 K' ~ M²). KV's 0806.2805 forms are used where
  a concrete choice is needed.

Two labelled variants:
- **QV1-cal.** As QV1, but mu = mu0[eps_inf] − rho_Lambda,obs/q0. This is a calibrated offset with relative size
  7e-121, so the calibration fixes the residual, and this must be said wherever it is used.
- **QV2 candidates.** Any extension needed to produce the observed residual. Each is a new labelled version, and QV1's
  result is kept beside it.

## Stage-2 tracks

- **scenario**: the complete solved example of QV1, including its gravitational response.
- **bound**: the class-wide relation for any q-sector remnant, and what screening would do (light q-mode as
  disappearing-Lambda f(R)).
- **vacua**: QV1 with more than one matter vacuum (metastability, degeneracy, the multiple-point principle) and the
  relation it implies between measurable Standard Model quantities.
- **dissipation**: can mu be formed dynamically (KV dissipative variants) rather than prepared?

Each track has a separate verifier that writes its own code. A relation agent assembles the gate and the inputs, and
two reviewers audit everything.

## Additional rules for stage 2

- Write only in tracks/<key>/ (verifiers: tracks/v_<key>/).
- Every number that compares with observation is retrospective. A comparison counts as held out only if the data do
  not yet exist; name the future measurement.
- A failed retrospective comparison is reported as a failure of that version. It is not repaired inside the same
  version.
- Approximations need an error estimate or a regime statement, validated numerically where feasible (for example, the
  adiabatic reduction against full integration at compressed hierarchies).
