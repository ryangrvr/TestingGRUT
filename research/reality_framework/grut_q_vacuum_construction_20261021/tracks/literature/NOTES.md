# Track `literature`: source-faithful extraction of Klinkhamer-Volovik q-theory

21 October 2026. This track belongs to the construction pass in `../../BRIEF.md`. It only extracts and translates the
sources, and it evaluates only where accuracy requires it. The construction tracks do the computing.

Labels: [D] derived here, [C] computed here (code and results in this directory), [Q] stated by the cited source
(paraphrased, with the source's equation numbers), [E] estimate, [A] argued, [R] recalled. No paper text or PDF is
stored in the repository. Downloads (arXiv TeX sources and PDFs, and pdftotext output) are in the scratchpad at
`/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/`.

Files:
- `code/kv_checks.py` contains sympy and numerical re-derivations of the source statements quoted below.
- `results/kv_checks.json` is that script's output.
- `results/sources.json` lists every source with its identifier, read level and role.

Equation numbers are the ones printed in the arXiv versions read: 0711.3170v6, 0806.2805v7 (TeX v6), 0905.1919v6,
0907.4887v6, gr-qc/0604062. I mapped them from the PDF text and checked them against the TeX labels.

---------------------------------------------------------------------------------------------------------------------

## 0. Summary of findings (details in §§1-8)

1. **The equations are transcribed correctly** [D]. The brief's (q), (g), (m) and rho_V = eps - mu q agree term by term
   with KV 0806.2805 eqs (2.6)-(2.7) and with 0905.1919 eqs (2.1)-(2.2), once three conversions are applied: Weinberg's
   curvature sign, K_brief = 2 K_KV(0905/0907) = G^-1/(8 pi), and the matter-Lagrangian sign. Section 2 checks this
   line by line, including the Bianchi identity. The equations contain no discrepancy. The mismatches are notational:
   the meaning of the letter K differs between papers, and 0905.1919 does not state its curvature convention.
2. **mu is fixed by a choice, and KV say so** [Q]. In every paper that has dynamics (0806.2805 onward), mu is a
   *global* integration constant. KV *choose* mu = mu0 so that the asymptotic matter-free state is Minkowski. They
   justify the choice by "self-sustained vacuum at zero external pressure" (0905.1919 §II, 0907.4887 eq (8b) and
   endnote) or by "memory of an original equilibrium" (0806.2805 §V.A). Later they concede that this moves the
   fine-tuning from Lambda to mu (0907.4887 §3; 1601.00601 §1; 1601.04676 §§I, II.B; 2207.02826 §I).
   **The eps in the condition** is the full effective eps(q), into which every constant matter term (the zero-point and
   condensate energy of the T = 0 Standard Model vacuum) has been absorbed (0806.2805 §II; 0907.4887 eqs (2)-(3)). It
   is evaluated in the **final, zero-temperature, post-transition Minkowski vacuum** (1102.3152 §6, which says so
   explicitly). The source does **not** use the pre-transition or the bare eps.
3. **The source contains an internal tension** [D, Q]. The static, thermodynamic papers (0711.3170 §II.C; Volovik
   2006 §7.1-7.2) let the zero-pressure condition be re-imposed after each phase transition, so mu changes. The
   dynamical equations conserve mu. At fixed mu, a constant shift Delta of eps gravitates in full: q does not move and
   rho_V -> rho_V + Delta. Readjustment would need Delta mu = Delta/q0 [D]. KV's dynamical papers accept the
   consequence: intermediate (pre-transition) epochs carry positive rho_V "plateaux" (1102.3152 §6). In the dynamical
   theory a remnant cosmological constant is equivalent to mu != mu0, with Lambda ~= -q0 (mu - mu0) [Q 0807.3896 eq
   (24); Volovik 2006 App. eq (37); D].
4. **The relaxation result of 0806.2805** [Q; C]. The oscillation is analytic (two-timing in 1/tau) and the numerics
   confirm it. The model uses eps = (1/2chi)(-F^2/F0^2 + F^4/3F0^4), G^-1 = s|F| and mu = mu0. The fast oscillation
   has omega^2 = 2 in microscopic units. The averages are <h> = 2/(3 tau) and <r_V> = (2/9) d_M^2 / tau^2. In
   dimensionful form that is **<rho_V> = (2/3) d_M^2 M^2 / t^2 [D]**, which is half the critical density of a matter-era
   universe of the same age [Q (5.17), D]. The effective equation of state is w = 0, a CDM-like a ~ t^{2/3}. I
   integrated their ODEs: tau<h> = 0.6655 and tau^2<r_V> = 0.2206 at tau 1800-3000, against KV's 2/3 and 2/9 [C].
5. **The electroweak kick** [Q; D]. The classical part is derived inside the model. For a slow kick the amplitude is
   rho_V0 = (9/8) kappa_M^2 H^4 K'(q0)^2/eps''(q0) [D], which reproduces KV's (1/8) kappa^2 h^4 and Klinkhamer's
   (9/8) kappa^2 h^4. The **remnant** E_ew^8/E_Planck^4 is not derived from the action. It needs a phenomenological
   relaxation equation, 0905.1919 eq (4.1), which Klinkhamer later called hypothetical and disconnected from q-theory
   (1001.1939 §I). It also needs an assumed kappa_M ~ 1 and a dimensional t_ew = E_P/E_ew^2. The "8" is 2 x 4.
   rho_V is quadratic in delta q, delta q is proportional to R, R is proportional to kappa H^2, and H is proportional
   to T^2/M. The numbers: (E_ew^2/E_P)^4 with E_ew = 1 TeV and E_P = 2.44e18 GeV gives (0.41 meV)^4, about 1/890 of
   rho_Lambda,obs [C]. KV write it as ~(1 meV)^4 by rounding E_P to 1e18 GeV. The follow-up papers *calibrate* E_ew
   to the observed Lambda: E_ew = 3-9 TeV (1001.1939 eq (5.2)) and E_ew = 4.7 TeV (1101.1281) [Q].
6. **Credit.** The scaling rho_Lambda ~ (M_EW^2/M_Pl)^4 is due to Arkani-Hamed, Hall, Kolda and Murayama,
   astro-ph/0005111 eq (3). They assumed an unknown mechanism that cancels the vacuum energy at the global minimum. KV
   0905.1919 present q-theory as that mechanism and credit AHHKM explicitly [Q]. GRUT's "closed whole has no bath, so
   its vacuum should not gravitate" is close to Volovik 2006 §2.2 eq (4): the natural vacuum energy is zero if the
   Universe does not interact with an environment [Q].
7. **Critiques.** I found no independent published critique aimed at the E_ew^8/E_P^4 remnant (§7). The relevant
   critiques are of three kinds:
   - Dolgov-type vector compensation, which KV adopt for their "Minkowski attractor", ruins Newtonian gravity
     (Rubakov-Tinyakov, hep-ph/9906239).
   - Global integration constants that cancel Lambda are a fine tuning, and the tuning must be redone loop by loop
     (D'Amico-Kaloper-Padilla, 1705.08950, who name q-theory as an analogue).
   - KV's own concessions, listed in §7.

---------------------------------------------------------------------------------------------------------------------

## 1. Realizations of q, the actions, and the q dependence of G or K

Conventions in the sources:
- The 4-form papers from 0806.2805 onward use signature (-+++) and **Weinberg's curvature sign**, so the action
  carries an overall minus sign: S = -∫sqrt(-g)[...].
- 0711.3170 uses signature (+---) with Einstein equation (21), (1/8piG)(R^mn - g^mn R/2) = -T_M^mn - Lambda g^mn.

| Realization | Paper, eq | Action (as written, schematic) | G or K dependence on q | Role of mu |
|---|---|---|---|---|
| Symmetric conserved tensor q^mn = q g^mn, nabla_m q^mn = 0 | 0711.3170 eq (1), §II.A | thermodynamic: E = ∫d^3r eps(Psi, Q/V), with fixed Q = qV and free V (eqs 2-4) | none (static) | mu = d eps/dq at equilibrium, a thermodynamic derivative (§II.B) |
| q^mn = <t_micro^mn> (vacuum expectation of the microscopic energy-momentum tensor) | 0711.3170 eq (49), §V.B | heuristic | none | |
| 4-form F = dA, F_klmn = q e_klmn sqrt(-g), q^2 = -(1/24) F F | 0711.3170 eqs (50)-(54), §V.C | S = ∫sqrt(-g)[R/16piG - eps(F)] (51) | constant G | d eps/dF is constant by the Maxwell equation (53) and plays mu (§V.C remark 1) |
| same 4-form with variable G | 0806.2805 eqs (2.1a-c) | S = -∫sqrt(-g)[R/(16pi G(F)) + eps(F) + L^M(psi)] | G^-1 = s|F|, k(f) = s f (4.14); motivated by 1/G ~ E_P(F)^2 ~ |F| (4.13); sigma = 3s/8pi set to 1 | integration constant of (2.5)-(2.6); effective action (2.8) has eps - mu Phi |
| several 4-forms F^(a) | 0807.3896 eqs (1)-(13) | S = -∫sqrt(-g)[K(F^(a)) R + eps(F^(a), phi) + (1/2)(dphi)^2] (2) | general K(F^(a)), K > 0 required (15) | n integration constants mu^(a) (8); rho_V = 0 fixes only one combination, so vacua can coexist (§3) |
| 4-form, K notation | 0905.1919 eqs (2.1)-(2.2); 0907.4887 eqs (1a-c), (5a,b) | S = -∫sqrt(-g)[K(q) R + eps(q) + L_SM] | 0905.1919: K linear, q0 K'(q0) = K(q0) = 1/(16pi G_N) = E_P^2/2 (below (2.5) and (2.6)); 0907.4887: K general, and K' is "not essential" (§3 end) | integration constant; mu0 fixed by (2.3b) or (8b) |
| 4-form, quadratic near q0 | 1001.1939 eqs (2.2), (3.9) | +∫sqrt(-g)[K R + eps + L^M] (overall sign differs; this does not change the field equations) | K = q/2, rho_V = (q - q0)^2/2, q0 = E_P^2 (reduced) | integration constant; the q-equation is modified by hand, (3.7a) |
| 4-form, temperature-dependent coupling | 1101.1281 eqs (2.1)-(2.2e) | ∫sqrt(-g)[K_T(q) R + eps_V(q) + L^M] | K_T = q/2 for T > T_c,K and q0/2 below, a step function (2.2e) | integration constant |
| 4-form plus a q kinetic term | 1609.03533 eq (2.1a); 1612.02326 eq (1a); 1612.04235 eq (2.1a) | adds (1/8) K(q)(nabla q^2)^2 or (1/2) C(q)(nabla q)^2 | G(q) general (1609); G_N constant (1612) | integration constant (4); propagating mode with M^2 = 1/(chi0 K q0^4) (1609 eq (3.4)) |
| 4-form, constant G, dissipative source | 1601.00601 eqs (2.1), (3.1)-(3.2), (4.1); 1601.04676 eqs (2), (7), (11), (13) | as 0806.2805 | 1601.00601: G = G_N; 1601.04676: G^-1 = sqrt(q^2)/(G_N|q0|) | mu_eff no longer constant: d mu_eff/dtau = s~ (1601.04676 eq (21)) |
| Aether or vector field: u^m_n = nabla_n u^m = u delta^m_n | 0711.3170 eqs (58)-(64), §V.D; 0907.4887 §3, eq (14) | S = ∫sqrt(-g)[R/16piG - eps(u^m_n)]; eps expanded with a "bare cosmological constant" called **K** (59) | constant G | **no integration constant** (0907.4887 §3); Dolgov-type attractor |
| Minimal 4-volume q = q^ l^4 | 0711.3170 §V.E | heuristic | | |
| Number density n of "spacetime atoms": q = n/sqrt(-g) (4D brane) | 1604.06060 eqs (11)-(14); 2207.03453 eqs (1)-(4) | S = -∫sqrt(-g)[R/16piG + eps(n/sqrt(-g)) + L^M] + mu ∫n | G_N | **mu is a constant in the action** (2207.03453, below eq (4)), not an integration constant |
| Metric or tetrad determinant (spacetime crystal, conserved lattice points) | 2207.02826 §§II-IV; 1812.07046 | -mu ∫n with n ~ M^4 sqrt(-g); unit-Jacobian diffeomorphisms | G_N | mu must be fine-tuned (§III) unless a nonstandard matter coupling is added (§IV) |
| Gluon condensate of QCD | 0811.4347 §II; 0904.3276 | q from <G^a_mn G^a_mn> (2.4); nonanalytic h~ ~ eta |R|^{1/2} |q|^{3/4} (0904.3276) | | Lambda from mu != mu0 induced by expansion: rho_vac ~ f |H| Lambda_QCD^3 (0811.4347 eq (6.3)) |
| Scalar "vacuum density" n with a kinetic term | Volovik 2006 App., eqs (26)-(30) | S_Lambda = -∫sqrt(-g)(eps(n) - mu n + gamma(n)(dn)^2); S_E with G(n) | dG^-1/dn = beta G^-1/n | mu is a Lagrange multiplier. The appendix adds that, in general, an equation letting mu relax is also needed |

Notation collisions to watch [D]:
- K means the gravitational coupling in 0905.1919, 0907.4887, 1001.1939 and 1101.1281.
- K means the bare cosmological constant in 0711.3170 eq (59).
- K means the q-kinetic coefficient in 1609.03533 and 1612.02326.
- q has mass dimension 2 in the 4-form realization and 4 in the brane realization, where mu is dimensionless
  (2207.03453 §II).
- E_Planck is non-reduced (1.22e19 GeV) in 0711.3170, 0806.2805, 0807.3896, 1601.00601 and 1612.02326.
- E_Planck is reduced (2.44e18 GeV) in 0905.1919, 0907.4887, 1001.1939, 1101.1281 and 1102.3152.

---------------------------------------------------------------------------------------------------------------------

## 2. Field equations: the KV forms against the brief's transcription, line by line

### 2.1 Conversion rules [D]

- Curvature: R_W (Weinberg, used by KV) = -R_MTW (brief), and G^W_mn = -G^MTW_mn. The metric signature is the same,
  (-+++). Evidence: KV 0806.2805 §II and 0907.4887 §2 cite Weinberg's conventions. Their de Sitter relation is R = -12
  H^2 (0806.2805 eq (3.1)), and their FRW Ricci scalar is R = -6(Hdot + 2H^2) (1601.00601 below (2.3)).
- Coupling: (1/2) K_brief = 1/(16 pi G) gives **K_brief = G^-1/(8 pi) = 2 K_KV(0905, 0907, 1001)**.
- Matter: L^M in the KV action is minus the usual Lagrangian. For example, 0807.3896 eq (2) has -∫sqrt(-g)(1/2)(dphi)^2,
  which is the correct kinetic sign in (-+++). T^M_mn is the physical tensor (T_00 = +rho) in both conventions.
- 0905.1919 states only the signature. Its curvature convention follows from the companion papers. It is also the only
  reading in which matter attracts: with MTW curvature, the signs of T and rho_V in its eq (2.2a) would be unphysical
  [A].

### 2.2 The q-equation

- KV 0806.2805 eq (2.6): d eps/dF + (R_W/16 pi) dG^-1/dF = mu.
- KV 0905.1919 eq (2.1) and 0907.4887 eq (5a): eps' + R_W K_KV' = mu.
- Translated with R_W = -R and K_brief' = (dG^-1/dF)/(8 pi) = 2 K_KV': **eps'(q) - (1/2) K'(q) R = mu**. This equals
  the brief's (q) [D]. mu arises by integrating d_nu(...) = 0 (0806.2805 eqs (2.2), (2.5)). The equation holds for
  every nu, so mu is a single constant over the whole connected spacetime: a global constant, not a field [D].

### 2.3 The Einstein equation

- KV 0806.2805 eq (2.7): (1/8piG) G^W_mn + (1/8pi)(nabla_m nabla_n - g_mn Box) G^-1 - (eps - mu F) g_mn + T^M_mn = 0.
- With G^W = -G and G^-1/(8pi) = K: -K G_mn + (nabla_m nabla_n - g_mn Box) K - rho_V g_mn + T_mn = 0. That is,
  **K G_mn = T_mn - (eps - mu q) g_mn - (g_mn Box - nabla_m nabla_n) K**, which equals the brief's (g) [D].
- 0905.1919 eq (2.2a), 0907.4887 eq (5b) and 1001.1939 eq (2.3) give the same result with K_brief = 2 K_KV [D].
- Before mu is used, the raw equation (0806.2805 eq (2.3)) has a term (1/16pi) F (dG^-1/dF) R_W g_mn. That term comes
  from the metric dependence of F. Eliminating it with (2.6) turns eps - F eps' into eps - mu F [Q, D].

### 2.4 The gravitating vacuum energy

- rho_V = eps - mu q [Q 0806.2805 eq (4.4); 0905.1919 eq (2.2b); 0907.4887 eq (6)].
- In equilibrium (R = 0), mu = eps'(q), so rho_V = eps - q eps' = eps~ [Q 0711.3170 eqs (11), (57); 0806.2805 eq
  (2.4)]. Away from R = 0 the two expressions differ by (1/2) q K' R.
- 1601.04676 eq (18a) writes the off-equilibrium expression explicitly: r~_V = eps - f eps' + f(hdot + 2h^2) [Q].
- The brief's "(in flat equilibrium: eps - q eps')" is correct as stated [D].

### 2.5 Matter conservation, the Bianchi identity, and vacuum-energy balance [D]

- The matter action does not depend on q (0806.2805 §II; 0907.4887 §2), so nabla^m T^M_mn = 0 [Q 0806.2805 eq (4.6)].
- Taking the divergence of the brief's (g), with nabla^m nabla_m nabla_n K = nabla_n Box K + R_nm nabla^m K, gives
  nabla_n rho_V = (1/2) R nabla_n K.
- Multiplying (q) by nabla_n q gives (eps' - mu) nabla_n q = (1/2) K' R nabla_n q = (1/2) R nabla_n K, which is the same
  relation. The system is consistent, and the Bianchi identity holds identically once (q) and (m) hold [D].
- In flat FRW, R = 6(Hdot + 2H^2), so rho_V-dot = 3 K-dot (Hdot + 2H^2). This equals KV 0806.2805 eq (4.5), (3/8pi)
  (dG^-1/dt)(Hdot + 2H^2) [D].
- The 00 component is 3K H^2 = rho_m + rho_V - 3H K-dot. This equals KV eq (4.2a) [D; checked symbolically, `A_FRW00`
  in `results/kv_checks.json`].

### 2.6 Verdict on the brief's transcription

**No discrepancy in (q), (g), (m), rho_V or the equilibrium conditions** [D]. The brief's equilibrium conditions,
eps'(q0) = mu0 and eps(q0) - mu0 q0 = 0, are 0905.1919 eqs (2.3a,b) and 0907.4887 eqs (8a,b) [Q]. The combined form
is eps - q eps' = 0 at q0, with mu0 derived (0907.4887 eqs (9)-(10)). The stability condition is
chi0^-1 = q0^2 eps''(q0) > 0 (0711.3170 eq (13); 0905.1919 eq (2.4); 0907.4887 eq (11)).

Four points of precision rather than error:
1. The brief's statement "mu: a spacetime constant (integration constant)" holds **only for the 4-form realization**.
   In the 4D brane realization mu is a parameter of the action (2207.03453 §II). In the aether realization there is
   no mu at all (0907.4887 §3).
2. KV allow no direct q dependence in L_matter in any paper read. The brief's S_m[psi, g] matches.
3. The constant matter vacuum energy is **absorbed into eps** (0806.2805 §II, the sentence after (2.1c); 0907.4887 eqs
   (2)-(3)). Physically, "eps" is therefore the full effective vacuum energy function, not a bare microscopic term.
4. KV's dynamics papers use one of three closures, and the conclusions depend on which:
   - linear K(q) (0806.2805, 0905.1919);
   - constant G (1601.00601, 1612.02326);
   - a T-dependent step in K (1101.1281).

---------------------------------------------------------------------------------------------------------------------

## 3. How mu is fixed (the central item)

### 3.1 What the equations alone imply [D]

1. mu is one global constant (§2.2). The equations hold for any value of mu.
2. **mu != mu0 gives de Sitter, not Minkowski**:
   - 0806.2805 §III eqs (3.2a,b), (3.7), §V.A, §V.D item (iv) and Fig. 6 [Q].
   - 0907.4887 §3 [Q]; 1601.04676 §II.B [Q].
   - Linearizing about (q0, mu0) gives rho_V ~= -q0 (mu - mu0) + O(delta mu^2) [D]. This agrees with 0807.3896 eqs
     (23)-(24), Lambda = -Sum_a q^(a) delta mu^(a) [Q], with Volovik 2006 App. eq (37), Lambda = -n0 delta mu [Q], and
     with 0806.2805 eq (3.7), delta mu = -3H^2/(8 pi G F0), so that rho_V = 3K H^2 = -q0 delta mu [D].
   - Sign check: KV's Fig. 6 uses delta u = -1/40000, and a de Sitter (positive Lambda) stage follows [Q].
3. **At fixed mu, a constant shift eps -> eps + Delta is not compensated** [D]. The q-equation contains only eps', so
   the equilibrium q is unchanged, and rho_V -> rho_V + Delta exactly in Minkowski. In de Sitter the change is
   (1 + O(Delta chi)).
4. The static "readjustment" of 0711.3170 eq (14), delta q0/q0 = chi0 eps_macro, re-imposes zero pressure after a
   transition. In the dynamical language it **requires mu to change by Delta mu = eps_macro/q0** (to first order;
   from mu' = eps'(q0 + delta q0) = mu0 + eps'' delta q0, with chi0 = 1/(q0^2 eps'')) [D]. Classical 4-form dynamics
   forbids that change, since mu is constant (§2.2).

### 3.2 Where each paper fixes mu, what it assumes, and which eps it refers to

| Paper, location | What is assumed | Which eps the condition uses | Status in the source |
|---|---|---|---|
| Volovik 2006 (gr-qc/0604062) §2.2, eqs (3)-(4); §3.3 and Fig. 1 (piston) | A many-body system isolated from its environment relaxes to T = 0 and p = 0. Applied to a Universe assumed not to interact with any environment, this gives Lambda_natural = eps_vac = -p_vac = 0 | the vacuum's thermodynamic potential <H - mu N>/V (eqs (8), (10)) in **whatever equilibrium vacuum** the system occupies | argued by condensed-matter analogy [A in source] |
| Volovik 2006 §7.1-7.2 | Every homogeneous equilibrium vacuum, true **or false**, has Lambda = 0. After an A->B transition the microscopic parameters re-adjust so that mu_B satisfies E_B - mu_B N = 0. Relaxation time and mechanism "not known" | **post-transition** (current) vacuum; mu changes | argued; the dynamics is explicitly not given |
| Volovik 2006 §7.3, eqs (16)-(20) | The electroweak energy -E_ew^4 is compensated by delta E_P/E_P ~ E_ew^4/E_P^4, so delta G/G ~ 1e-65 | post-transition | dimensional [Q] |
| 0711.3170 §II.A, eqs (2)-(4), (7) | Isolated portion of vacuum: P = 0 -> -eps + q eps' = 0 fixes q0 (fixed charge Q = qV, free volume V) | full eps(Psi0, q) of the current equilibrium; Psi0 included (eqs (5)-(6)) | thermodynamic argument [A]; the authors call the gravitational step heuristic while the microscopic theory is unknown (§III.A, first remark) |
| 0711.3170 §II.C, eqs (6), (14)-(16) | After a phase transition q0 -> q0' re-adjusts to zero pressure of the new state | **post-transition** eps_micro + eps_macro(Psi0) | thermodynamic; implies mu changes (§3.1 item 4) |
| 0711.3170 §III.A eq (30), §III.B eq (31) | Lambda = eps~(q0); q0 must be fixed by an extra condition, pressure equilibrium being the example given; with no environment, P = 0 so Lambda = 0 | equilibrium vacuum | [A]; supported by a Lorentz-invariance argument, eq (32): a vacuum with zero 3-momentum in every frame must have total eps~ energy E~ = 0 |
| 0711.3170 §V.C (after eq (57)) | In the 4-form model q0 is "determined by condition (7)". Numerical example: eps = 1/5 + (1/2) sin^2 F gives q0 ~= 3.20479 and chi0^-1 ~= 10.1887 (reproduced [C]) | eps of the 4-form action | chosen |
| **0806.2805 §III, last paragraph** | Instead of fixing H, **fix mu**. "From now on" the physical situation is fixed mu (grand-canonical) | eps(F) **with matter constants absorbed** (§II, after (2.1c)) | stated choice |
| **0806.2805 §V.A** | The Minkowski end state needs mu = mu0. Since mu is an integration constant, a physical reason "may" exist: if the nonequilibrium start came from a large perturbation of an **initial equilibrium vacuum**, the constant "would remember" it | the same function eps(F) **before** the perturbation. This silently assumes the perturbation changes initial data, not the function eps | argued [A in source] |
| 0806.2805 §V.D item (iv), Fig. 6; footnote in §V.F | mu = mu0 + delta mu with u = -1/3 - 1/40000 gives an oscillating stage, then de Sitter. Vacuum-matter energy exchange might be modelled as a small perturbation of mu | | numerical |
| **0905.1919 §II, eqs (2.3a,b)** | "We choose" mu0 so that, without matter or perturbations, the solution is the full-equilibrium Minkowski vacuum. The CC problem "would be completely solved" if the origin of mu0 were explained. Their assumption: the Minkowski vacuum is self-sustained (isolated, P = 0), and the Gibbs-Duhem relation P = -eps + mu q then gives (2.3b) | eps of the **asymptotic matter-free Minkowski vacuum**. The EW epoch is a crossover, so the function eps does not change (§V); this is assumed | stated assumption [Q] |
| **0907.4887 §2, eqs (2)-(3), (8a,b), (9)-(10), endnote "chemical potential"** | Every constant in L_SM (including SM zero-point energies) is absorbed into eps, so L_SM[psi0, eta] = 0 at the final constant fields psi0. eps = Lambda_bare + eps_var = Lambda_SM + Lambda_UV + eps_var. Then (8a,b) give (9) for q0 and mu0 is **derived** from (10). The endnote says the vacuum "can be viewed" as self-sustained at P_ext = 0, with (8b) read as the integrated Gibbs-Duhem relation | **full effective eps of the final T = 0 Minkowski vacuum** (all zero-point and condensate energies included) | stated |
| **0907.4887 §3, first paragraph** | The fine-tuning of Lambda is replaced by fine-tuning of the integration constant mu (or of q0). mu != mu0 gives de Sitter. Remedy: the aether/vector realization has **no** integration constant (Dolgov attractor, eq (14), H -> 1/t). This is "not realistic" (they cite Rubakov-Tinyakov) | | conceded [Q] |
| **1102.3152 §6** (penultimate paragraphs) | mu takes its "natural" value mu0, the true vacuum of global equilibrium. In cosmology that is the **final Minkowski state with massive fermions at zero temperature**, and **not** the value of an intermediate false vacuum. mu0 is either fixed globally as a conserved integration constant [0711.3170, 0806.2805] or emerges from an attractor [0907.4887]. **Consequence**: before each transition there is a positive plateau of rho_V, released stepwise | **asymptotic-future, post-all-transitions, full eps** | stated, with the consequence drawn |
| 1601.00601 §1 | 0806.2805 relaxes to zero "provided" mu = mu0, "but then" the CC problem becomes "why mu0?" (refs 4, 5) | | conceded |
| 1601.00601 §§3-7, eqs (5.3)-(7.1) | A dissipative backreaction from particle production, S = (gamma/36)|H|R^2 with gamma = 1/(32 pi^2) for gravitons (5.5), lets the effective mu relax. A Planck-scale bare lambda in {-1, 0, +1} E_P^4 (6.2b) is cancelled. Minkowski is an attractor for a **finite domain** of initial data (Fig. 3). The approach is slow: h ~ (6 gamma tau)^{-1/3}, r_V ~ (6 gamma tau)^{-2/3} (7.1) | any eps (the bare constant included) | numerical, plus an analytic asymptote (App. B, not read) |
| 1601.04676 abstract, §§I, II.B, IV | Reversible dynamics reaches Minkowski only if the Big Bang began in equilibrium. Minkowski is also reached when the catastrophe is local and leaves mu at spatial infinity unchanged. With dissipation S = Gamma_q qdot^2 + Gamma_H Hdot^2 there is **no attractor**: the outcome is de Sitter or collapse, and Minkowski needs fine-tuned initial data (a separatrix). De Sitter decay is argued to be **necessary** | | numerical [Q] |
| 2207.02826 §I | The original papers did not show that mu relaxes to mu0. The CC problem is to cancel **all contributions over the whole history**, not one number | | conceded [Q] |
| 2207.03453 §§II-III | Brane realization: mu is in the action. With dissipation Gamma ~ |h| R^2, Minkowski is an attractor whenever some xi with r_V(xi) = 0 exists (14) | any eps with a zero of r_V | numerical |

### 3.3 The tension between the thermodynamic and dynamical readings [D, Q]

- **Thermodynamic reading** (Volovik 2006 §7; 0711.3170 §II.C; Volovik's later thermodynamic papers 2312.02292 and
  2605.21047, abstracts only). Zero pressure is re-established in **each** equilibrium vacuum, false or true, so mu
  follows the current vacuum. Under this reading a phase transition leaves no gravitating residue after relaxation,
  and before the transition there is none either.
- **Dynamical reading** (0806.2805 onward). mu is conserved. The only value that gives a Minkowski end state is mu0 of
  the **final** eps. A constant shift of eps then gravitates in full until the epoch where it disappears (§3.1 item
  3). KV draw exactly this conclusion as "plateaux" of positive rho_V (1102.3152 §6).
- The two readings agree only if mu can change, through dissipation (1601.00601, 1601.04676, 2207.03453) or through
  the closed-droplet picture with fixed charge Q and free volume (0711.3170 §II.A). In the latter, mu = d eps/dq
  re-equilibrates because V changes. **In FRW the 4-form conserves mu, not Q = q a^3** [D]: q tends to q0 while a
  grows.
- This is a genuine inconsistency between the narratives, not a refinement. KV acknowledge the dynamical half
  (0907.4887 §3; 1601.04676 §II.B; 2207.02826 §I). I found no KV text that reconciles the Volovik 2006 §7 claim (zero
  Lambda in every false vacuum) with conserved-mu dynamics, except by appeal to unmodelled relaxation.

### 3.4 Answer to "which eps?"

| Candidate | Used by KV? |
|---|---|
| bare (microscopic-only) eps | No. Every paper absorbs matter constants into eps (0806.2805 §II; 0907.4887 eqs (2)-(3)). In 1601.00601 the bare lambda is an *input* that dissipation cancels; mu0 is not set by it |
| full renormalized, effective eps | **Yes.** All zero-point and condensate contributions are in eps (0907.4887 eq (3); Volovik 2006 §9) |
| pre-transition eps | Only via 0806.2805 §V.A's "memory of an original equilibrium". There it is harmless because the same function eps is assumed before and after |
| post-transition eps | **Static papers only** (0711.3170 §II.C; Volovik 2006 §7), with mu changing |
| asymptotic-future (final T = 0 Minkowski) vacuum | **Dynamical papers** (0905.1919 §II; 0907.4887 §2; explicitly 1102.3152 §6) |

### 3.5 Mapping to the brief's proposed formalization (source facts only)

The brief proposes "the closed whole has no exterior, so its equilibrium vacuum carries no external pressure", that
is, mu = mu0[eps].
- This is KV's own justification, nearly verbatim: Volovik 2006 §2.2 eq (4); 0711.3170 §II.A; 0905.1919 §II; and the
  1604.06060 film, which is "freely suspended" with no external forces [Q]. It must be credited to them.
- In the dynamical theory the condition is a choice of a global integration constant. KV's explicit
  dynamical-consistent choice is mu0 of the **final** full eps (1102.3152 §6).
- Under that choice, intermediate eps-shifts gravitate, and any present-day Lambda must come from an extra mechanism
  equivalent to mu_eff != mu0 (§3.1 item 2). In KV this mechanism is dissipation or freezing (§5).

---------------------------------------------------------------------------------------------------------------------

## 4. The dynamical relaxation result of 0806.2805 (read at full level, TeX source)

### 4.1 Inputs [Q]

- eps(F) = (1/2chi)(-F^2/F0^2 + F^4/(3F0^4)) (4.7), so mu0 = -1/(3 chi F0) (4.8), and chi(F0) = chi. Verified
  symbolically: mu0, rho_V(F0) = 0 and chi(F0) = chi [C].
- G^-1(F) = s|F|, k(f) = s f (4.14); sigma = 3s/(8 pi) = 1 (§IV.C: sigma is absorbed by rescaling h and tau).
- Units (4.9)-(4.10): F = f F0; y = f - 1; mu = u/(chi F0); G^-1 = k|F0|; H = h/sqrt(chi|F0|); t = tau sqrt(chi|F0|);
  rho = r/chi.
- At u0 = -1/3: r_V = y^2/2 + 2y^3/3 + y^4/6 (5.3), verified [C].
- ODEs (4.15a-c), or equivalently (5.2a,b) with r_M eliminated. Matter has a constant w_M; the analysis assumes w_M > 0
  (5.4), with w_M = 0 treated separately.
- Microscopic scales: |F0| ~ E_P^2, chi ~ 1/E_P^4, |mu0| ~ |F0| (§IV.B), with E_P non-reduced = 1.22e19 GeV.
- **Boundary conditions** at tau = 1 (figure captions):
  - Figs. 1-2: (a, h, f, r_M) = (1, 1/2, 1/2, 1/20), w_M = 1/3.
  - Figs. 3-4: (1, 1/2, 1/2, 1/200), w_M = 0, giving d_M ~= 0.986 from r_M,inf ~= 0.0125.
  - Fig. 5: (1, 1/2, 1/2, 1/2), w_M = 0, giving d_M ~= 0.590 from r_M,inf ~= 0.290.
  - Fig. 6: u = -1/3 - 1/40000.
- KV flag behaviour near tau ~ 1 as "only indicative", because of expected quantum corrections (§V.A, §V.D).

### 4.2 Analytic results (§§V.B-C) [Q]

- Linearized: ydd + 2hdot = 0 and hdot = y (5.5), so y = y0 sin(omega tau) and h = h0 - (y0/omega) cos(omega tau),
  with **omega^2 = 2** (5.6). The period of y and h is tau0 = pi sqrt(2) ~= 4.44 (5.7). The period of r_V is half that.
- Two-timing expansion in 1/tau (5.8)-(5.13). Harmonic matching fixes l0 = 2/3 (5.12), and the zeroth harmonic fixes
  b0 = omega l0 = 2 sqrt(2)/3 for w_M != 0 (5.13).
- For w_M = 0: b0 = d_M omega l0, with d_M = sqrt(1 - (9/4) r_M,inf), where r_M ~ r_M,inf/tau^2 and r_M,inf < 4/9 is
  assumed (5.14b). r_M,inf comes from the numerics.
- Asymptotes (5.15):
  - y ~ (2/3) d_M (sqrt(2)/tau) sin(omega tau);
  - h ~ (2/(3 tau))(1 - d_M cos(omega tau));
  - **r_V ~ (4/9) d_M^2 tau^-2 sin^2(omega tau)**.
- So <r_V> = (2/9) d_M^2/tau^2, <h> = 2/(3 tau), <a> ~ tau^{2/3}, and <rho_V> ~ a^-3. The q-field itself has w = -1,
  but the oscillating (rho_V, G) system acts as **w = 0, CDM-like** (§V.E). The effective Lagrangian
  (1/2) ydot^2 - (1/2) omega^2 y^2 is "induced by gravity", and the pressure averages to zero.
- KV note the same oscillation of h as in Starobinsky 1980 [Q].

### 4.3 Numerics [Q] and the independent check [C]

- KV's items (i)-(iv) (§V.D):
  - f -> 1;
  - <a> ~ tau^{2/3} even for ultrarelativistic matter, instead of tau^{1/2};
  - the same for nonrelativistic matter at both low and high initial density;
  - with u slightly below u0, oscillation and then a de Sitter stage (Fig. 6).
- Check [C]: I integrated (5.2a,b) with the Fig. 1-2 boundary conditions (`code/kv_checks.py` §C, DOP853, rtol
  1e-11). Over tau in [1800, 3000]:
  - tau<h> = 0.6655 (KV 2/3) and tau^2 <r_V> = 0.2206 (KV 2/9 = 0.2222). The remaining 0.7% is the expected O(1/tau)
    correction.
  - w_M = 0 with the Fig. 3-4 conditions: measured r_M,inf = 0.0114 gives d_M = 0.987 (KV: 0.986).
    tau^2 <r_V> = 0.21655 against (2/9) d_M^2 = 0.21650.
  - w_M = 0 with the Fig. 5 conditions: r_M,inf = 0.2890 gives d_M = 0.591 (KV: 0.590). tau^2 <r_V> = 0.07782
    against 0.07774.
- **The asymptotic law is reproduced** [C].

### 4.4 Dimensionful late-time law [D]

- With sigma = 1: s = 8 pi/3 and |F0| = G^-1/s = 3/(8 pi G) = **3 K = 3 M^2** (brief units). Then
  **<rho_V> = (2/9) d_M^2 |F0| / t^2 = (2/3) d_M^2 M^2 / t^2** [D].
- The critical density with <H> = 2/(3t) is (4/3) M^2/t^2, so the ratio is 1/2 for w_M > 0. This reproduces KV's
  eq (5.17), including their factor d_M^2 for w_M = 0 [D].
- The source of the "right order of magnitude" is that <rho_V> **tracks the critical density** of a t^{2/3} universe of
  age t [D]. The Planck scale cancels out of the coefficient.
- **KV's own estimate** (5.16) [Q]: <rho_V(t_now)> ~ |F0|/t^2 ~ E_P^2/t^2 ~= (4e-3 eV)^4 (1e10 yr/t)^2, with E_P
  non-reduced and t_P = 5e-44 s. KV stress that a power t^-n with n != 2 would destroy the agreement.
- Recomputed [C]: E_P^2/t^2 at t = 1e10 yr is (5.0 meV)^4. That is a minor discrepancy with KV's 4 meV in an
  order-of-magnitude statement. With the full coefficient, (2/3) M^2/t^2 is (2.04 meV)^4 at 1e10 yr and (1.73 meV)^4
  = 0.36 rho_Lambda,obs at 13.787 Gyr [C]. This is illustrative only, because the model's history (a ~ t^{2/3}
  throughout) is not ours.
- **KV's own caveats** (§V.F) [Q]:
  - the 1/t^2 behaviour contradicts the observed acceleration;
  - a ~ t^{2/3} contradicts big-bang nucleosynthesis, which needs t^{1/2};
  - particle production (parametric resonance) and other processes "intervene".

### 4.5 G(t) [Q, D]

- Eq. (5.18): G^-1(t) ~ G^-1_inf (1 + c0 (t_UV/t) sin(t/t_UV)), with t_UV = sqrt(chi|F0|/2) ~ the Planck time [Q].
- From (5.15a), y = delta k/k, so c0 = (4/3) d_M [D] (KV: "of order unity").
- The amplitude today is ~ t_P/t0 ~ 1e-61 at a Planck frequency. KV judge the solar-system bounds satisfied, while
  they do not rule out surprises [Q].

### 4.6 What is analytic, what numerical, and what was later revised

- Analytic: the linear oscillation, l0 and b0 (for w_M > 0), and the asymptotic forms.
- Numerical: r_M,inf and hence d_M for w_M = 0; the transient behaviour; the mu != mu0 behaviour.
- Assumed: mu = mu0; linear G^-1(F); the quartic eps; no particle production.
- KV later replace the 1/t^2 relaxation by **exponential decay** once particle production is included, with Gamma ~
  omega ~ E_P (0905.1919 eq (1.1b); 1102.3152 §6). With the R^2-type backreaction and constant G they replace it by a
  **slow t^{-2/3} attractor** (1601.00601 eq (7.1); also for nonconstant G, §7). For that attractor the decay
  constant is gamma = 1/(32 pi^2) for gravitons (5.5).
- Extrapolated naively to t0 in Planck units, the 1601.00601 attractor would leave rho_V ~ 1e36 GeV^4 [E]. KV present
  it as an early, inflation-like phase, with "possibly" a standard FRW phase afterwards (§8), not as late cosmology.
- The **CDM-like oscillation** is revived in 1612.02326 and 1612.04235 with a q kinetic term and constant G:
  - rho_q-DM = (1/2) chi0^-1 a_xi^2 (21a), with P ~ 0 on average;
  - rho ~ a^-3 numerically (5.2a in 1612.04235);
  - frequency omega^2 = 1/(q0 chi0) ~ E_P^2 (17b).
  The amplitude, and so the DM/DE ratio, is **not calculated**. §5 says no definitive calculation exists yet.

---------------------------------------------------------------------------------------------------------------------

## 5. The electroweak kick (0905.1919 full; 0907.4887 full; follow-ups 1001.1939 and 1101.1281 at section level)

### 5.1 The derivation step by step, with status

| # | Step (0905.1919) | Equation | Status in source |
|---|---|---|---|
| 1 | Before the kick, Planck-scale oscillations have already decayed rho_V to zero, by dissipative decay at Gamma ~ omega ~ E_P (refs Starobinsky 1980; Vilenkin 1985) | (1.1b); §IV | **assumed** (needs mu = mu0 or dissipation) |
| 2 | Classical equations of 0806.2805 with linear K (q0 K' = K), quadratic eps near q0, and mu = mu0 | (2.1)-(2.7) | derived within the model; mu0 **chosen** (§3.2) |
| 3 | kappa_M = 4 - 3(1 + w_M) = 1 - 3w_M | (2.8) | definition |
| 4 | For kappa_M = 0 (pure radiation) there is an exact solution y = 0, h = 1/(2 tau): rho_V stays zero | (2.10) | **derived** (R = 0 for radiation) |
| 5 | Linearization in y and h^2: ydd + 3h ydot + omega^2 y = kappa_M h^2 and hdot + 2h^2 - y = 0, with omega^2 = 2 | (2.11a,b) | derived (first order) |
| 6 | Smooth kick (Delta tau_kappa >> 1/omega): drop ydot and ydd, so y = kappa_M h^2/2 | (3.1a,b) | derived (adiabatic). At t_ew, omega t_ew ~ 1e30, so this approximation is exact for practical purposes [E] |
| 7 | h(tau) = [2∫(1 - kappa/4) dtau']^-1; after the kick h ~= 1/(2(tau - tau0)) with tau0 = (1/4)∫kappa dtau | (3.2)-(3.3) | derived |
| 8 | r_V = (1/8) kappa_M^2 h^4; dimensionful rho_V ~ kappa_M^2 H^4, the precise constant is left to microphysics and taken to be O(1) | (3.4a,b) | (3.4a) derived; the O(1) constant **assumed** |
| 9 | t_ew ~ E_P/E_ew^2 with an "effective" E_ew ~ 1e3 GeV (not the crossover temperature, T_ew ~ 1e2 GeV, which §V quotes as ~300 GeV) | (1.2) | **dimensional** (radiation Friedmann) |
| 10 | kappa_M(t_ew) ~ 1. This needs massive particles dominating; TeV-scale new physics "may enhance" it | below (3.5); §III | **assumed** |
| 11 | rho_V0(t_ew) ~ H^4(t_ew) ~ t_ew^-4 ~ E_ew^8/E_P^4 | (3.5) | dimensional |
| 12 | Classically rho_V returns to zero after the kick | text after (1.3); §IV | derived and acknowledged |
| 13 | Phenomenological relaxation rho_V-dot = -Gamma(t)(rho_V - rho_V0(t)) | **(4.1)** | **assumed** ("rough estimate", awaiting a theory of irreversible vacuum dynamics) |
| 14 | Gamma peaks at t_ew with Gamma(t_ew) ~ 1/t_ew. Argument: Zel'dovich-Starobinsky n-dot ~ R^2, which vanishes in pure radiation; produced quanta have energy ~ 1/t_ew ~ meV, so only gravitons and neutrinos qualify | (4.3) and text | **argued** heuristically |
| 15 | Solution (4.4): lim rho_V ~ rho_V0(t_ew), a "time-lag" effect, since the response time is comparable to the kick width | (4.4)-(4.5) | estimate |
| 16 | Later kicks (QCD at H^4(t_QCD); matter era below H^4(t_rec)) are negligible. A crossover far above a TeV would give too large a remnant, which would argue against a GUT scale | §IV end | estimate |
| 17 | Lambda ~ (E_ew^2/E_P)^4 ~ (1e-3 eV)^4 with E_P ~ 1e18 GeV (abstract; (1.4)), against observed ~(2 meV)^4. Heavier TeV relics could raise E_ew (§V) | (1.4); §V | numerical illustration |

### 5.2 The general coefficient of the classical kick [D, C]

- Use the brief's (q) at mu = mu0 with K ~= K0 and the trace of (g): -K R = T^m - 4 rho_V, so R ~= (rho_m - 3P_m)/K
  = 3 kappa_M H^2.
- Then delta q = K'(q0) R/(2 eps''(q0)) and rho_V = (1/2) eps'' delta q^2 = **(9/8) kappa_M^2 H^4 K'(q0)^2/eps''(q0)**
  [D; sympy `B_kick`].
- In the 0806.2805 units (K' = 1/3, eps'' = 1/(chi F0^2)) this gives r_V = (1/8) kappa^2 h^4, which equals KV (3.4a)
  [C]. With K = q (brief) and eps'' = 1 it gives (9/8) kappa^2 h^4, which equals 1001.1939 (3.18b) [C].
- With linear K: K'^2/eps'' = K0^2 chi = M^4 chi. So the O(1) constant KV assume is **(9/8) kappa_M^2 M^4 chi**.
  It is O(1) only if chi^-1 ~ M^4 (reduced). If chi^-1 = E_P^4 (non-reduced, as in 0806.2805 §IV.B), the constant
  is smaller by (8 pi)^2 ~= 632 [D].
- Also: with t_ew = E_P/E_ew^2 and H = 1/(2t), H(t_ew)^4 = E_ew^8/(16 E_P^4) [D].

### 5.3 What sets the power 8, and what would change it [D]

rho_V is proportional to (delta q)^2. The factor delta q is proportional to K' R/eps'', R is proportional to
kappa_M H^2, and H is proportional to T^2/M in the radiation era. That gives
rho_V ~ kappa^2 (K'^2/eps'') T^8/M^4 at the freezing epoch.

The exponent 8 = 2 (quadratic minimum at mu = mu0) x 2 (R is quadratic in H) x 2 (H is quadratic in T). It changes
if any of the following holds:
1. **mu != mu0.** A linear term -q0 delta mu dominates, and Lambda is a free constant.
2. **K' = 0.** Constant G (the setting of 1601.00601 and 1612.02326) gives no curvature coupling. The q-equation
   reduces to eps' = mu, and there is no kick.
3. **The freeze is not at H(t_ew).** The remnant is H^4 at whatever epoch Gamma ~ H. 1001.1939 §III.E shows a finite
   remnant needs both 0 < zeta and zeta < infinity: Gamma >> H tracks kappa back to 0, and Gamma << H leaves
   x = 0.
4. **The kick is not from matter-trace anomalies.** For example, KV's QCD variant has rho ~ |H| Lambda_QCD^3, giving
   Lambda ~ E_QCD^6/E_P^2 (0811.4347 eq (6.3); 1102.3152 §1).

Changes in chi or kappa_M change only the prefactor, not the exponent.

### 5.4 Numbers [C; inputs E_ew = 1e3 GeV, M = 2.435e18 GeV, rho_obs = 2.52e-47 GeV^4]

- (E_ew^2/E_P)^4 with E_P = 2.44e18 GeV gives 2.82e-50 GeV^4 = (0.41 meV)^4, which is **rho_obs/893**. KV's
  "(1e-3 eV)^4" uses E_P rounded to 1e18 GeV (rho_obs/25).
- With unit coefficient, (E_ew^2/M)^4 = rho_obs requires **E_ew = 2.34 TeV**.
- Radiation era: H = 3.42 T^2/M at g* = 106.75. **H^4 = rho_obs at T = 1.26 TeV**. At the SM crossover temperature
  (~160 GeV [R]), H^4/rho_obs = 6.5e-8.
- Follow-ups, **calibrated to the observed Lambda** [Q]:
  - 1001.1939 eq (5.1): Lambda_theory = r_V^num E_ew^8/E_P^4. Eq (5.2) gives E_ew = 3.2 TeV for r_V^num = 0.051 and
    Lambda_exp = (2.0 meV)^4.
  - Range: E_ew = 3-9 TeV, (5.6). Lower bound: E_ew >= 2 TeV if r_V^num <= 1, (5.3).
  - 1101.1281: r_V(inf) ~ 2.4e-3, so E_ew ~ 4.7 TeV (§IV).
  - My recomputation with rho_obs = (2.24 meV)^4: 3.39, 4.97, 9.04 and 2.34 TeV for r_V^num = 0.051, 2.4e-3, 2e-5 and
    1 [C].
  - In these papers E_ew is a fitted parameter, **calibrated, not predicted**.
- The modifications in the follow-ups are themselves phenomenological:
  - 1001.1939 eqs (3.7a,b): a dissipation function gamma(t/t_ew) and coupling zeta. The q-equation is modified "to
    match" a standard Lambda-FRW at late times. The ODEs are deliberately singular after t_freeze.
  - 1101.1281 eq (2.2e): a discontinuous step in K_T(q) at T_c,K ~ E_ew, described as the single key ingredient still
    to be derived (§IV).
  - Both papers state that a frozen nonzero rho_V is **inconsistent with the unmodified q-equation** in a de Sitter
    end state (1001.1939 §III.E, fourth remark; 1101.1281 §I: "potential mismatch"). This is consistent with §3.1
    item 2: a lasting Lambda means mu_eff != mu0 [D].

### 5.5 Arkani-Hamed, Hall, Kolda and Murayama, astro-ph/0005111 (PRL 85, 4434 (2000); full read) [Q]

- Two scales are assumed, M_EW ~ 1 TeV and the reduced M_Pl ~ 1e18 GeV.
- WIMP relic density: rho_M ~ (M_EW^2/M_Pl) T^3 (1). Matter-radiation equality: T_eq ~ M_EW^2/M_Pl (2). A "triple
  coincidence" therefore requires **rho_Lambda ~ (M_EW^2/M_Pl)^4** (3).
- Models:
  - hidden-sector SUSY QCD whose quark and gluino masses are ~ M_EW^2/M_Pl, giving false vacua split by
    ~ (M_EW^2/M_Pl)^4 (5)-(9);
  - a brane construction whose n = 3 term gives the same combination (10);
  - a quintessence variant with m_phi^2 ~ M_EW^8/M_Pl^6.
- AHHKM explicitly **assume an unknown mechanism cancels the vacuum energy at the global minimum**.
- Order-one factors matter: 1/(pi alpha^2) ~ 1e3 moves T_eq to ~1 eV.
- **Credit**: the TeV^2/M scaling and the coincidence argument are AHHKM's. KV 0905.1919 (abstract, §I, §V) and
  Klinkhamer 1001.1939 §I credit them, and KV claim to supply AHHKM's missing cancellation mechanism [Q].

---------------------------------------------------------------------------------------------------------------------

## 6. Observable consequences KV claim, with numbers [Q unless marked]

| Claim | Paper, location | Number or form | Nature |
|---|---|---|---|
| Lambda = 0 in the perfect vacuum; Lambda != 0 only from Lorentz-violating perturbations (matter, boundaries) | 0711.3170 §§III-IV | Nongravitating matter: eps~_vac = P_M (40), (42). Droplet: eps~ = P_sigma + P_M (47) | thermodynamic |
| Vacuum response to matter | Volovik 2006 §8.1-8.2 | eps_vac = w eps_matter (22)-(23) (no gravity); Einstein static universe eps_vac = (1/2) eps_matter (1 + 3w) (25); App. (31) | thermodynamic and GR |
| Natural macroscopic scales of Lambda | Volovik 2006 §8.4 | ~ H^2/G, eps_matter, T^4, Hdot/G, 1/(G R^2) | dimensional |
| Change of constants from a transition | 0711.3170 §II.C | delta q0/q0 ~ 1e-68 (16); delta alpha/alpha ~ 1e-68 alpha (18); delta G/G ~ E_ew^2/E_UV^2 ~ 1e-34 (20) | dimensional |
| Same | Volovik 2006 §7.3 | delta G/G ~ E_ew^4/E_P^4 ~ 1e-65 (20) | dimensional. **Differs from 0711.3170 eq (20)** (1e-34), which adds a low-energy E_ew^2 term [D] |
| Relaxing vacuum energy is the right order today | 0806.2805 (5.16)-(5.17) | (4e-3 eV)^4 (1e10 yr/t)^2; rho_V/rho_c -> 1/2 (times d_M^2 for w = 0) | analytic and numerical; the source says it conflicts with acceleration and BBN |
| Oscillating vacuum energy acts as CDM | 0806.2805 §V.E; 1612.02326; 1612.04235 | rho ~ a^-3, P ~ 0; Planck frequency, so **direct detection fails** (1612.02326 §5) | analytic and numerical; amount not computed |
| Planck-scale oscillation of G | 0806.2805 (5.18) | relative amplitude ~ (4/3) d_M t_UV/t | analytic |
| Effective f(R) near equilibrium | 0807.3896 (20)-(22) | f~ = -R + 16 pi G chi~ R^2 - 2 lambda (Weinberg sign); scalaron M^2 = 1/(96 pi G chi~), chi~ = chi0 (F dK/dF)^2, of Planck scale; lambda = -8 pi G Sum q delta mu (23) | derived (near equilibrium) |
| q-waves | 1609.03533 (3.4) | M^2 = 1/(chi0 K q0^4) ~ E_P^2 (K = kinetic coefficient) | derived |
| q-balls | 1609.03533 §4 | stable droplet of physical vacuum in the q = 0 "empty vacuum"; a present-size Lambda needs R ~ E_P/H^2, far beyond the horizon | heuristic |
| Remnant Lambda from the EW kick | 0905.1919 (1.4), (4.5) | ~ E_ew^8/E_P^4 ~ (1e-3 eV)^4 for E_ew ~ 1 TeV; the present epoch is observationally the same as LambdaCDM (w = -1) | dimensional plus assumed relaxation |
| New TeV physics required | 0905.1919 §V; 1001.1939 (5.2)-(5.6); 1101.1281 §IV | E_ew = 3-9 TeV (>= 2 TeV); 4.7 TeV | **calibrated to Lambda_obs** |
| No grand unification below the Planck scale | 0905.1919 §IV end | a crossover at T >> TeV would leave a much larger remnant | argued |
| Neutrino plateau | 1102.3152 §5, eqs (10)-(11) | rho_vac = c_nu [tr(M^2) det(M^2)]^{1/3} = (2.3 meV)^4. With c_nu = 1: normal hierarchy (2.8e-6, 8.8, 49.0) meV, inverted (49.0, 49.8, 1.8e-7) meV; KATRIN (0.2 eV design) would not see them. If KATRIN finds >= 0.4 eV, the neutrino plateau must already be released | heuristic; **calibrated to rho_vac,obs** |
| QCD alternative | 0811.4347 (6.3)-(6.4); 0904.3276 (3.17b) | Lambda = f H Lambda_QCD^3, giving Lambda ~ K_QCD^3/E_P^2; w_X(t_p) ~= -0.66, tending to -1 from above (one model) | heuristic plus numerical; 1102.3152 §6 says it may already be ruled out |
| Minkowski attractor with dissipation | 1601.00601 (7.1)-(7.2); 2207.03453 (18) | h ~ (6 gamma tau)^{-1/3}, r_V ~ tau^{-2/3}; inflation-like growth of the particle horizon | numerical plus analytic asymptote |
| Big bang as a topological transition between q = -q0 and +q0 | 2111.07962 (abstract only) | kink solution; conformal at q = 0 | not read beyond the abstract |
| Local energy-momentum nonconservation in generalized unimodular gravity | 1604.03065 (abstract) | negligible for meV vacuum scale | abstract |

---------------------------------------------------------------------------------------------------------------------

## 7. Published critiques and independent analyses

Search method:
- arXiv API, with queries au:Klinkhamer AND abs:q-theory, au:Klinkhamer AND au:Volovik, abs:"q-theory" AND
  (vacuum OR cosmological), abs:"vacuum variable".
- INSPIRE citation lists ("refersto:recid") for 0711.3170 (121 citing records), 0806.2805 (95), 0905.1919 (30) and
  0907.4887 (44). I sampled the citation contexts of the non-KV records listed below.

1. **Rubakov and Tinyakov, hep-ph/9906239 (PRD 61, 087503 (2000))** (read §§1-3).
   - Target: Dolgov's vector and tensor compensation. KV 0907.4887 §3 uses it as their "existence proof" of a
     Minkowski attractor without an integration constant.
   - Finding: the compensating field grows to A0 ~ sqrt(eps_vac)/H, which is far above M_Pl. Newton's law becomes
     unacceptable (the effective Planck mass is too large), and gravitational waves differ completely from GR.
   - KV acknowledge this in 0907.4887 §3 (third remark) and 1601.00601 §7 [Q].
   - Responses:
     - Emelyanov and Klinkhamer, 1107.0961 and 1109.4915 (abstract plus grep of the text): modified two-vector
       models said to restore Newtonian dynamics and an asymptotic FRW.
     - Santillan and Scornavacche, 1705.04897 (abstract): an independent check that agrees G_N and gravitational
       waves behave well at late times in such vector adjusting mechanisms.
2. **D'Amico, Kaloper and Padilla, 1705.08950** (read the section around their eqs (17)-(19)).
   - A sourceless 4-form flux is a global degree of freedom, and they compare it to q-theory models. Its value is arbitrary
     and can be used to cancel the vacuum energy, but picking an integration constant to cancel Lambda is fine
     tuning.
   - In global-multiplier proposals that do not treat the CC counterterm as dynamical, the tuning must be redone
     order by order in loops.
   - The paper is not aimed at KV; it names q-theory as an analogue. The point applies directly to the
     mu0-as-integration-constant reading (§3).
3. **KV's own concessions.** These are the most specific critiques of their construction in the literature.
   - 0907.4887 §3: the fine-tuning moves to mu.
   - 1601.00601 §1 and 1601.04676 abstract and §II.B: Minkowski needs mu = mu0 or fine-tuned initial data; de Sitter
     decay is necessary.
   - 2207.02826 §I: the original papers did not show that mu relaxes; the problem is the whole history, not one number.
   - 1001.1939 §I: eq (4.1) of 0905.1919 was hypothetical and disconnected from q-theory.
   - 1101.1281 §I: a potential mismatch with GR plus Lambda; modification by hand is "unsatisfactory".
   - 0806.2805 §V.F: the 1/t^2 relaxation conflicts with acceleration and BBN.
   - 0711.3170 §III.A: the gravitational identification is heuristic.
4. **Osorio Morales and Santillan, 1808.07756** (abstract plus the gauge-condition passage). The Cauchy problem of
   q-models is well posed for a finite time for suitable data. Caveat: if the gauge condition is violated during
   evolution, the solution is not a q-theory solution. This is a technical analysis, not a critique of the vacuum
   claim.
5. **Chung and Long, 1108.5193** (abstract plus introduction). They list self-tuning mechanisms, KV 0711.3170 among
   refs 15-27, as alternatives that would remove or strongly modify the effect of vacuum energy during an EW-scale
   phase transition. That effect is testable through a shift in the dark-matter relic abundance. For
   non-first-order transitions (crossovers) the effect is suppressed, O(1/g_E). This is relevant to §3.3: the
   thermodynamic and dynamical readings differ on whether the pre-transition vacuum energy gravitates.
6. **Emelyanov, 1907.07938** (abstract, introduction, conclusion line). In lambda phi^4, a nonzero renormalized
   zero-point energy in Minkowski is in tension with the field equation at two loops. This independently supports a
   zero Minkowski vacuum energy. It is not a critique.
7. **Xie, Hsiang and Hu, 2312.09047** (abstract plus the citation context). They compute vacuum compressibility for
   dynamical spacetimes from QFT. They cite KV 2008 as motivation but do not test q-theory. Neutral.
8. **Selch, Miller and Zubkov, 2301.02914** (abstract plus introduction). Gravastar-like black-hole solutions in
   q-theory with a minimally coupled q-scalar. This is an extension, not a critique.
9. **The Weinberg no-go** (RMP 61, 1 (1989) §6) [R]. KV claim to evade it: the q-equation requires only
   d/dx(dL/dq) = 0, so dL/dq = mu with arbitrary mu (0907.4887 §2) [Q]. Point 2 is the standing objection to that
   evasion.
10. **On the E_ew^8/E_P^4 remnant specifically.** I found **no independent published critique or reanalysis**. Of the
    30 INSPIRE records citing 0905.1919, the contexts sampled are neutral list-citations: Urban-Zhitnitsky
    0906.2162, Basilakos-Plionis-Sola 0907.4555 and Poplawski 1005.0893. The only critical assessments are
    Klinkhamer's own (item 3). [Limit: not all 30 contexts were opened.]

---------------------------------------------------------------------------------------------------------------------

## 8. Claim table

Read levels:
- **full**: the whole text, including TeX source where downloaded.
- **sect**: the named sections.
- **abs**: the abstract only.

Status is the status in the source, judged by this track.

| Claim | Paper / section | Read level | Status in source |
|---|---|---|---|
| Vacuum is a self-sustained medium; P = 0 in isolation gives Lambda = 0 | Volovik 2006 §2.2 eqs (3)-(4); 0711.3170 §II.A, §III.B eq (31) | sect (Volovik §§1-3, 7-9, App. part); full (0711) | argued (thermodynamic analogy) |
| Gravitating vacuum energy is eps - q eps' (= -P), not eps | 0711.3170 eqs (11), (28)-(30); 4-form check (55)-(57) | full | derived for the 4-form and aether models; the general identification is "heuristic" |
| Readjustment after a transition: delta q0/q0 = chi0 eps_macro ~ 1e-68 | 0711.3170 eqs (14)-(16) | full | derived (static) |
| delta G/G ~ 1e-34 (0711) or ~1e-65 (Volovik 2006) | 0711.3170 eq (20); Volovik eq (20) | full / sect | dimensional |
| Field equations (q), (g), rho_V = eps - mu q | 0806.2805 eqs (2.2)-(2.8); 0905.1919 eqs (2.1)-(2.2); 0907.4887 eqs (4)-(6) | full | derived (verified here, §2) |
| mu = mu0 (Minkowski end state) | 0806.2805 §V.A; 0905.1919 §II; 0907.4887 eqs (8)-(10) | full | **assumed / chosen** |
| mu0 is the final T = 0 vacuum; intermediate plateaux gravitate | 1102.3152 §6 | sect | stated |
| mu != mu0 gives de Sitter, with Lambda = -q0 delta mu | 0806.2805 §III, Fig. 6; 0807.3896 eq (24); Volovik App. eq (37) | full / full / sect | derived plus numerical |
| Relaxation <r_V> = (2/9) d_M^2/tau^2, CDM-like | 0806.2805 §§V.B-E, eqs (5.5)-(5.15) | full | analytic (two-timing) plus numerical (d_M); reproduced [C] |
| Present-day order of magnitude (4 meV)^4 | 0806.2805 eq (5.16) | full | dimensional extrapolation; conflicts with data per source |
| G oscillation (5.18) | 0806.2805 §V.F | full | analytic |
| Effective f(R) = R + R^2 term; Planck scalaron | 0807.3896 eqs (16)-(22) | full | derived near equilibrium |
| Classical EW kick r_V = kappa^2 h^4/8 | 0905.1919 eqs (2.11), (3.1)-(3.4a) | full | derived (adiabatic, linear) |
| rho_V0(t_ew) ~ E_ew^8/E_P^4 | 0905.1919 eq (3.5) | full | dimensional (kappa ~ 1 and O(1) constant assumed) |
| Remnant freezes at rho_V0(t_ew) | 0905.1919 eqs (4.1)-(4.5) | full | **assumed** relaxation law; Gamma(t_ew) ~ 1/t_ew argued |
| E_ew = 3-9 TeV | 1001.1939 eqs (5.1)-(5.6) | sect (§§I-V; App. A not read) | numerical within phenomenological ODEs; **calibrated to Lambda_obs** |
| E_ew ~ 4.7 TeV; a step in K_T(q) is needed | 1101.1281 §§II-IV | sect | phenomenological model; calibrated |
| Neutrino plateau, masses (11) | 1102.3152 §5 | sect | heuristic; calibrated |
| Fine-tuning shifted from Lambda to mu | 0907.4887 §3; 1601.00601 §1; 1601.04676 §II.B; 2207.02826 §I | full / full / sect / sect | conceded |
| Aether/vector: no mu, Minkowski attractor H -> 1/t | 0907.4887 §3 eq (14); Dolgov 1997 | full | numerical existence proof; phenomenologically excluded in its simple form (Rubakov-Tinyakov) |
| Dissipative attractor, r_V ~ tau^{-2/3}, gamma = 1/(32 pi^2) | 1601.00601 §§5-7, App. A | full main text; App. B not read | numerical plus asymptote; gamma computed (App. A) |
| No attractor for S = Gamma_q qdot^2 + Gamma_H Hdot^2 | 1601.04676 §§III-IV | sect | numerical |
| Brane q = n/sqrt(-g); film analogy; generalized unimodular gravity | 1604.06060 | full | argued plus derived equations |
| Brane attractor | 2207.03453 §§II-III | sect | numerical |
| q-waves M^2 = 1/(chi0 K q0^4); q-balls | 1609.03533 §§2-4 | sect | derived / heuristic |
| q-oscillation DM, rho ~ a^-3, no direct detection | 1612.02326 (full); 1612.04235 §§I-V | full / sect | derived plus numerical; amount not computed |
| QCD remnant ~ K_QCD^3/E_P^2; w_X ~ -0.66 | 0811.4347 §VI, App. A; 0904.3276 §III | sect | heuristic; numerical in one model |
| Equilibrium boundary conditions plus d rho_V/dt ~ rho_M (closed FRW) | 0803.0281 §§I-III | sect | phenomenological "Keplerian" ansatz |
| Metric determinant as q; mu fine-tuned unless unimodular-type coupling | 2207.02826 §§I-III | sect | argued plus numerical (not read beyond §III) |
| Tetrad or crystal realization; IIB-matrix origin; q-bubble decay; Ostrogradsky stability | 1812.07046; 2212.00709; 1901.05938; 1704.05436 | abs | — |
| AHHKM triple coincidence, rho_Lambda ~ (M_EW^2/M_Pl)^4 | astro-ph/0005111 eqs (1)-(10) | full | dimensional plus toy models; cancellation at the global minimum assumed |

---------------------------------------------------------------------------------------------------------------------

## 9. Pointers for the construction tracks (source facts and derivations only)

1. Use the brief's equations as written. They match KV (§2). If K'(q) != 0, the vacuum sector also enters the
   expansion through -3H K-dot and through an oscillating G. KV's "CDM-like" result needs both (§4.2, §4.4).
2. "mu = mu0[eps]", in the only dynamically consistent KV form, means mu0 of the **final** full eps (§3.2, 1102.3152
   §6). Every constant shift of eps that happens at or after the epoch of interest then gravitates in full until it
   disappears (§3.1). Re-imposing zero pressure after a transition is the static papers' version. It needs a change
   Delta mu = Delta/q0 that the 4-form equations forbid (§3.1, §3.3).
3. A present-day Lambda in any KV variant equals -q0 (mu_eff - mu0). It comes from a non-conservative, phenomenological
   ingredient: Gamma in 0905.1919 (4.1); gamma and zeta in 1001.1939; the K_T step in 1101.1281. None is derived from
   an action.
4. Classical-kick coefficient (§5.2): (9/8) kappa_M^2 H^4 K'^2/eps'' = (9/8) kappa_M^2 H^4 (M^4 chi) for linear K. The
   factor M^4 chi is free in KV. A reduced-versus-non-reduced Planck choice in chi alone moves it by 632.
5. Credit lines: AHHKM for the scaling and the coincidence argument; KV and Volovik for the self-sustained-vacuum,
   zero-pressure argument, the 4-form realization with K(q), the equations, the relaxation, the kick and the remnant
   hypothesis; Dolgov, and Emelyanov-Klinkhamer, for the vector attractors.

---------------------------------------------------------------------------------------------------------------------

## 10. Limits of this extraction

- Figures in 0806.2805, 1601.00601, 1601.04676 and 1001.1939 were not viewed as images. Their captions and text were
  read.
- Not read: Volovik 2006 §§4-6; 1001.1939 App. A; 1601.00601 App. B; 1101.1281 §III.C in detail; 1609.03533 §§5-6;
  1612.04235 beyond §V; 2207.02826 beyond §III.
- Read at abstract level only (see `results/sources.json`): 1704.05436, 1812.07046, 1901.05938, 2212.00709,
  2111.07962, 1108.1995, 1204.5085, 1205.7072 (§1 read), 1604.03065, 2205.04703, 2312.02292, 2605.21047.
- The critique search covered the arXiv API and INSPIRE citation lists. Journals not indexed by INSPIRE were not
  searched. Not every citing record's context was opened.
- The equation numbers are those of the arXiv versions listed at the top. Published versions may number differently.
- KATRIN and SM crossover-temperature figures are [R], not checked here.
