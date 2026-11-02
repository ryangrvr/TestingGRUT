# Track qm_transfer: linear matter power of Q-matter + a kinetic warm vapour, against Lyman-α, satellites and lensing

2 November 2026 · Claude, decisive pass (BRIEF_DECISIVE.md of this directory; issue #13 context). Same-session AI work on
one supplied model (M1r, Q-matter branch), not outside replication. Nothing here is observational evidence, a
validation of GRUT, an exclusion of the real world or a new theory; all data comparisons are retrospective and the
probe criteria are proxies. **Every input (M1r, the two-phase vapour model, sticking, the kinetic closure) is a model
choice, not a consequence of closure.** Git read-only, nothing committed. At most 2 worker processes, foreground.
Labels: [D] derived, [C] computed with our own code, [E] estimate, [A] argued, [Q] quoted from a reviewed track,
[R] recalled/literature (arXiv numbers given; the bound values were checked by web search on 2 Nov 2026).

## Result in one paragraph

A new linear Boltzmann code (synchronous gauge; photons with polarisation, baryons, massless neutrinos; Planck 2018
ΛCDM) with the vapour X as a 24-bin × ℓ ≤ 10 momentum hierarchy coupled to a cold object fluid by elastic scattering,
absorption with sticking s and energy-balanced evaporation, was validated against EH98 (±3%), against the thermal-WDM
transfer (within 4–6% of the CLASS-based fit of arXiv:2210.10753 in k_hm; the Viel fit is ~15–20% warmer), against the
fx_observables two-fluid code in its tight-coupling regime (≤1%), and on CDM/adiabatic controls (≤6×10⁻⁴) [C]. It
changes the vapour physics in three ways [C]:
1. Energy balance on zero-heat-capacity objects forces re-emission after capture freeze-out. This heats the comoving
   momenta (σ_p ≈ 0.038 ζT_f, 3× the "sudden freeze-out" value) and stops net condensation.
2. The vapour temperature lags T_o before Γ/H = 1.
3. Where tight coupling fails, the kinetic suppression is much weaker than the two-fluid estimate.
Over 102 models (f_v = 1–50% by mass, ζ = 10⁻⁷–10⁻², T_f = 0.004–0.014, s = 1 and 0.01) the observables collapse to
**m_eff,Lyα ≈ 19 keV (σ_p/10⁻⁸)^−0.49 (f_v/0.1)^−0.45** (rms 0.065 dex) [C fit]. So ζ_max ≈ 1×10⁻³ (f_v = 3%,
T_f = 0.009) or 2–3×10⁻³ (the physical S/Q = 1.235 line, f_v ≈ 3%, T_f ≈ 0.005) today. A 20 keV-class probe would see
ζ ≳ 0.7–1.6×10⁻⁴ (f_v = 3%). **At qm_micro's computed point (ζ_eff ≈ 0.6–1.0×10⁻⁷, or 5.5×10⁻⁹ with Beliaev
trapping; f_v = 8–26% by mass; T_f ≈ 0.009–0.013) the relic is CDM-like on every scale computed:
|ΔP/P| ≤ 10⁻³ (the numerical floor) at k ≤ 136/Mpc, m_eff > 100 keV for all three probes** [C]. The observable band
needs ⟨Q⟩ ≈ 4×10⁷–3×10¹² at s ≈ 1 (10× less at qm_micro's s ≈ 0.47), i.e. 10⁸–10¹² times smaller objects than
the ripening result.

## 0. Verdict table

| # | Question | Answer | Label |
|---|---|---|---|
| 1 | Kinetic linear Boltzmann solver, Q-matter + vapour on Planck-2018 ΛCDM | built: `code/qt_bg.py`, `qt_hidden.py`, `qt_boltz.py` (§1) | C |
| 2 | Validation | LCDM vs EH98 ±3% (0.05–100/Mpc); thermal WDM within 4–6% (k_hm) of the CLASS-based fit; two-fluid agreement ≤1% where ε_TC ≪ 1; cold-fraction control ≤6×10⁻⁴; adiabatic low-k ≤7×10⁻⁵; Murgia δA reproduced (0.392/0.221 vs 0.38/0.21) (§2) | C |
| 3 | Is the freeze-out picture of fx_formation/fx_observables kinetically right? | Only roughly. (a) Objects cannot store latent heat (thermal energy ~ 10 m_X vs ~10¹⁰ m_X of latent heat per Hubble time at Q = 10¹² [E]), so every absorbed X is re-emitted at T_o and evaporation does not stop at Γ_a = H. After freeze-out this heats σ_p by ×3.1 over "sudden" and ×1.4–1.7 over fx's T_v = T_eq-then-free history, and it blocks post-freeze-out condensation. (b) T_vapour lags T_o (0.57 T_o at Γ/H ≈ 2) (§3) | C/A/E |
| 4 | Tight-coupling failure (fx_observables SKEPTIC §2) | confirmed and quantified on fx's own history. Kinetic vs two-fluid: Q = 10⁸, k = 50: **0.73 vs 0.52**; Q = 10¹⁰, k = 100: 0.80 vs 0.72, k = 300: 0.60 vs 0.23. The fluid codes over-suppress where ε_TC ≳ 1 (§2b) | C |
| 5 | P_mix/P_LCDM(k) grid | 108 runs (102 vapour models + 6 thermal relics), k = 1.5–136/Mpc, z = 0 and 3: `results/transfer_grid.json`. Control parameter σ_p ≈ 0.038 (ζT_f)^0.99; sticking at fixed ζ is secondary (s = 0.01 vs 1: m_eff 19 vs 18 keV at ζ = 10⁻⁴) (§4) | C |
| 6 | Allowed region today (Lyα 5.7 keV area criterion, MW satellites 6.5 keV, lensing 6.1 keV) | ζ_max ≈ 5.8×10⁻³ (f_v = 1%, T_f = 0.004) … 9×10⁻⁴ (3%, 0.009) … 3×10⁻⁴ (10%, 0.009) … 9×10⁻⁵ (35%, 0.009). Physical M1r lines: S/Q = 1.235 (f_v ≈ 3%): 3.0×10⁻³ (Lyα), 3.8×10⁻³ (satellites); S/Q = 3.4 (f_v ≈ 9%): 7×10⁻⁴. Lyα and satellites bind about equally; the lensing proxy is weaker (§5) | C/E |
| 7 | Observable band (allowed now, but beyond a 10–20 keV-equivalent reach of Rubin satellites, lensing samples or next-generation Lyα) | for f_v = 3%: ζ ∈ [0.7–1.6×10⁻⁴ (20 keV) or 3–7×10⁻⁴ (10 keV), 0.9–3×10⁻³]; for 10%: [2–5×10⁻⁵, 3–7×10⁻⁴]; for 35%: [0.7–1.5×10⁻⁵, 0.9–2×10⁻⁴] (§5, `allowed_region.json`, `scaling_fit.json`) | C/E |
| 8 | At qm_micro's values | CDM-like everywhere: ζ_eff = 5.9–9.6×10⁻⁸ (base) / 5.5×10⁻⁹ (Beliaev), σ_p = 1.4–2.2×10⁻¹¹ / 1.7–1.8×10⁻¹², f_v = 8.1–18.7% / 12–26% by mass, max\|ΔP/P\| ≤ 1.2×10⁻³ at k ≤ 136/Mpc. That is 3 decades in ζ below the 20 keV reach. f_v comes out 1.2× qm_micro's at-T_f value and 1.3–1.7× its post-freeze-out value (point 3a) (§6) | C |
| 9 | Mixed-DM literature cross-check | Boyarsky et al. 0812.0010 (any mass for f < 0.35) and Baur et al. 1706.03118 (f ≲ 0.14 (m/keV)^1.37, as quoted) agree with the area criterion for m_eq ≳ 0.5 keV. For the hottest allowed vapours (m_eq ≈ 0.2 keV, ζ ≈ 3×10⁻³, f_v ≈ 3%) the extrapolated Baur relation is stricter, so quote ζ_max ≈ 2–3×10⁻³ there | E/R |

## 1. The calculation [C]

**Background.** Planck 2018 TT,TE,EE+lowE+lensing (arXiv:1807.06209) [R]:
- h = 0.6736, ω_b = 0.02237, ω_c = 0.1200 (all DM: objects + vapour), N_eff = 3.046 massless, T_cmb = 2.7255 K,
  Y_p = 0.2454, n_s = 0.9649;
- recombination from Saha He plus Peebles H (RECFAST fudge 1.14), with Compton-coupled T_m (z* = 1086);
- `qt_bg.py`.

**Perturbations** (`qt_boltz.py`). Synchronous gauge, Ma & Bertschinger 1995 conventions (astro-ph/9506072) [R]:
- photons F_ℓ, G_ℓ ≤ 10, with Thomson coupling;
- baryons with c_s²(T_m);
- neutrino hierarchy ℓ ≤ 12;
- MB95 adiabatic initial conditions at kτ = 0.05.

Approximations (CLASS-like, Blas, Lesgourgues & Tram 2011 [R]):
- zeroth-order photon–baryon tight coupling until kτ_c or τ_c/τ > 0.005;
- RSA for photons and neutrinos after kτ = K_rsa = 200, with δ_r = 4(ℋh′ − k²η)/k² and θ_r = −(h′ + 6η′)/2 solved
  implicitly with the constraints;
- before recombination the baryons are tied to the slow photon velocity, then released with drag (R/τ_c ℋ < 300);
- after a_end = max(kτ = 400, end of vapour collisions) a sub-horizon Newtonian phase: objects and baryons by
  predictor–corrector, the vapour by the exact collisionless Volterra (Gilbert) integral. Its free term continues the
  q-resolved state on the same q grid (Σ_ℓ (2ℓ+1)(−1)^ℓ ∫ δN_ℓ j_ℓ(kqS)); its kernel is W(kS) from the background N(q).

BDF with the exact Jacobian (the system is linear), rtol = 10⁻⁷.

**DM sector.** Variables:
- objects: charge contrast δ_oc, number contrast δ_N, velocity θ_o;
- vapour: per-charge momentum distribution N(q) on a fixed comoving log-q grid (24 bins, Δln q ≈ 0.38, covering the
  distribution from a_kbg (Γ/H = 60) to the hot re-emission tail at 25 q_f), as a hierarchy δN_ℓ(q_j), ℓ ≤ 10.
  The state holds bin populations w_j δN_j.

Collision operator in the object rest frame:
- elastic Γ_el(q)[⟨f⟩ − f] (hard sphere, isotropic);
- absorption −Γ_a(q) f, emission +Γ_a(q) f_eq(q; T_o), with Γ_a = s πR² n_obj v = ζ K n_Q C_o^{2/3} v and
  K = π(3/4πn_B)^{2/3};
- σ_el = (1 − s)πR², so that πR² is the total momentum-transfer cross-section;
- the object temperature T_o comes from the objects' instantaneous energy balance, Σ w Γ_a (Δ + E_k)(N − N_eq(T_o)) = 0
  (zero heat capacity: §3).

The operator is linearised exactly. This includes the boost terms (dipole sources −(θ_o/3kv)[Γ_el D + Γ_a D_eq +
q∂_qΓ_a (N_eq − N)]), the rate perturbation δΓ/Γ = δ_N/3 + 2δ_oc/3, and the momentum, charge and energy exchange with
the objects. The objects take (1 − Δ) m_X per absorbed charge. The derivative D = q∂_qN − 2N is moment-corrected:
Σw q^n D = −(n+3) Σw q^n N for n = 0, 1, 2. Without that correction, a cold control is off by 1.45% × f_v (§2c).

Before a_kin(k) = earliest of (Γ_m/H = 30, ε_TC = kv̄/aΓ_m = 0.03) the DM is one fluid with the two-phase adiabatic
c_eff² (fx_observables). It is mapped to the kinetic state by local equilibrium (δN = N[rδ + (E_k/T − 3/2)δ/D_f]).
The adiabatic IC in the kinetic phase is δN = −[∂N/∂ln a + D]δ_Q/3.

**Hidden background** (`qt_hidden.py`). The quasi-static two-phase equilibrium:
- n_sat = g(T/2π)^{3/2}e^{−Δ/T}, Δ = 1 − Ω(Q), g = 1 for fx_formation's ideal vapour;
- the adiabat N_v(5/2 + Δ/T) = S/Q. It is conserved for any g [D: d[N(5/2 + Δ/T)] = 0 follows from energy
  conservation with n_v = n_sat];
- n_Q fixed by ω_c and the mass per charge (1.2× the old M1 Y_Q).

From Γ/H = 60 onward the full kinetic background N(q, a) is integrated with T_o(a). The modes are:
- 'thermo': physical, the energy balance throughout;
- 'sudden': exchange switched off at Γ̄_a = H, fx's picture;
- 'fx': fx_observables' exact history, for validation;
- 'frozen': collisionless Fermi-Dirac or Maxwellian, for the WDM validation.

Optional energy-dependent sticking s(E_k) uses qm_micro's stick_thermal.json table.

## 2. Validation [C] (`results/validation.json`, `validation_lcdm.json`)

**a. ΛCDM** (`val_lcdm.py`). T(k)/T_EH98,no-wiggle (astro-ph/9709112) [R]:

| k (1/Mpc) | 0.05 | 0.2 | 1 | 3 | 10 | 30 | 100 |
|---|---|---|---|---|---|---|---|
| K_rsa = 800 | 1.020 | 0.986 | 1.007 | 1.004 | 1.032 | 1.030 | 1.017 |
| K_rsa = 200 (production) | same | same | 1.008 | 1.028 | 1.044 | 1.121 | 1.173 |

- At k = 0.05–0.2 the residuals are the BAO wiggles that the no-wiggle fit omits.
- RSA during radiation domination is not converged at kτ = 50–200, because there the oscillating part of η is as large
  as its slow part. That is why CLASS uses RSA only after recombination.
- The model/ΛCDM ratio is insensitive to it. The run-to-run change in P_mix/P_LCDM between K_rsa = 200 and 400 is:
  - ≤ 4×10⁻⁵ at k ≤ 8;
  - 0.5% (ζ = 10⁻³) and 1.2% (3×10⁻³) at k = 24;
  - 3% and 9% at k = 68, with production slightly under-suppressed.
- Tolerance: rtol 10⁻⁷ vs 10⁻⁸ agree to ≤ 10⁻⁴ (controls c4).

**b. Two-fluid code on fx_observables' exact history** (`val_twofluid.py`). It imports skeptic/code/two_fluid_check.py
read-only and compares with their two_fluid_transfer.json. Kinetic | two-fluid (P at z = 0):

| case | k = 3 | 10 | 20 | 50 | 100 | 300 |
|---|---|---|---|---|---|---|
| Q = 10¹⁴, s = 1 (ζ = 2×10⁻⁵) | 1.000 \| 1.000 | 0.9997 \| 0.9998 | 0.9987 \| 0.9989 | 0.9916 \| 0.9910 | 0.970 \| 0.965 | 0.919 \| 0.904 |
| Q = 10¹² (10⁻⁴) | 0.9997 \| 0.9998 | 0.9958 \| 0.9965 | 0.9835 \| 0.9840 | 0.946 \| 0.934 | 0.913 \| 0.903 | 0.821 \| 0.785 |
| Q = 10¹⁰ (4.6×10⁻⁴) | 0.996 \| 0.997 | 0.963 \| 0.962 | 0.932 \| 0.922 | 0.874 \| 0.861 | **0.800 \| 0.720** | **0.604 \| 0.230** |
| Q = 10⁸ (2×10⁻³) | 0.971 \| 0.970 | 0.907 \| 0.901 | 0.851 \| 0.819 | **0.731 \| 0.521** | **0.592 \| 0.210** | 0.272 \| 0.140 |
| Q = 2×10⁵, s = 0.01 (332/m_X² elastic) | 0.999 \| 1.000 | 0.985 \| 0.995 | 0.957 \| 0.979 | 0.846 \| 0.870 | 0.641 \| 0.593 | 0.040 \| 0.053 |

- Where tight coupling holds (the skeptic's ε_TC < 1), agreement is ≤ 1%.
- Where it fails, the kinetic vapour leaks out of the coupled fluid and the acoustic and damping suppression of the
  fluid codes is not there. This is the skeptic's open item: a "BGK or Boltzmann treatment ... is the calculation that
  would decide the k ≥ 10/Mpc values for ζ ≳ 5×10⁻⁴".
- The Q = 2×10⁵, s = 0.01 case (long elastic-only coupling, diffusive regime) gives about 2× more suppression at
  k = 10–20 than the fluid. This is the diffusion damping that the fluid with momentum exchange alone misses.

**c. Thermal WDM** (Fermi-Dirac, g = 2, relativistic early, background e(a) included; scan family W). Half-mode
k_hm (P = 1/4), in 1/Mpc:

| m (keV) | own code | CLASS-based fit (Vogel & Abazajian 2022, arXiv:2210.10753) | Viel et al. 2005 fit (astro-ph/0501562) |
|---|---|---|---|
| 3.5 | 29.8 | 28.6 | 26.3 |
| 5.3 | 48.7 | 46.3 | 41.7 |
| 5.7 | 53.9 | 50.9 | 44.8 |
| 6.5 | 62.7 | 59.5 | 52.3 |

- Our relic is colder than the Viel fit, in the same direction and within 4–6% of the CLASS-based fit. That paper finds
  CLASS "colder" than earlier fits, with mass corrections of ≈ 14% at 10 keV [R].
- Internal convergence (3 keV): Nq 16–44, ℓ 8–24, K_end 100–3000 agree to ≤ 0.5%.
- A Maxwellian with the same ⟨q²⟩ suppresses more than the FD relic (P(20) = 0.352 vs 0.387 at 3 keV equivalent).
- The constraints use the Viel curves at the published bound masses as the allowed suppression, since that form is
  behind the bounds. The own-code curves are reported as the systematic (`judge_own`).

**d. Controls** (`controls.py` → validation.json): all pass, except that 16/8 is too coarse.

| Control | Result |
|---|---|
| Cold collisionless fraction (σ = 10⁻¹¹; f = 0.03, 0.3) | P = 1 + (0.4–2)×10⁻⁵ at k = 1, 10; +6.3×10⁻⁴ at k = 68 |
| Full kinetic model at low k (ζ = 10⁻⁴) | 1 − 7×10⁻⁵ (thermo), 10⁻⁶ (sudden) at k = 0.3–1 |
| Resolution: Nq/ℓ = 24/10 vs 36/14 | ≤ 4×10⁻⁴ |
| Resolution: Nq/ℓ = 16/8 | fails (0.947 vs 0.762 at k = 48) |

**e. Area criterion** (Murgia et al. 2017, arXiv:1704.07838 [Q, read]). Integrating P_1D = (1/2π)∫_k^{k_up} k′P dk′
up to k_up = 75 h/Mpc reproduces their δA_REF = 0.38 (3.5 keV) and 0.21 (5.3 keV) for the Viel spectra: we get 0.392
and 0.221. With k_up = ∞ the values are 0.51 and 0.38. All model/reference comparisons use the same implementation.

## 3. What the kinetic treatment changes in the vapour physics [C/A]

1. **Re-emission after freeze-out is forced by energy conservation** [A + C].
   - An absorbed X releases Δ + E_k ≈ 0.22 m_X into the object. Q-matter has no light neutral channel and a negligible
     heat capacity (phonon gas ≈ 14 m_X of thermal energy at T = 0.007 for Q = 10¹², against ≈ 10¹⁰ m_X of latent
     heat absorbed per Hubble time near freeze-out).
   - So the heat must leave as evaporated X: about one X re-emitted at T_o per X absorbed. The net condensation is only
     2T_o/(Δ + 2T_o) ≈ 6% of the absorption.
   - After Γ_a = H the vapour therefore keeps being re-thermalised at a rate (a_f/a)² per e-fold (RD). That adds a
     constant amount a_f²·3T_o per e-fold to ⟨q²⟩, a log-growing hot tail.
   - Result: **σ_p,thermo = 0.038 (ζT_f)^0.99 vs σ_p,sudden = 0.0124 (ζT_f)^0.99** [C fit, rms 0.07 / 0.009 dex].
     With fx_observables' own inputs (Q = 10¹², S/Q = 1.235, their n_Q) its history_B gives σ_p = 1.51×10⁻⁸, between
     our thermo 2.13×10⁻⁸ (×1.41) and sudden 0.59×10⁻⁸ (×0.39).
   - The vapour fraction does **not** fall after freeze-out (f_v,charge = 0.0287 at a_f, 0.0282 today, for Q = 10¹²).
   - The brief's item (iii), "evaporation until capture freeze-out", is therefore replaced by evaporation balanced
     against absorption at all times.
2. **Temperature lag before freeze-out** [C]. Exchange at Γ/H of a few cannot keep the vapour at T_o against 2H
   adiabatic cooling, and fast particles are absorbed preferentially (Γ ∝ v). The ratio T_vapour/T_o is 1.0 at
   Γ/H = 60, 0.93 at 21, 0.80 at 7, 0.57 at 2 and 0.30 at 0.5 (Q = 10¹², S/Q = 1.235). The kinetic freeze-out
   (Γ̄_a = H with the actual distribution) comes 1.7× earlier in a than the equilibrium estimate (1.27 vs 2.18×10⁻⁷).
3. **Tight-coupling failure** (§2b): the kinetic vapour gives far less suppression at k ≳ 1/r_s than the fluid models
   for ζ ≳ 5×10⁻⁴.
4. **Sticking at fixed ζ is secondary** [C]. With s = 0.01 (elastic (1 − s)πR², momentum decoupling 10× later in a) and
   s = 1 at the same ζ, m_eff,Lyα = 19.1 vs 18.2 keV (ζ = 10⁻⁴) and 22 vs 12 keV (3×10⁻⁴). The long elastic phase
   lowers σ_p.

## 4. The P_mix/P_LCDM grid [C] (`results/transfer_grid.json`; all at z = 0, z = 3 also stored)

Families (scan.py):
- **P**: physical g = 1, thermo, s = 1; S/Q = 0.5, 1.235, 3.4.
- **S**: sudden.
- **E**: s = 0.01.
- **X**: parametric, with target vapour charge fraction f_c ∈ {0.03, 0.1, 0.35} and T_f ∈ {0.004, 0.010}, realised by a
  saturation factor g and adiabat S/Q. g spans 10⁻⁸–10⁵, i.e. a latent-heat shift −T ln g. f_c = 0.35 at T_f = 0.004
  needs S/Q ≈ 20 and is two-phase only below T ≈ 0.012; it is outside M1r's survival window.
- **M**: qm_micro points.
- **W**: thermal relics 3.5–20 keV.

Physical line S/Q = 1.235, s = 1 (f_v,mass = 2.6–6.4%), P at k = 5 / 12 / 24 / 48 / 96 per Mpc:

| ζ | σ_p (c) | P(5) | P(12) | P(24) | P(48) | P(96) | δA | N(>10⁸)/CDM | m_eff Lyα/sat/lens (keV) |
|---|---|---|---|---|---|---|---|---|---|
| 10⁻⁷ | 4.2×10⁻¹¹ | 1.000 | 1.000 | 1.000 | 0.999 | 1.001 | 0.001 | 1.000 | >100 |
| 10⁻⁵ | 2.9×10⁻⁹ | 1.000 | 1.000 | 0.999 | 0.995 | 0.987 | 0.002 | 0.999 | 58/56/59 |
| 3×10⁻⁵ | 8.2×10⁻⁹ | 1.000 | 0.998 | 0.994 | 0.983 | 0.963 | 0.005 | 0.996 | 31/30/37 |
| 10⁻⁴ | 2.5×10⁻⁸ | 0.998 | 0.990 | 0.976 | 0.951 | 0.910 | 0.017 | 0.986 | 18/18/24 |
| 3×10⁻⁴ | 7.1×10⁻⁸ | 0.990 | 0.970 | 0.940 | 0.892 | 0.812 | 0.042 | 0.967 | 12/13/18 |
| 10⁻³ | 2.1×10⁻⁷ | 0.968 | 0.923 | 0.862 | 0.762 | 0.602 | 0.099 | 0.927 | 8/9/13 |
| 3×10⁻³ | 5.6×10⁻⁷ | 0.926 | 0.843 | 0.730 | 0.559 | 0.335 | 0.194 | 0.854 | 5.7/6.8/9.5 |
| 10⁻² | 1.3×10⁻⁶ | 0.867 | 0.728 | 0.557 | 0.336 | 0.117 | 0.321 | 0.744 | 4.1/5.3/7.3 |

- At S/Q = 3.4 (f_v ≈ 7–18%) the same ζ gives about twice the deficit: P(24) = 0.936 at 10⁻⁴ and 0.659 at 10⁻³.
- ζ ≥ 3×10⁻³ shows the coupled-fluid regime: the deficit extends to k = 1.5/Mpc (4% at ζ = 10⁻²).
- **Scaling** (`results/scaling_fit.json`, 63–72 models with 1.5 < m_eff < 80 keV):
  - m_eff,Lyα = 19.3 keV (σ_p/10⁻⁸)^−0.49 (f_v/0.1)^−0.45 (rms 0.065 dex, max 0.16);
  - satellites: 19.5 (…)^−0.44 (…)^−0.43;
  - lensing: 23.2 (…)^−0.36 (…)^−0.38;
  - with σ_p ≈ 0.038 ζT_f these give the limits in §5.

## 5. Constraints, allowed region and observable band [C on E proxies]

Criteria (`constraints.py`), each applied identically to the model and to the thermal reference:

| Probe | Criterion | Bound |
|---|---|---|
| Lyα | area criterion δA, z = 0 linear, 0.5–20 h/Mpc | ≤ δA(Viel 5.7 keV) = 0.195, Iršič et al. 2024, arXiv:2309.04533 [R]; also 5.3 keV (1702.01764: 0.221) and the conservative 3.5 keV (0.392) |
| MW satellites [E] | N(>10⁸ M_⊙) ratio, sharp-k + Sheth–Tormen | ≥ that of 6.5 keV (Nadler et al. 2021, arXiv:2008.00022), 0.837 |
| Strong lensing [E] | dn/dlnM at 10^7.5 M_⊙ | ≥ 6.1 keV (Keeley et al. 2024, arXiv:2405.01620), 0.397; 5.2 keV (Gilman et al. 2020, arXiv:1908.06983) as the 2020 bound |
| Future [R forecasts used only as thresholds] | the same metrics at 10 and 20 keV | Rubin-era satellites forecast 10–20 keV (Nadler et al. 2024); larger lensing samples; a DESI-class high-resolution Lyα follow-up |

Villasenor et al. 2023 (arXiv:2209.14220) give m > 3.1 keV (95%), so the 5.7 keV bound used here is the tighter
choice. The half-mode criterion is never triggered for f_v < 0.5 (P ≥ (1 − f_v)² asymptotically); the area and
abundance criteria are used instead.

Results (`results/allowed_region.json` from the grid, `scaling_fit.json` from the fit; ζ at s = 1):

| f_v (mass) | T_f | ζ_max today (most constraining probe) | 10 keV reach | 20 keV reach |
|---|---|---|---|---|
| 1% | 0.004 / 0.009 | 5.8×10⁻³ / 2.6×10⁻³ (Lyα) | 1.8×10⁻³ / 8×10⁻⁴ | 4.3×10⁻⁴ / 1.9×10⁻⁴ |
| 3% | 0.004 / 0.009 | 2.1×10⁻³ / 9×10⁻⁴ (Lyα ≈ sat) | 6.6×10⁻⁴ / 2.9×10⁻⁴ | 1.6×10⁻⁴ / 7×10⁻⁵ |
| 10% | 0.004 / 0.009 | 6.9×10⁻⁴ / 3.1×10⁻⁴ | 2.2×10⁻⁴ / 9.6×10⁻⁵ | 5.2×10⁻⁵ / 2.3×10⁻⁵ |
| 20% | 0.004 / 0.009 | 3.5×10⁻⁴ / 1.6×10⁻⁴ (sat) | 1.1×10⁻⁴ / 5×10⁻⁵ | 2.6×10⁻⁵ / 1.2×10⁻⁵ |
| 35% | 0.004 / 0.009 | 2.0×10⁻⁴ / 9×10⁻⁵ (sat) | 6.8×10⁻⁵ / 3.0×10⁻⁵ | 1.5×10⁻⁵ / 6.8×10⁻⁶ |

Physical M1r lines (direct grid interpolation):
- S/Q = 1.235: ζ_max = 3.0×10⁻³ (Lyα 5.7 keV), 3.8×10⁻³ (satellites), > 10⁻² (lensing proxy). The 20 keV reach is
  8×10⁻⁵ and the 10 keV reach 5–7×10⁻⁴.
- S/Q = 3.4: ζ_max = 7.0×10⁻⁴ / 7.5×10⁻⁴ / 2.6×10⁻³; reach 3×10⁻⁵ (20 keV).
- S/Q = 0.5 (f_v ≈ 1–2.5%): nothing in the grid is excluded today.
- With the colder own-code relics as references (same mass label, less suppression allowed), the limits tighten by
  ≈ 30% in ζ (S/Q = 1.235: 2.1×10⁻³ instead of 3.0×10⁻³ for Lyα; `judge_own`).

The observable band is the interval between the reach and ζ_max in the table: about 1.1–1.6 decades in ζ at
20 keV and 0.5 decade at 10 keV. In objects at s ≈ 1 that is ⟨Q⟩ ≈ 4×10⁷–3×10¹². fx_observables placed the edge between 5×10⁻⁴ (allowed) and 2×10⁻³ (marginal in the two-fluid
re-run). The kinetic edge is 2–3×10⁻³ at f_v ≈ 3% (point 3 of §3, partly offset by the thermo heating) and
7×10⁻⁴ at f_v ≈ 9%.

## 6. Evaluated at qm_micro's values [C] (`results/qm_micro_eval.json`)

Inputs from tracks/qm_micro/NOTES.md §3 [Q]: S/Q = 1.5 / 2.5 / 3.2, ⟨Q⟩ = 5.5×10¹⁹ / 1.2×10²⁰ / 2.7×10²⁰, sticking
s(E_k) from stick_thermal.json (angle-averaged, threshold law 5.31√E_k). Beliaev variant: trapping factor 0.08–0.15
on s, ⟨Q⟩ = 2.3×10²⁰ / 1.6×10²¹. Thermo kinetics:

| point | T_f (equilibrium; theirs) | ζ_eff(T_f) (theirs) | f_v charge / mass (theirs, at T_f; after) | σ_p today | max\|1 − P\|, k ≤ 136 | m_eff |
|---|---|---|---|---|---|---|
| base S/Q 1.5 | 0.0084 (0.0083) | 9.6×10⁻⁸ (1.2×10⁻⁷) | 0.064 / 0.081 (0.054; 0.038–0.050) | 2.2×10⁻¹¹ | 1.0×10⁻³ | > 100 keV all probes |
| base 2.5 | 0.0089 (0.0087) | 7.5×10⁻⁸ (9.5×10⁻⁸) | 0.113 / 0.141 (0.094; 0.067–0.086) | 1.7×10⁻¹¹ | 7×10⁻⁴ | > 100 |
| base 3.2 | 0.0092 (0.0091) | 5.9×10⁻⁸ (7.3×10⁻⁸) | 0.151 / 0.187 (0.126; 0.089–0.115) | 1.4×10⁻¹¹ (sudden 6.0×10⁻¹²) | 7×10⁻⁴ | > 100 |
| Beliaev 1.5 / 3.2 | 0.0121 / 0.0127 (0.0119 / 0.0126) | 5.5×10⁻⁹ (6×10⁻⁹) | 0.095 / 0.217 charge | 1.8 / 1.7×10⁻¹² | ≤ 1.2×10⁻³ | > 100 |

- **The qm_micro relic is CDM-like on every scale computed here** (to the 10⁻³ numerical floor at k ≤ 136/Mpc;
  m_eff > 100 keV). It is a factor ≈ 10³ (base) or 10⁴ (Beliaev) in ζ below the 20 keV reach.
- Its warm vapour carries 8–26% of the DM mass but has σ_p ≲ 2×10⁻¹¹ (λ_fs ~ 10 pc).
- The vapour fraction is higher than qm_micro's post-freeze-out values because the object energy balance forbids
  capture-only condensation after freeze-out (§3.1). qm_micro's variant (b), "decoupled vapour cooling, capture only",
  is the case affected.
- This changes M/m_X by +0.224 Δf_v ≈ +0.004–0.011 m_X (0.5–1.4%, more vapour), not the observables.

## Limits

- **Linear physics only.** The criteria are proxies:
  - the area criterion on linear z = 0 P(k);
  - sharp-k Sheth–Tormen abundances at z = 0 for satellites and lensing [E]. These are not the published likelihoods,
    and subhalo stripping and galaxy formation are not modelled.
  - The future thresholds are forecasts quoted from abstracts [R].
- **Vapour model.** Ideal NR vapour, no X–X scattering (n_v σ_XX/n_obj πR² ≈ 10⁻⁴); single object size Q (no
  spectrum); hard-sphere isotropic elastic scattering; zero object heat capacity; constant s except in the M family;
  Gibbs–Thomson shift of Δ frozen. Sudden/thermo bracket the post-freeze-out physics.
- **Before T_h < 0.05** the hidden sector is assumed to be in equilibrium (fx_formation [Q]).
- **Numerics.** The RSA switch at kτ = 200 is not converged in absolute amplitude at k ≳ 30 (ratios to ≤ 1.2% at
  k ≤ 24, ≤ 9% at k = 68 for the most collisional cases). The 10⁻³ ratio floor at k ≥ 68 limits statements about
  ζ ≲ 10⁻⁵.
- **Background.** Massless neutrinos; constant-N_eff radiation before e⁺e⁻ annihilation (H high by ≤ 10% at
  a ≲ 10⁻⁹, affecting ζ ≲ 10⁻⁶ freeze-out timing only).

## Reproduction (OMP_NUM_THREADS=1; scratch $SP/wt2/qm_transfer; from code/)

```
python3 -I scan.py $SP/wt2/qm_transfer/runs/scan 0 2 585 ; python3 -I scan.py ... 1 2 585   # repeat until done (~14 chunks of 10 min, 2 workers)
python3 -I val_twofluid.py $SP/.../val_twofluid.json 0,1,2,3,4        # ~10 min
python3 -I controls.py $SP/.../controls.json                            # ~12 min
python3 -I val_lcdm.py ../results/validation_lcdm.json                  # ~3 min
python3 -I summary.py $SP/.../scan $SP/.../val_twofluid.json $SP/.../controls.json   # writes results/*.json
python3 -I fit_scaling.py                                               # results/scaling_fit.json
```
`val_wdm.py` (early thermal-relic checks) is superseded by scan family W.
