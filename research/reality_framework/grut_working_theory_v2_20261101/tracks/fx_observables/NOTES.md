# Track fx_observables: what M1 shows an observer, and the smallest extension that would show more

1 November 2026 · Claude, issue #13 (BRIEF.md of grut_working_theory_v2_20261101) · same-session AI work on one supplied
model; not outside replication. Nothing here is a validation, an exclusion of the real world, or a new theory. All data
comparisons are retrospective. **Every extension below is a model choice, not a consequence of closure.** Git was used
read-only, nothing was committed, and no Codex code was run. Labels: [D] derived, [C] computed with our code,
[E] estimate, [A] argued, [Q] quoted, [R] recalled from the literature.

## Result in one paragraph

On every halo-scale and particle observable, conservative M1 (zero SM portals) is dark matter that interacts only
through gravity. The objects are 8–13 fm, 17 TeV composites with σ/M = 3–7×10⁻⁴ cm²/g (Q = 2×10⁵), which is 2–4 decades
below any SIDM effect. Leftover X̄ are burned out on the objects (≈10⁵ e-folds after 2 MeV), and mergers release ≲10⁻³ of
the local DM mass, only into hidden quanta. The sector has no light hidden state, so ΔN_eff ≲ 10⁻⁸. **The one
cosmological channel that is not CDM is new to this pass.** In the branch where objects survive (tracks/fx_formation),
latent heat holds a 2–6% X vapour at v ≈ 0.1c in equilibrium with the objects until vapour capture freezes out. Until
then that vapour couples the whole DM into one fluid with c_eff ≈ 0.01c, and afterwards it free-streams as a warm
component. One parameter controls the size of the effect, ζ = s·Q^(−1/3), where s is the sticking probability. The
result is CDM at k ≤ 50/Mpc to better than 1% for ζ ≲ 2×10⁻⁵ (e.g. ⟨Q⟩ ≳ 10¹⁴ at s = 1). It gives a 7–11% deficit at
50/Mpc for ζ ≈ 10⁻⁴ (⟨Q⟩ ≈ 10¹², s = 1). It is in conflict with the Lyman-α/satellite bounds for ζ ≳ 10⁻³. That last
case includes the reference Q = 2×10⁵ with s ≳ 0.1. Codex's pulse gas, which is purely elastic, stays CDM to
|ΔP/P| < 2×10⁻⁵ at k ≤ 100/Mpc. A direct-detection signal needs a portal. The φ²|H|² portal fails: at the coupling
XLZD could see (b ≈ 2×10⁻¹⁰), its own freeze-in adds ΔS/Q ~ 10³–10⁴ of hidden entropy per charge, against a survival
budget of ≈ 2. Objects survive only for b ≲ 1.5×10⁻¹² (≈ 5×10⁻⁴ XLZD events). The pure-mixing trilinear A φ|H|², which
the M1 tilt allows, is the one extension that could be detectable (A ≈ 3×10⁻¹⁰ GeV at XLZD) and still roughly fit the
entropy budget (ΔS/Q ≈ 0.4–4 [E]). Its signal is a soft, light-WIMP-like Xe spectrum with a form-factor zero at
34.1 keV, and a deficit in argon. Even then it requires ⟨Q⟩ ≲ 10⁶, a narrow Q spread, and s ≲ 10⁻². fx_formation's
ripening argues against all three.

## 0. Verdict table

| # | Question | Answer | Label |
|---|---|---|---|
| 1 | Object properties, Q = 2×10⁵ | M = 17 324.83 GeV (E/Q = 0.937400 m_X); charge radii r₅₀/r₉₀/r₉₉ = 3.79/5.94/8.42 m_X⁻¹, φ-half radius 4.25, bulk-density sphere 4.74 m_X⁻¹ (1/m_X = 2.135 fm: R ≈ 8–13 fm); object–object σ = π(2R)² = 1.8–4.4×10² m_X⁻² = 0.8–2.0×10⁻²³ cm²; **σ/M = 2.7–6.6×10⁻⁴ cm²/g** (1.1×10⁻³ at Q = 10⁵, 5.7×10⁻⁵ at 10⁸, 2.9×10⁻⁶ at 10¹², 6×10⁻⁹ at 10²⁰) against SIDM 0.1–10 cm²/g. Halo collisions are classical (pR ≈ 10³) and inelastic: they merge and release 3.7% of the pair mass (1.29 TeV), or fission (Codex's 3D pairs re-split 1:2) | C |
| 2 | Gas coupling and CDM status | (A) Codex's gas with the elastic 332/m_X² coupling: Γ/H = 1.2×10⁶ at 2 MeV, decouples at T_vis = 2.6 keV, then free-streams 25 pc rms; sound horizon of the coupled fluid 3.7 pc; P/P_CDM = 1 − O(10⁻⁵) at k ≤ 100, 0.9986 at 10⁴, 0.986 at 3×10⁴/Mpc. (B) Equilibrium vapour (NEW): a table in ζ = sQ^(−1/3), §2. CDM at k ≤ 50 (≥ 99%) only for ζ ≲ 2×10⁻⁵ | C / E |
| 3a | Leftover X̄ | burned out on the objects: X̄ + core → X + core(Q−2), σv ≈ 1.95/m_X², Γ/H = 1.2×10⁵ at 2 MeV (Q = 2×10⁵), so p_today ≈ e^(−10⁵). If p = 0.005 survived anyway: 2.7% would annihilate per t_H at the Sun, each event putting 0.163 GeV into one X at 0.93c (invisible, escapes the halo) | C/E |
| 3b | Mergers today | 6×10⁻³ (Sun), 0.09 (inner kpc), 0.07 (cluster core) collisions per object per t_H at Q = 2×10⁵; mass released ≤ 1.6×10⁻³, into hidden quanta only, so no dark radiation | C |
| 3c | Hidden radiation | no hidden state lighter than X (92 MeV); vector heavy (3×10⁶–3×10⁷ GeV) and decays to XX̄ at t ≈ 3×10⁻⁷–3×10⁻⁴ s. At BBN the vapour is non-relativistic with T_h/T_vis = 0.85 (1 MeV) and 6.4 (0.1 MeV), ΔN_eff = 2×10⁻⁹–2×10⁻⁸; even counting all DM energy as radiation, < 4×10⁻⁶ at 1 MeV | C |
| 4a | Higgs portal b φ²\|H\|² | zero of the Xe recoil spectrum at 27.56 keV (Ge 48.80, Ar 90.27; reproduces Codex); σ_n0 = 1.40×10⁻²⁶ b² cm²; LZ-2024-equivalent b < 6.5×10⁻¹⁰; XLZD 200 t·yr: 10 events at b = 2.2×10⁻¹⁰. At that b: BR(h→inv) = 2×10⁻¹⁶, f_rest = 181 (erasure 99.994%), **ΔS/Q ≈ 3×10³–4×10⁴ ≫ 2.2**, so objects evaporate. Survival needs b ≲ 1.4–1.7×10⁻¹² | C/E |
| 4a′ | Trilinear A φ\|H\|² (pure mixing) | zero at 34.11 keV Xe (Ge 60.4, Ar 111.7); LZ A < 8.2×10⁻¹⁰ GeV; XLZD 10 events at A = 2.7×10⁻¹⁰ GeV (θ = 4.3×10⁻¹²); freeze-in ΔS/Q ≈ 0.4–4 (top-Compton dominated), marginal against the budget | C/E |
| 4a″ | λ_HX \|X\|²\|H\|² | zero at 41.68 keV Xe; LZ λ < 2.7×10⁻¹⁰; XLZD 10 events at 9×10⁻¹¹; h → XX̄ freeze-in ΔS/Q ≈ 3×10³: fails like b | C/E |
| 4b | Kinetic mixing | zero at all orders in conservative M1. With Option 1A (hidden C broken by α) and a Higgs portal it needs ≥ 2 portal insertions plus SM CP violation: negligible | A |
| 4c | Gravity, stars | M = 3×10⁻²⁰ g (Q = 2×10⁵) to 1.3×10⁻⁵ g (10²⁰): 27+ decades below microlensing. Only Q ≳ 10²⁰ (≈ M_Pl) reaches the proposed gravitational direct detection (flux ≈ 0.3 m⁻² yr⁻¹). No portal means no stellar or NS capture; with a portal, coherence at NS momentum transfer is \|F/F₀\|² = 10⁻⁸ (0.3 GeV) to 10⁻¹³ (0.5 GeV) | C/R |
| 5 | Verdict | conservative M1 is secluded and CDM-like on every accessible scale **if** ζ ≲ 2×10⁻⁵ (big objects or poor sticking); otherwise the late vapour is its one observable (§2), and it is generic mixed/acoustic DM. Direct detection needs a portal; only the pure-mixing trilinear survives the entropy budget, marginally | E |

## 1. Object properties [C] (`code/objects.py` → `results/objects.json`)

The profile solver is `code/qball_m1_copy.py`, an unchanged copy of tracks/m1_vacuum `qball_m1.py` (Numerov plus
Newton). The driver raises the Newton tolerance above the round-off floor for Q ≥ 3×10⁷, and continues the solution in
small steps. Every E/Q of m1_vacuum and fx_formation is reproduced to its quoted digits (2×10⁵: Ω = 0.881720,
E/Qm = 0.937400; 10⁸: 0.787772 and 0.793978). The charge check Ω∫F²/g²Q equals 1 to 10⁻⁶.

| Q | M (GeV) | r₅₀ / r₉₀ / r_S / R_bulk (m_X⁻¹) | σ/M (r₉₀), cm²/g | merger release E(Q)+E(Q)−E(2Q) |
|---|---|---|---|---|
| 1×10⁵ | 9 046.5 | 3.30 / 5.62 / 3.17 / 3.76 | 1.12×10⁻³ | 4.2% (768 GeV) |
| **2×10⁵** | **17 324.8** | **3.79 / 5.94 / 4.25 / 4.74** | **6.6×10⁻⁴** (2.7×10⁻⁴ with r₅₀) | **3.7% (1.29 TeV)** |
| 1×10⁶ | 80 121 | 6.14 / 8.48 / 7.77 / 8.10 | 2.9×10⁻⁴ | 2.3% |
| 1×10⁷ | 7.54×10⁵ | 13.4 / 16.9 / 17.2 / 17.5 | 1.2×10⁻⁴ | 1.1% |
| 1×10⁸ | 7.34×10⁶ | 29.3 / 35.9 / 37.4 / 37.6 | 5.7×10⁻⁵ | 0.48% |
| 10¹² / 10¹⁴ / 10²⁰ (bulk spheres, fx_formation fit) | 7.2×10¹⁰ / 7.2×10¹² / 7.2×10¹⁸ | R = 810 / 3 760 / 3.8×10⁵ | 2.9×10⁻⁶ / 6.3×10⁻⁷ / 6.3×10⁻⁹ | 2×10⁻⁴ / 5×10⁻⁵ / 5×10⁻⁷ |

- **SIDM.** σ/M ∝ Q^(−1/3) at large Q. Across the whole bracket it is 10²–10⁴ below the range that changes halos:
  0.1–10 cm²/g in dwarfs and low-surface-brightness galaxies (Tulin & Yu 2018, arXiv:1705.02358; Kaplinghat, Tulin & Yu
  2016, arXiv:1508.03339) [R], and below the cluster bounds of ≲ 0.1–1 cm²/g (Randall et al. 2008, arXiv:0704.0261,
  < 1.25; Harvey et al. 2015, arXiv:1503.07675, < 0.47; Sagunski et al. 2021, arXiv:2006.12515, < 0.35) [R]. Neither the
  Yukawa ranges (φ: 1/m_φ = 0.41 m_X⁻¹ < R; vector: ~10⁻⁸ fm) nor gravity enlarge the cross section. In halos pR ≈ 1100
  (Q = 2×10⁵) to 3×10⁶ (10⁸), so the scattering is classical and geometric.
- **Inelasticity.** At v ~ 10⁻³ the kinetic energy per charge is ~10⁻⁷ m_X, while the merger release is 1.29 TeV
  (13 975 m_X) for 2×10⁵ + 2×10⁵. The collision is therefore always inelastic. The energy can only go into hidden quanta:
  free X (each costs 1 − Ω(4×10⁵) = 0.143 m_X), φ → XX̄, or fission. The pair runs of Codex (ed73791: merge, then re-split
  1:2) and the literature (Battye & Sutcliffe 2000, hep-th/0003252; Axenides et al. 2000, hep-th/9910388 [R]: charge
  exchange, merger or bounce depending on relative phase) agree that halo-velocity collisions are not elastic. At these
  rates (§3b) it makes no difference.

## 2. The free-X gas and the matter power spectrum (`code/vapour_transfer.py`, `controls_vt.py`, `summary.py`, `fs_lengths.py`)

**Linear model [E].** The model is the sub-horizon Newtonian limit on a Planck-like background with g*(T) steps.
Baryons are smooth until z = 1100 and pressureless after. The initial condition is Hu–Sugiyama log growth at kτ = 10 [R].
Only the ratio P_M1/P_CDM at z = 0 is reported, and the CDM reference runs through the same pipeline. While the gas is
tightly coupled, the DM is one fluid with c_eff² = γ n_v T/ρ_DM (γ = 1.05 in two-phase equilibrium, 5/3 after chemical
freeze-out). After momentum decoupling the objects are cold, and the gas is a collisionless Maxwellian solved exactly
with the Newtonian Volterra (Gilbert) equation.

Controls (`results/controls_vt.json`):
- A cold gas through the Gilbert branch reproduces CDM to 1.7×10⁻⁵, and to 4×10⁻⁶ at half step.
- The marching scheme against a pure ODE agrees to ≤ 3×10⁻⁴.
- A Maxwellian with the velocity dispersion of a 3 keV relic suppresses 20–30% less in k than the Viel fit
  (astro-ph/0501562). The "m_WDM-equivalent" column below is therefore an overestimate by roughly that much.
- A fluid gas in place of the Gilbert one changes ΔP by up to ×2 near k_fs.
- No diffusion damping is included. The oscillating tail beyond the first suppression (values > 1 at k ≳ 1/r_s) is not
  physical, and only the onset scales are used.

**(A) Codex's gas (f8186fa, sticking 0) [C].** Inputs: 6.595% of the charge free with v_rms = 0.0684 at 2 MeV,
elastic σ = 332/m_X² on Q = 2×10⁵ cores. Results:
- Γ/H = 1.2×10⁶ at 2 MeV; momentum decoupling at T_vis = 2.64 keV.
- Coupled-fluid sound horizon 3.7 pc; rms free-streaming after decoupling 25 pc (Codex's collisionless value: 40 pc).
- Equivalent thermal mass 1.9 MeV.
- P/P_CDM = 1.000005 (k = 100/Mpc), 0.99999 (10³), 0.9986 (10⁴), 0.986 (3×10⁴).

**M1 with this gas is CDM on every accessible scale** (|ΔP/P| < 2×10⁻⁵ at k ≤ 100/Mpc, at the integrator floor), which
confirms REVIEW_M1. This holds only if the cores never absorb or emit X. If they absorb X with any sticking above about
10⁻⁶, cold cores would swallow this gas (Γ/H ~ 10⁶·s), heat up and evaporate, which is case B.

**(B) The equilibrium vapour (NEW in this pass; conditional on tracks/fx_formation's thermodynamics) [C/E].**
fx_formation finds exchange rates ≫ H, and the objects survive only in its low-heat branch. In that branch the hidden
entropy S/Q sits in an X vapour at μ = Ω(Q) m_X. Then n_v(5/2 + Δ/T_h) = (S/Q) n_Q with Δ = 1 − Ω(Q). This is our own
reconstruction, and it agrees with fx_formation's table: T_h = 0.0109 and 1 − η = 0.053 at T_vis = 3.7 MeV, against
their 0.0105 and 0.052. Its consequences:
- T_h falls only logarithmically, from 0.014 m_X at 17 MeV to 0.002–0.008 m_X at freeze-out. Latent heat keeps the
  hidden sector hotter than the photons late on (T_h/T_vis = 0.85 at 1 MeV, 6.4 at 0.1 MeV).
- About (S/Q)/(5/2 + Δ/T_h) = 2–6% of the charge stays in vapour at v_rms ≈ 0.1–0.2c.
- It stays there until capture freezes out, Γ_c = n_obj·sπR²·v̄ = H, which happens at T_vis ∝ 1/ζ with
  ζ = s·Q^(−1/3).

The scan covered Q = 2×10⁵–10²⁰, S/Q = 1.235 (R0 heavy) and 2.07 (medium), and s = 1, 10⁻², 10⁻⁴; 72 histories and
48 transfer functions are in `results/summary_transfer.json`. Representative rows (S/Q = 1.235, then 2.07 in brackets):

| ζ = sQ^(−1/3) | example | T_vis,f | vapour mass fraction | m_WDM-equiv | r_s | λ_fs | P(10) | P(20) | P(50) |
|---|---|---|---|---|---|---|---|---|---|
| ≤ 2×10⁻⁶ | Q ≥ 10¹⁷, s = 1; most s ≤ 10⁻² | ≥ 44 keV | 4–35% | ≥ 50 keV | ≤ 0.05 kpc | ≤ 4 kpc | 1.0000 | 1.0000 | ≥ 0.9997 |
| 2×10⁻⁵ | Q = 10¹⁴ s = 1; Q = 10⁸ s = 10⁻² | 5 keV | 3.8% (6.4%) | 10 keV | 0.36 kpc | 29 kpc | 0.9998 | 0.999 | **0.991 (0.985)** |
| 1×10⁻⁴ | Q = 10¹² s = 1; Q = 10⁶ s = 10⁻² | 1.1 keV | 3.4% (5.8%) | 3.4–4 keV | 1.4–4 kpc | 63–104 kpc | 0.997 | 0.984 (0.973) | **0.934 (0.891)** |
| 1.7×10⁻⁴ | Q = 2×10⁵, s = 10⁻² | 0.56 keV | 2.8% (4.7%) | 2.7 keV | 7.7 kpc | 86 kpc | 0.995 | 0.978 | 0.86 (0.77) |
| 4.6×10⁻⁴ | Q = 10¹⁰, s = 1 | 0.25 keV | 3.1% (5.2%) | 1.2 keV | 5.6 kpc | 360 kpc | 0.962 (0.936) | 0.92 (0.87) | 0.81 (0.70) |
| 2×10⁻³ | Q = 10⁸, s = 1 | 62 eV | 2.8% (4.8%) | 0.43 keV | 20 kpc | 1.1 Mpc | 0.88 (0.81) | 0.71 (0.56) | 0.06 |
| 1.7×10⁻² | **Q = 2×10⁵, s = 1** | 6.7 eV (z ≈ 2.8×10⁴) | 2.2% (3.7%) | 0.11 keV | 180 kpc | 3.8 Mpc | **0.05 (0.41)** | oscillating | — |

For the comparison with data we take the thermal-WDM curves at the bounds, from the Viel fit [R]:
- 5.3 keV (Lyman-α; Iršič et al. 2017, arXiv:1702.01764): P(10, 20, 50) = 0.94, 0.74, 0.13.
- 6.5 keV (MW satellites; Nadler et al. 2021, arXiv:2008.00022): 0.96, 0.84, 0.27.

Requiring P_M1 ≥ P_WDM at k = 10–50 gives three regimes:
- ζ ≲ 5×10⁻⁴ is allowed.
- ζ ≈ 2×10⁻³ (P(10) = 0.81–0.88) is excluded.
- ζ ≳ 10⁻² is excluded outright: the whole DM oscillates acoustically below 0.2 Mpc.

The mixed-DM Lyman-α bound (Garcia-Gallego et al. 2025, arXiv:2504.06367: f_WDM < 0.16 at 2σ for 1 keV) does not
constrain a 3–6% warm part [R]. The suppression is that of a mixed cold + warm fluid with a dark-acoustic onset. In shape
it is the generic ETHOS-like or mixed-DM class (Cyr-Racine et al. 2016, arXiv:1512.05344; Boyarsky et al. 2009,
arXiv:0812.0010 [R]), so it is not distinctive in form. Its size, though, is computable from M1 inputs: S/Q from the heat
budget, Q from formation, s from the object's internal dynamics.

What this changes:
1. "M1 is observationally CDM" is a statement about ζ, not a property of the model.
2. The reference object Q = 2×10⁵ as the whole population, with s ≳ 0.1, would be ruled out by the Lyman-α forest
   [E, conditional].
3. fx_formation's ⟨Q⟩ = 10¹²–10²⁰ gives ζ ≤ 10⁻⁴, so the population ranges from CDM to a 7–11% deficit at k = 50/Mpc.
   That is the largest conservative-M1 signal, generic in kind, and at the edge of present Lyman-α and satellite
   sensitivity.

The deciding inputs are the sticking s and the formed ⟨Q⟩. Neither is computed anywhere yet. fx_formation's 48³
bath run would measure s.

## 3. Late-universe channels of conservative M1 (`code/late_and_portal.py` → `results/late_and_portal.json`)

- **X̄ [C/E].** The charge-exchange cross section of bf8ab3b is 2.76/m_X² at k = 1 [Q Codex, via m1_transport SKEPTIC].
  The channel is exothermic (2ω = 1.76 m_X), so σv stays constant at 1.95/m_X² down to low k [D]. With
  n_core(2 MeV) = 1.21×10⁻¹⁸ m_X³ this gives Γ/H = 1.2×10⁵ at 2 MeV, falling ∝ T to unity only at 16 eV. Other cases:
  - At fixed per-hit probability the burn-out scales as Q^(−1/3): Q = 10¹⁴ still has 150 e-folds.
  - In fx_formation's equilibrium the X̄ fraction is 10⁻⁵⁰ by 17 MeV.
  - So p_today ≈ 0 (consistent with fx_formation's "p ≈ 0"). A survival of p = 0.005 to halos would need a per-hit
    annihilation probability below 10⁻⁵.
  - At face value (p = 0.005, Sun, ρ = 0.4 GeV/cm³) it would give 6.2×10⁻²⁰ s⁻¹ per X̄ (2.7% per t_H), 0.163 GeV per
    event into one X at 0.93c, and 2.5×10⁻⁴ of the DM mass converted per t_H. All of it is invisible and escapes halos.
- **Mergers [C].** The rates per object per t_H for Q = 2×10⁵ are in the verdict table. Converted mass fractions:
  1.2×10⁻⁴ (Sun), 1.6×10⁻³ (inner kpc), 1.3×10⁻³ (cluster cores), 8×10⁻⁵ (dwarf cores). For larger Q, both σ/M and the
  release fraction fall. Decaying-DM-to-dark-radiation fits allow percent-level conversion (Poulin, Serpico & Lesgourgues
  2016, arXiv:1606.02073 [R]), and the products here are massive anyway. Dissipative-SIDM bounds need σ/M ~ 0.1–1 cm²/g
  [R]; we have 10⁻³.
- **Hidden radiation [C].** No massless or light hidden state exists:
  - lightest X 92.4 MeV, φ 226 MeV (it decays to XX̄);
  - the vector is heavy (m_A = 6eh) and decays to XX̄ well before BBN, with ε_BX = 0 so no SM injection;
  - Y and Σ are at ~10¹² GeV, and there is no Goldstone (the U(1) is gauged).
  The hidden entropy is Y_Q·S/Q = 6×10⁻⁹ of the visible. The vapour's kinetic energy gives ΔN_eff = 2.4×10⁻⁹ (1 MeV)
  and 2×10⁻⁸ (0.1 MeV), even though T_h > T_vis late on. Planck 2018 gives N_eff = 2.99 ± 0.17; Simons Observatory or
  CMB-S4-class σ ≈ 0.03–0.06 [R]. The answer is **zero hidden radiation.**

## 4. Minimal extensions and their size (`code/dd.py` → `results/dd.json`; obligations in `late_and_portal.json`)

**Physics [D].** A portal with source J(x) (b v_H(φ² − v²), A v_H(φ − v), or λ_HX v_H|X|²) gives a coherent object–nucleon
amplitude σ_n(q) = (μ_n²/π)(f_N m_N c F̃(q)/m_h²)², with F̃(q) = ∫J̃ e^{iqx} d³x. Its point limit reproduces the
standard Higgs-portal scalar formula.

**Rates [C].** The halo model is SHM (ρ = 0.3, v₀ = 238, v_esc = 544, v_E = 250.5 km/s), with Helm nuclear form factors
and f_N = 0.30 (0.26–0.35 [R]). For LZ we used the WS2024 numbers [Q]: 4.2 t·yr, 50% NR efficiency at 5.4 keV, S1c < 80
phd. The upper roll-off at 45/55/65 keV changes the results by < 0.5%. The calibration fixes N₉₀ = 1.86 at 40 GeV.
For point WIMPs the implied limits are 2.5×10⁻⁴⁷ (1 TeV) and 4.2×10⁻⁴⁶ cm² (17.3 TeV).

| Q = 2×10⁵ | Xe / Ge / Ar zero (keV) | σ_n0 per coupling² (cm²) | LZ 2024 limit | 10 events XLZD 200 t·yr | 10 events Ar 200 t·yr |
|---|---|---|---|---|---|
| b φ²\|H\|² | **27.56 / 48.80 / 90.27** | 1.40×10⁻²⁶ | b < 6.5×10⁻¹⁰ (σ_n0 < 6.0×10⁻⁴⁵, 14× the point limit) | b = 2.2×10⁻¹⁰ | 8.5×10⁻¹⁰ |
| A φ\|H\|² | 34.11 / 60.38 / 111.7 | 7.0×10⁻²⁷ GeV⁻² | A < 8.2×10⁻¹⁰ GeV | A = 2.7×10⁻¹⁰ GeV | 9.8×10⁻¹⁰ GeV |
| λ_HX \|X\|²\|H\|² | 41.68 / 73.79 / 136.5 | 5.3×10⁻²⁶ | λ < 2.7×10⁻¹⁰ | 9.0×10⁻¹¹ | 3.1×10⁻¹⁰ |

The zero moves as E_R ∝ Q^(−2/3): Xe 45 keV (10⁵), 9.7 keV (10⁶), 2.2 keV (10⁷), 0.5 keV (10⁸). Required couplings rise
slowly with Q (b_LZ = 6.9×10⁻¹⁰, 9.6×10⁻¹⁰, 2.2×10⁻⁹ at 10⁵, 10⁶, 10⁸). For ⟨Q⟩ ≥ 10⁸ every recoil lies beyond the first
zero, and the coherent advantage is lost.

**Is the signal distinctive? [C/E].** Asimov estimates, background-free:
- Xe alone: the composite is mimicked by a 28 GeV point WIMP, and 31 events separate them at 3σ.
- Xe plus Ar: 9.5 Xe events suffice. Argon sees 6% of the xenon count, so the target dependence E_R ∝ 1/m_A is the
  discriminator (Butcher, Kirk, Monroe & West 2017, arXiv:1610.01840; Hardy, Lasenby, March-Russell & West 2015,
  arXiv:1504.05419 [R]).
- The zero itself sits where the rate is already 10⁻³ of the peak (secondary maximum 1.06×10⁻³). It is visible only for
  a Q spread σ_lnQ ≲ 0.1 (contrast 8.8), and it is gone by σ_lnQ = 0.3 (0.9).
- Neutron-star capture through the portal loses coherence (|F/F₀|² = 10⁻⁸–10⁻¹³ at 0.3–0.5 GeV).

**Obligations [C/E].**
- φ²|H|² portal (θ = 0.0291 b), at the XLZD b = 2.2×10⁻¹⁰:
  - Γ(h→φφ) = 9.3×10⁻¹⁹ GeV, BR_inv = 2.3×10⁻¹⁶, far below ATLAS < 0.107 (arXiv:2301.10731) and CMS < 0.15 [Q/R].
  - Freeze-in Y_pair = 4.3×10⁻⁷, i.e. 85 pairs per charge, f_rest = 181. REVIEW_REHEATING's erasure obligation becomes
    99.994%. The pairs are in fact erased (objects absorb X̄), but their energy stays in the hidden sector.
  - Entropy injected at T_vis ≈ 20–150 GeV, where T_h = 0.0059 T_vis: ΔS/Q = 3.6×10⁴ (linear) to 3.3×10³ (saturated).
  - fx_formation's survival condition is S/Q < 3.47, leaving a budget of 2.2 (R0) or 1.4 (medium). Hence
    **b_max = 1.7×10⁻¹² / 1.4×10⁻¹², i.e. 4–6×10⁻⁴ XLZD events.**
  - A detectable φ² portal is therefore incompatible with the object population in M1. It needs a second hidden light
    state to carry the entropy, which brings its own ΔN_eff ledger.
- λ_HX: the same, with ΔS/Q ≈ 3×10³ at the XLZD value (scaled from h → XX̄).
- Trilinear A φ|H|²:
  - allowed because the ε tilt already breaks φ → −φ; A = 0 is stable in zero-portal M1;
  - there is no tree h → φφ (b_eff = θ·6(λ_φ − ε)v/v_H = 2.7×10⁻¹⁵);
  - φ is produced only through mixing, like a light Higgs-mixed scalar. Our IR freeze-in estimate (q g → q φ and
    q q̄ → g φ, effective Yukawa θm_q/v, Γ = c α_s y²T with c = 0.1–1) gives Y_φ = 4×10⁻¹²–4×10⁻¹¹ (top-dominated) and
    ΔS/Q = 0.4–4 at the XLZD value. This is marginal against 1.4–2.2.
  - The bound A ≲ 2–6×10⁻¹⁰ GeV allows about 5–50 XLZD events, uncertain by ×10 [E].
  - Higgs, beam-dump, supernova and BBN constraints on a 226 MeV scalar with θ ~ 10⁻¹¹ do not apply: φ decays invisibly
    to XX̄ in 10⁻²⁰ s, and the bounds start near θ ~ 10⁻⁷–10⁻⁶ [R].
- **Kinetic mixing [A].** The stated tree Lagrangian with real (rephased) κ and μ is even under hidden C (X, Y, Σ → c.c.,
  A → −A), and B_μνF_X^μν is odd, so ε_BX is not generated at any order (rp_portals item 9, confirmed). The Option 1A
  phase α breaks hidden C, but no field carries both charges. A non-zero ⟨J_Y J_X⟩ then needs a portal bridge with at
  least two insertions, and a C-odd (hence CP-odd) piece on the SM side, i.e. CKM at high loop order. With zero portals
  only gravity connects the sectors. Either way ε_BX is negligible.
- **Gravity [C/R].** Mass ranges from 3.1×10⁻²⁰ g (2×10⁵) to 1.3×10⁻⁵ g (10²⁰). Microlensing starts at ~10⁻¹² M_⊙
  (≈ 2×10²¹ g) (Niikura et al. 2019, arXiv:1701.02151 [R]). Gravitational direct detection (Carney et al. 2020,
  arXiv:1903.00492 [R]) targets ~M_Pl ≈ 2×10⁻⁵ g, so only the top of the fx_formation range qualifies, with a flux of
  0.3 m⁻² yr⁻¹. Capture in the Sun or neutron stars needs non-gravitational scattering, so the rate is zero without a
  portal.

## 5. Verdict on the observable status

- **[C/E] Conservative M1 is secluded dark matter.** No particle, collider, BBN, CMB-radiation, SIDM, lensing or capture
  signature exists: ΔN_eff ≲ 10⁻⁸, σ/M ≲ 10⁻³ cm²/g, X̄ burned out, merger output invisible.
- **Its only observable is gravitational clustering.** It is CDM on accessible scales (|ΔP/P| ≤ 1% at k ≤ 50/Mpc) when
  ζ = s·Q^(−1/3) ≲ 2×10⁻⁵. Codex's purely elastic gas meets this trivially (|ΔP/P| < 2×10⁻⁵).
- **The single most promising conservative signature** is the latent-heat vapour of §2, in the large-object branch with
  good sticking. For ⟨Q⟩ ≈ 10¹²–10¹³ and s ≈ 1 it gives a 7–11% deficit in P(k) at 50/Mpc (2–3% at 20/Mpc). That is
  just below present Lyman-α and satellite reach, and a target for higher-redshift Lyman-α and 21-cm work.
  - It is generic in shape and distinctive only in that its size is computable from M1 inputs.
  - It is conditional on fx_formation's equilibrium and on the uncomputed s and ⟨Q⟩ [E].
  - The reference population, Q = 2×10⁵ with s ≳ 0.1, would conflict with the Lyman-α forest [E].
- **The single most promising extended signature** is the pure-mixing trilinear A φ|H|² ≈ 3×10⁻¹⁰ GeV
  (θ ≈ 4×10⁻¹²) [E]. It needs ⟨Q⟩ ≲ 10⁶, σ_lnQ ≲ 0.1 and ζ ≲ 2×10⁻⁵ at the same time (so s ≲ 10⁻³ at Q = 2×10⁵), all
  against fx_formation's ripening. The signal would be a soft xenon spectrum mimicking a ~30 GeV WIMP, a deficit in
  argon, and with ≳10³ events a zero at 34.1 keV in Xe (60.4 Ge, 111.7 Ar).
  - The φ²|H|² portal Codex costed, and λ_HX, cannot be detectable without evaporating the objects
    (b_max ≈ 1.5×10⁻¹², 10⁴ below XLZD reach).

## Limits

- Section 2 is a linear toy: Newtonian, sub-horizon, smooth radiation, sudden decoupling, no diffusion damping.
- The vapour history is a non-relativistic two-phase model, started once T_h < 0.05 m_X. Equilibrium before that is
  assumed (fx_formation).
- Sticking s and ⟨Q⟩ are scanned, not computed. Momentum transfer is geometric (332/m_X² at 2×10⁵).
- Detector models are calibrated toys. Projections are 10 signal events without backgrounds or neutrino fog.
- Freeze-in for the trilinear is an order-of-magnitude estimate. The entropy budget 3.47 is fx_formation's [Q own].
- Literature numbers are [R] unless the text above says they were read (LZ 2024 PDF text: exposure, ROI, efficiency
  threshold, 2.2×10⁻⁴⁸ cm²; XLZD, Garcia-Gallego et al., ATLAS/CMS, Butcher et al., Hardy et al. via search abstracts).

## Reproduction (scratch $SP/wt2/fx_observables; OMP_NUM_THREADS=1; from this directory)

```
python3 -I code/objects.py                 # 110 s  -> results/objects.json, $SP/.../profiles.npz
python3 -I code/dd.py && python3 -I code/dd.py extra    # 100 s -> results/dd.json
python3 -I code/vapour_transfer.py A       # 2 s
python3 -I code/vapour_transfer.py 1e20,1e17; ... 1e14,1e12; ... 1e10,1e8; ... 1e6,2e5   # 55 s each, parts in $SP
python3 -I code/vapour_transfer.py merge && python3 -I code/summary.py && python3 -I code/fs_lengths.py
python3 -I code/controls_vt.py             # 60 s -> results/controls_vt.json
python3 -I code/late_and_portal.py         # 2 s  -> results/late_and_portal.json (reads fx_formation eq_history.json)
```
`results/claims.json` holds the verdict table in machine-readable form.
