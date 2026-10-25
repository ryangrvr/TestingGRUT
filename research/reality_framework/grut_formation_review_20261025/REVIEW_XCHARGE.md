# Independent review of Codex 8986d34 (C* audit and X-charge breaking), with an Affleck–Dine window test

25 October 2026. Claude, issue #13. Three tracks wrote their own code and ran none of Codex's:
- `tracks/xcharge`: leakage rates;
- `tracks/affleck_dine`: one coefficient for both supply and survival;
- `tracks/cstar_lit`: the literature and Codex's corrections.

All are same-session AI agents, not outside replication. Every number is conditional on supplied fields and inputs.
None is a GRUT prediction, a validation or an exclusion.

## 1. X-charge leakage (xcharge)

Codex's rates are confirmed for the channel it computed: X⁵ 1.53600e-12 /s, X⁶ 3.91960e-50 /s. Two
forced-response solvers and a nonlinear time-domain run agree with them to 0.03–0.5%. The mechanism is confirmed
too: the first-order term averages to zero on a single-frequency core, and the secular loss is second order, ∝ |c|².

**Correction: pair emission dominates.** Expanding X^d to second order about the condensate emits anti-X pairs at
total energy (d−2)ω, which is above threshold. Back-to-back pairs avoid the edge suppression that limits single
quanta.
- The pair rate is larger by about 45× for d = 5 and about 2.7×10³× for d = 6 [E, free-Born, with two local-density
  variants agreeing within 10%].
- Total rates: 7.1e-11 /s (d = 5) and 1.04e-46 /s (d = 6).
- The |c₅| bound for under 1% loss in 13.8 Gyr tightens from 1.22e-4 to about 1.8e-5. It loosens to about 5.3e-5
  against the decaying-dark-matter limit of Poulin et al.
- X⁶ still passes at unit coefficient.
- The precedent is Hertzberg (PRD 82, 045022): quantum emission beats exponentially suppressed classical emission.

**Other results.**
- The formation era does not give a stronger bound than today's object.
- The required wormhole suppression is S ≥ 9.0 for the coherent channel alone and S ≥ 10.9 with pairs, against
  S ≈ 10 in Einstein gravity: marginal.
- A discrete gauge Z₆ (Krauss–Wilczek; Banks–Seiberg) removes X⁵ exactly.
- Codex's "1.23e-4" should read 1.2227e-4 (rounding).

**Open.** A Bogoliubov pair calculation in the medium, summed over partial waves; depleting profiles; a population
average.

## 2. Can one coefficient supply the charge and keep it? (affleck_dine)

**X⁶: conditionally yes.** A non-empty window for c₆ exists under a narrow set of supplied inputs. Examples with
instant reheating:

| Inputs | c₆ window |
|---|---|
| c_H = 1, H_inf = 4.7e13 GeV | [137, 1.6e3] |
| c_H = 100 | [0.010, 16] |

Within those inputs:
- The upper edge is the potential staying bounded during the Affleck–Dine epoch. Late survival never binds
  (c₆ ≤ 1.8e9 against 7.65e14).
- There is no window if T_R = 1e9 GeV, or if the Hubble-induced mass lasts until reheating.

**X⁵: no window** except a marginal edge at c_H = 100 and H_inf ≳ 2.4e13 GeV.

**New constraints** (all retrospective where they compare with data):
1. **Phase.** A random Affleck–Dine phase is excluded by CMB isocurvature, about 3300× over the Planck limit at
   c_H = 1. Only a heavy-angle, supplied-phase branch of Dine–Randall–Thomas type survives.
2. **Walls.** The large condensate pins φ = 0. When the condensate falls below X_c ≈ 1.13 GeV, φ picks its sign at
   random, forming Z₂ walls with σ^(1/3) = 505 MeV. That tension is excluded as a surviving network, so the walls
   must decay before about 1 MeV. That needs a sign bias ΔV ≳ 5.8e-26 GeV⁴, which is not in the action.
3. **Late fragmentation.** Onset comes at T ≈ 40–86 GeV, depending on the charge fraction. The estimated fragment
   sizes span 2.6e5 to about 1e15 charges, so the reference Q = 2e5 sits near the low end but is not selected.
4. **Charge-symmetric remainder.** Where the charge fraction is small, any surviving X X̄ or anti-objects could add up
   to 0.9(ε^(−4/3) − 1) of the calibrated matter. Unless hidden processes erase it, the window shrinks to a tuned
   point.

**Bearing on Codex's 7737487 gauge model.**
- The boundedness edge is the same physics as the induced radial term. Integrating out Y gives −|y|²|X|⁶/m_Y², a
  negative sextic, so the window needs a stabilizer or c₆ ≲ 0.3 λ_X (M/X_max)².
- The isocurvature, wall and remainder constraints apply unchanged.
- The pair channel raises the survival rates. X⁶ still passes.
- The gauged-Q-ball stability and kinetic-mixing items raised earlier remain to be checked at a stated g and h.

## 3. Literature and corrections (cstar_lit)

Codex's three corrections to Claude's earlier C* statements are right:
- the Hamiltonian constraint is not a measured zero sum;
- a one-dimensional fundamental space and e^(10¹²²) effective states need different spaces and a dictionary;
- fixed de Sitter capacity and a vacuum relaxing to zero are alternatives, and closure supplies neither.

All 41 arXiv identifiers and 5 DOIs check out.

**Refinements.**
1. **Observer dictionaries.** The published dictionaries disagree (Akers et al. 2503.09681). The Harlow–Usatyuk–Zhao
   rule gives at most e^(S_Ob) effective states, set by the observer's own entropy. All of its calculations use Λ < 0
   or 2D models. Sizes of order e^(S_dS) exist only in toy de Sitter models.
2. **Clause 5.** "The vacuum residue is a departure from equilibrium" is Klinkhamer–Volovik's own published proposal.
   It is prior-art conjecture, not a GRUT bet. Only its link to closure is GRUT's.
3. **Clause 1.** Semiclassical relational clocks are published: Chandrasekaran–Longo–Penington–Witten; Chen–Penington
   2406.02116; Kaplan–Marolf–Yu–Zhao 2410.00111.
4. **Missing from the audit.**
   - McNamara–Vafa 2004.06738.
   - The dissent of Antonini et al. 2507.10649 and Liu 2509.14327, and the opposing view of Engelhardt–Gesteau.
   - Nomura–Ugajin.
   - The covariant-bound proofs of Flanagan–Marolf–Wald and Bousso–Casini–Fisher–Maldacena.
5. **Z₆ global charge.** The total Z₆ charge must vanish modulo 6 on a compact space [A].

**Clause status.**

| Clause | Status |
|---|---|
| 1 | Semiclassically realized; nonperturbative 4D is conjecture |
| 2 | GRUT bet |
| 3 | 2D or Λ < 0 toy realizations only; contested |
| 4 | Conjecture; the capacity law is a missing GRUT bet |
| 5 | Klinkhamer–Volovik model-level; the link to closure is a GRUT bet |
