# S2-7 RESULT — consistent histories: does consistency select a unique quasi-classical set?

**Charter:** `probes/PROBE_CHARTERS.md` §S2-7. Pre-registered firewall: if selection needs a supplied
coarse-graining, split or time sequence → **A-PRICED**.

**Files:** `probes/S2-7/s2_7_histories.py`, with log `probes/S2-7/s2_7_histories.log`.

**Wall attacked:** A.

**Model.**
- A system qubit S plus m environment qubits, H = ω X_S + Σ g_k Z_S Y_k.
- Two-time histories: projectors onto a qubit basis at angle θ₁ (time t₁ = 1.3) and θ₂ (time t₂ = 2.9) on a chosen
  slot.
- Consistency: max off-diagonal |D| / √(p p′) < ε = 0.05.

**Labels:**
- **NONUNIQUE**;
- **A-PRICED** (slot, basis class, time sequence);
- quasi-classicality of the pointer set is **D-priced** (a large, persistent environment);
- KNOWN RESULT IMPORT: Griffiths / Gell-Mann–Hartle; Dowker–Kent — SECONDARY, with numerics ✓.

## 0. Verdict

**1. Consistency does not select a unique set.** Fraction of the 37×37 (θ₁, θ₂) grid that is consistent:

| model | slot S | slot E₁ |
|---|---|---|
| closed (g = 0) | 0.001 | 0.080 |
| monitored (g ~ 1) | 0.014 (19 sets, **all** > 0.15π from the pointer Z) | 0.037 |
| strongly monitored (g ~ 4) | 0.028 (39 sets, all away from Z) | 0.011 |

- Consistent sets exist on the **environment** slot as readily as on the system slot. Consistency does not say which
  slot is "the system".
- With a small environment, the Z–Z pointer set is **not** consistent: |D| = 0.25 at g ~ 1 and 0.24 at g ~ 4. The
  cause is recoherence in a small environment.

**2. Pointer consistency needs a large persistent environment.** With ω = 0.3 and g_k ~ U(2,6):

| m | 1 | 2 | 3 | 5 | 7 |
|---|---|---|---|---|---|
| Z–Z max\|D\| | 0.232 | 0.489 | 0.376 | 0.092 | 0.051 |

At m = 7 there are 42 consistent sets on a 13×13 grid, and the first-time bases among them cover **every** θ₁ on
the grid. Pointer consistency is approached only as the environment grows. Even then it is one consistent set among
many. Quasi-classicality needs a supplied D-item (a large, non-recurring environment) and a supplied A-item (which
set is used).

**3. The final-time basis is never constrained.** Off-diagonal terms with b ≠ b′ vanish identically, because the final
projectors are orthogonal. In every consistent set the last basis is a free choice. → **A**.

**4. Dowker–Kent-type sets.** Take projectors onto a randomly tilted *global* vector at t₁, and onto its evolved image
at t₂:
- max|D| ≈ 3·10⁻¹⁶ (exactly consistent) in 3/3 trials;
- this needs no TPS and no records.

Consistency admits non-local, non-record sets.

## 1. Structure inserted vs forced

| Item | Status |
|---|---|
| consistency of a given set | **forced** (computed from H, ψ₀) |
| which set (slot, basis class, time sequence) | **inserted: A** |
| quasi-classicality of the pointer set | **needs D** (a large persistent environment) **+ A** (choosing it) |
| the final-time basis | **free** in every set |

## 2. Accounting against T2-1′

Consistency is a **filter, not a selector**. Each consistent set it passes is a possible "what is recorded"; it does
not choose one. Choosing the quasi-classical family requires:
- a split (which slot), shared with S2-8 and S2-1b;
- a basis / coarse-graining class;
- a time sequence.

All of these are A-data. **A is not eliminated.**

**Status: S2-7 COMPLETE — NONUNIQUE; A-PRICED.**
- Consistent sets abound on system and environment slots.
- The pointer set becomes consistent only with a large environment (m = 7: |D| = 0.051).
- The final-time basis is always free.
- Exactly consistent global (Dowker–Kent-type) sets exist.
