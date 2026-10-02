# NR-D — independent-code review of the S2-D-arrow key numbers

**Files:** `review/nr_darrow.py` (+ `.log`). No D-arrow code is used.

**Independent routines:**
- Krylov evolution (`expm_multiply`);
- n-general partial traces;
- an exact fresh-ancilla channel written separately from the joint simulation (the two agree);
- a QR product frame;
- a fresh H / S / CZ Clifford.

**Independence level:** *independent code path, not independent reviewer.* Each phenomenon was run at the original
parameters and at a variant set.

| ID | frozen claim | original number | reproduction (original params; variant) | grade |
|---|---|---|---|---|
| **NR-D1** | the exact time-reversal mirror; the reversed state runs an anti-arrow | S 1.238 → 0; mirror 8.9·10⁻¹⁶ | **1.238 → 0.000**, mirror 1.7·10⁻¹³ (Krylov tolerance); variant n = 8, tilt 1.1, S = {3}: 0.600 → 0.000, 9.8·10⁻¹⁴ | **REPRODUCED** |
| **NR-D3/D6** | same marginals: the correlated state recoheres, the product of marginals does not; bath ≈ maximally mixed | C 0.955 → 0; D ≥ 0.968; bath 6.902 bits | C **0.955 → 0.000**; D 0.955 → 0.968 (≥ 0.955); bath **6.902** / 7; marginals identical. Variant n = 7, τ = 2: C 0.999 → 0.000, D ≥ 0.705 (bath 5.445 / 6) | **REPRODUCED** |
| **NR-D4** | a fresh maximally mixed bath still gives Markovian contraction, an entropy rise and records; a correlated bath with the same marginals is non-Markovian | TD 0.482 … 0.048; S 0.633 → 0.997; record 0.225 | TD **0.482 … 0.048** monotone; S **0.633 → 0.997**; record **0.225**; the exact channel equals the joint simulation. Classically correlated bath: TD 0.482, 0.263, 0.170, **0.202**, 0.203 … (non-monotone) | **REPRODUCED** |
| **NR-D5** | bath reuse destroys semigroup behaviour; a single reused ancilla re-purifies | M = 1 → 0.007 | M = 1: non-monotone, S **0.633 → 0.007**; M = 2: non-monotone, min 0.126 | **REPRODUCED** |
| **NR-D7** | records form forward and un-form from the reversed boundary state | redundancy 0 → 6; reversed 6 → 0 | seed A: 0, 0, 0, 0, 4, 6 / reversed 6, 4, 0, 0, 0, 0 (MI 0.975 → 0.000); seed B: 0 → 6 / 6 → 0 (fresh couplings) | **REPRODUCED** |
| **NR-D8** | the same global state is "special" in one TPS and equilibrium-like in another | Householder frame 0 → rising; Clifford ≈ 1.9 flat | QR product frame **0.000 → 1.933**; fresh Clifford 1.91 – 1.99 flat; computational 0.93 – 1.26 | **REPRODUCED** |
| **NR-D2** | typicality: P(increase) ≈ ½, states already near maximal | 0.499 ± 0.013 | Δt = 1: **0.527 ± 0.013**; Δt = 0.3: **0.479 ± 0.013**; ⟨S_S(0)⟩ = 1.977 / 2. Both within ~2σ of ½; the exact value ½ is forced by the measure's U_t- and K-invariance (ΔS ↦ −ΔS) | **REPRODUCED** (statistical; ~2σ fluctuations noted) |

**S2-D-arrow numerical verdict: REVIEW-CONFIRMED.** The theorem-level scope errata (IR-03, IR-04) are separate and still
required.
