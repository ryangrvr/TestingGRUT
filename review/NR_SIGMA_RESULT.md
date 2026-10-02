# NR-Σ — independent-code review of the S2-Σ load-bearing claims

**Files:** `review/nr_sigma.py` (+ `.log`). No S2-Σ code is imported. It reuses the review routines from `nr_sigmah.py`:
- Möbius support weights;
- eigvalsh entropies;
- fresh H / S / CZ Cliffords.

**Independence level:** *independent code path, not independent reviewer.*

**Setup:**
- 202 groupings × 6 fresh frames = 1212 candidates;
- 4 Hamiltonians:
  - MFI chain (positive control);
  - translation-symmetric TFIM ring;
  - clustered pairs + weak bridges;
  - a **conflict-generating** H = MFI + (Clifford-rotated MFI), which is local in no single tested frame;
- 5 states: product, mid-spectrum eigenstate, Haar, record-forming, Gibbs β = 1 (20 runs);
- criteria (L, Q, R, P, M), plus graph sparsity G for NR-Σ-2.

**Historical counts (58/58) were not targeted**, by design.

**Review-code note (not a SCOUT-2 issue).**
- The first run showed a spurious H-only mismatch for one commutant frame (M only).
- Cause: Möbius inclusion–exclusion leaves ~10⁻¹⁰ cancellation noise on zero-weight supports, which an absolute 10⁻¹²
  "presence" threshold counted.
- Fixed with a relative threshold (10⁻⁹ × total weight). After the fix the commutant frames match the identity frame to
  2·10⁻¹⁶.
- The original S2-Σ uses Pauli transforms with an absolute threshold on coefficients and does not have this cancellation
  route.

## Grades

| # | historical S2-Σ claim | independent result | grade |
|---|---|---|---|
| 1 | no universal Pareto-dominating TPS in the tested families | **0 / 20** runs have a weakly dominant candidate, including the positive control under the full 5-criterion vector. Fronts hold 10 – 78 candidates spanning 3 – 10 local-dimension types | **REPRODUCED** |
| 2 | local d / grouping not selected | clustered model: L, P and Q are best at **16×4**; graph sparsity G and description length M are best at **qubits (2⁶)**. MFI: L and P are best at 2-cuts (32×2, 16×4, 8×8); G and M at qubits; Q at 16×4 | **REPRODUCED** |
| 3 | the scalar winner is objective-dependent | per run: 1 – 4 distinct winners over 4 monotone transforms × 3 normalizations; **3 – 21** Dirichlet-simplex winners (top share 0.23 – 0.88); **2 – 8** distinct winner sets over all 120 lexicographic orders. **The Pareto set is invariant under every monotone transform in 20 / 20 runs** | **REPRODUCED** — scope: *the tested natural criteria do not determine a winner without a further selection principle*; this does **not** prove no fundamental objective exists |
| 4 | the winner is scale-dependent | the equal-weight winner changes with the **time horizon τ** (clustered 4×4×2×2 → 4×4×4; MFI 8×4×2 → 4×4×4) and with the **perturbation strength ε** (both → 16×2×2 at ε = 0.3). **Stable window:** the redundancy threshold δ = 0.1 vs 0.3 changed no winner in any tested condition | **REPRODUCED (τ, ε), with a recorded stable window for δ** |
| 5 | the winner is state-dependent | for the same H, the **H-only** (L, P, M) winner is fixed (MFI: id 8×4×2; clustered: id 4×4×2×2), while the **H + ψ** winner changes across product / eigenstate / Haar / Gibbs (MFI: 8×4×2 → 4×4×4 → 16×2×2 in a Clifford frame) | **REPRODUCED** (H-only and H + ψ kept separate) |
| 6 | the supplied-locality (CPR-type) positive control succeeds | in the qubit class (n = 6, d = 2), the id frame and its commutant images tie at the top on (L, P, M) = (−1.3614, 0.6386, −0.01538). The fresh Cliffords (L ≤ −1.95) and the Haar frame (L = −4.47, M = −0.9998) are strictly dominated | **REPRODUCED** — positive control only; (n, d, k) is supplied, not derived |
| 7 | the commutant is gauge for H-only selection | e^{−iH·0.7} and e^{−iH·1.9}: H-only criteria identical to id, to 2.2·10⁻¹⁶ | **REPRODUCED** |
| 8 | epoch dependence scoped to exp(−iHs) | state criteria in the e^{−iHs} frame equal the id frame at ψ(−s), to ≤ 9·10⁻¹⁵ (and differ from the id frame at ψ by 1.2 – 2.4) | **REPRODUCED**. *Exploratory only:* a nonlinear commutant element e^{−i0.05H²} is gauge for H-only criteria (4·10⁻¹⁶), but **no** time shift s ∈ [−6, 6] reproduces its state criteria (min max\|ΔQ\| = 0.169). This is consistent with **not** generalizing the epoch statement to the full commutant. The full-commutant question remains **NOT ADJUDICATED** |

**Interpretive note (REPAIR 02 consistent).** Several scalar "winners" sit in the commutant frame. For state-sensitive
criteria this is the ψ(−s) time shift (claim 8), an epoch effect, not a physically distinct TPS. For H-only questions it
is gauge.

## S2-Σ review verdict

**REVIEW-CONFIRMED-WITH-SCOPE.**
- All eight load-bearing claims reproduce on an independent code path, with fresh frames and a different model set that
  includes a deliberately conflict-generating Hamiltonian.
- The scopes that must be carried are the ones already in REPAIR 02 / IR-01:
  - objective dependence is shown for the tested natural criteria only;
  - CPR-type selection needs (n, d, k) supplied;
  - the epoch statement holds only for the time-evolution subgroup.
- The δ-threshold stability window is a new, minor qualification of "scale-priced". Not every scale moves the winner.
