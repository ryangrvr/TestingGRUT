# NR-ΣH — independent-code reproduction of the S2-ΣH key claims

**Files:** `review/nr_sigmah.py` (+ `.log`), `review/nr2_diag.py` (+ `.log`), `review/nr2_certified.py` (+ `.log`).

**Independence level.**
- Same agent, **independent code path**: no shared code with `probes/S2-Sigma` or `probes/S2-SigmaH`.
- The one exception is that `nr2_diag.py` imports the original module **only** to rebuild the original conflict
  state for diagnosis.
- Different numerical routines (README); fresh random frames over three seeds.
- Sanity check: the Möbius support weights sum to ‖H − tr‖² exactly (885.478880 both ways), with maximum support size 2.

## Results

| ID | claim (S2-ΣH) | reproduction | status |
|---|---|---|---|
| **NR-ΣH-1** | positive compatible case: a dominant set exists and, after the LU quotient, lies in the native (id) frame | 3/3 seeds: dominant set yes; front 6 – 7 bare → **5 LU classes, frame id only**; rules A and B agree on the frame (extra labels are LU copies) | **REPRODUCED** |
| **NR-ΣH-2** (as originally constructed) | conflict "H local in Σ_A, ψ product in Σ_B": no dominance; both frames on the front; rules A and B disagree | fresh weak-conflict states: 0/3 matched the full claim (seed 11: no dominance, but the front holds only the id frame; seeds 22, 33: dominance in id) | **NOT REPRODUCED as constructed → diagnosed below (IR-05)** |
| **NR-ΣH-2-certified** | the same claim for a certified-incompatible pair (in the H-local frame, every grouping has orbit-minimum C ≥ 0.5 bit) | 3/3 seeds (id-frame min C_min = 1.67 – 1.68 bits): **no dominance; the id and conflict frames are both on the front; rule A → id, rule B → conflict frame** | **REPRODUCED** |
| **NR-ΣH-3** | Haar: no dominance; the correlation-first rule returns the trivial ψ-product frame | 3/3: no dominance; rule B → W_ψ (QR-built) | **REPRODUCED** |
| **NR-ΣH-4** | local-dimension tie in the compatible case | 3/3: the dominant LU classes span 32×2, 16×4, 8×8 with identical (L, C_min) | **REPRODUCED** |
| **NR-ΣH-5** | epoch: C₀ loses dominance under a time shift; the C_min dominant set is unchanged, with t* shifted | 3/3: C₀ dominance yes → NONE after s = 2.25; the C_min set is identical; t* 0.00 → −2.25 | **REPRODUCED** |

## IR-05 diagnosis: the original H1 conflict instance was not certified incompatible

`nr2_diag.log`: in the H-local id frame, the **original** S2-ΣH conflict state CLIFFS[0]·ψ_prod is **exactly product
across the 2-cut (01235)|(4)**, with C_min = 0 at t = 0. The fresh weak-conflict states are likewise product across some
2-cut:

| seed | product cut | effect |
|---|---|---|
| 11 | (1)\|rest (inner) | no dominance |
| 22 | rest\|(5) (end) | dominance in id |
| 33 | (0)\|rest (end) | dominance in id |

**Consequences:**

1. The original H1 non-dominance was real. But it came mainly from a **trade-off between cuts inside the id frame**:
   the inner product cut (01235)|(4) has C = 0 but crosses two bonds, while the end cuts cross one bond but have C > 0.
   It was **not** a clean "H-local frame vs ψ-product frame" conflict. The pair was **coarse-grouping-compatible**.
2. When a weak-conflict state's product cut coincides with H's locality-optimal end cut, the pair is compatible and
   **dominance correctly appears**. This is consistent with the dichotomy, not a counterexample to it.
3. The **intended** claim, tested on certified-incompatible pairs, **reproduces 3/3**.

**Net:** the S2-ΣH qualitative conclusion stands. Its H1 *description* needs an erratum (E-06). The compatibility
condition must be checked **at every grouping**, not only at the finest qubit TPS. That makes "compatibility" itself
grouping-relative, which further supports the Σ ⊗ H_corr coupling.
