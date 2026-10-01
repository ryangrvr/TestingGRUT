# S2-1 RESULT — system individuation from (Hilbert space, H) with no declared TPS

**Charter:** `probes/PROBE_CHARTERS.md` §S2-1 (pre-registered).

**Files:** `probes/S2-1/s2_1_individuation.py`, with log `probes/S2-1/s2_1_individuation.log`.

**Labels:**
- **EMERGENT BUT NONUNIQUE** (small systems) + **SYSTEM INDIVIDUATION SELECTOR, criterion-priced** (large systems);
- KNOWN RESULT IMPORT: Cotler–Penington–Ranard "Locality from the spectrum" (2019); Zanardi observable-induced TPS
  (2001) — SECONDARY, with the core dimension count reproduced ✓.

## 0. Verdict

> **(a) No-go baseline.** (Hilbert space, H) alone individuates nothing. Every TPS of a given total dimension is
> unitarily equivalent, so H fixes only its spectrum. *Any* individuation must come from an added **criterion**.
>
> **(b) The locality criterion selects, but only above a size threshold.** Criterion: "H is 2-local for some qubit
> TPS". The Jacobian rank of the map (2-local couplings) → (spectrum), modulo local unitaries, gives:
>
> | n | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
> |---|---|---|---|---|---|---|---|---|
> | fibre dimension | 20 | 39 | 59 | 72 | 62 | **0** | **0** | **0** |
>
> - **n ≤ 7:** the 2-local TPS is **non-unique**. Continuous families of inequivalent 2-local TPSs share the spectrum.
> - **n ≥ 8:** the 2-local TPS is **locally unique** (rank = P − 3n). There is no continuous ambiguity; discrete
>   ambiguities are not excluded by this test.
>
> **(c) Explicit hostile witnesses.** At n = 3 and n = 4, two 2-local Hamiltonians with **identical spectra**
> (difference 2·10⁻¹⁵ and 7·10⁻¹⁵) have **different** local-unitary invariants:
> - local-field norms {1.30, 2.03, 2.24} vs {0.75, 1.60, 1.88} at n = 3;
> - the pair-coupling singular values also differ.
>
> The same abstract H is therefore 2-local in two inequivalent tensor-product structures.

## 1. Structure inserted vs structure forced

| Item | Status |
|---|---|
| total dimension `2^n` and its factorization into **qubits** | **inserted**: the local dimension is a choice |
| the **locality criterion** (k = 2) | **inserted**: a criterion, not derived from H |
| uniqueness of the TPS **given** the criterion, for n ≥ 8 | **forced** (dimension count + rank; locally) |
| the fact that a 2-local TPS **exists** at all | **not forced**: for n ≥ 8 the 2-local image has dimension P − 3n < 2^n − 1, a Lebesgue-null set of spectra |

So a *generic* H has **no** local TPS. That locality exists is a contingent property of special Hamiltonians.

## 2. Measure ledger entry

"Generic spectra admit no 2-local TPS" is relative to **Lebesgue measure on spectra** (supplied). Under a measure
concentrated on local Hamiltonians the statement inverts. **MEASURE-PRICED.**

## 3. Information accounting

| Input | Output |
|---|---|
| spectrum (2^n − 1 numbers) + criterion (k-locality, qubit dimension) + a *promise* that the spectrum lies in the null set admitting a local TPS | a TPS (subsystem individuation), locally unique for n ≥ 8 |

- **Erased:** nothing. **Selected:** the TPS.
- **Price:** the criterion + the local dimension + the existence promise + (for the "generic" claim) a measure.
- **Relabelling of C5-A:** "which factorization" is converted into "which locality criterion + which local dimension
  + does H satisfy it".
- **The boundary moves partially.** For large systems, individuation is not an independent supplied datum *once a
  locality criterion is supplied*. The criterion is the irreducible remainder.

## 4. Hostile notes

- Local uniqueness from the Jacobian rank at one random point is generic in the coupling space (Gaussian couplings:
  a **measure-priced** sampling statement). Global uniqueness (no discrete alternative) is not tested.
- Level gaps shrink with n (1.5 → 3.7·10⁻³). Hellmann–Feynman stays valid because there are no exact degeneracies.
- The criterion could be replaced by others: k-local with a different k, geometric locality on a graph, or
  qudits instead of qubits. Each is another supplied choice. Whether the *criterion itself* can be selected is the
  open sub-question (S2-1b, queued).

**Status: S2-1 COMPLETE — No-go baseline confirmed. The locality criterion individuates subsystems uniquely (locally)
for n ≥ 8, but is non-unique for n ≤ 7. Criterion, local dimension and the existence of locality are supplied.**
