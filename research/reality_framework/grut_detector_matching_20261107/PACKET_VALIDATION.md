# Standalone packet validation

The archive builder was executed from a clean local source checkpoint `7842051b37b9d1d291ed7df26e955960faa5d483`. The ZIP contained 147 tracked source/result files plus its manifest and instructions. It was extracted into a new directory; every SHA-256 matched. The documented local Git bootstrap was then performed. This local bootstrap is not the historical physics registration.

A new virtual environment was created outside the extracted tree. NumPy 2.3.5 and SciPy 1.17.0 were installed there from the exact requirements file; no source path from the original repository or pre-existing TeNPy directory was added. The following documented commands were executed from `/tmp` using only the extracted source and that environment:

| Command | Observed result |
|---|---|
| `verify_packet.py --output ...` | Source hashes, three reconstructed analyses, 128 distributions and D1 bounds verified; all scientific misses retained. |
| `bosonized.py --verify --output ...` | Vertex and interval normal-order checks reproduced. |
| `bosonized.py --mass 0.5 --cutoff 4 --output ...` | Maximum history TV from archive: 0.0 at stored precision. |
| `reduce_d1.py --output-dir ...` | All full-parent, dephased and reduced distributions reproduced with maximum TV 0.0; all gate decisions identical, including the misses. |

Every command exited zero and the extracted tree remained clean. D1 completed in 32.4 seconds under its 300-second cap. [Machine summary](packet_validation/summary.json), [commands and output](packet_validation/events.json), [replayed D1 result](packet_validation/d1_replay.json).

The archive being distributed also includes this validation record, a clarification of the analytic clock embedding, and an optional post-result field-moment diagnostic. None changes the frozen M0/M1/D1 production implementation or analysis. The final archive has its own manifest and accessible source commit. Integrity is checked again after building it; the expensive calculation is not rerun merely for documentation changes.

This validation establishes portability and reproducibility of the stated calculations. It does not improve the physical accuracy, establish independent external replication, or certify the research evaluator on unfamiliar theories. The original DMRG four-case and deliberate-failure checks are recorded separately in [VALIDATION_REPAIR.md](VALIDATION_REPAIR.md).
