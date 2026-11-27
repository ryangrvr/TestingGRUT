# Reproduction

Start at the F4 draft PR head in a checkout of TestingGRUT. Run these commands from the repository root with Python 3 and Git. The calculation uses the Python standard library, exact integers and rational arithmetic. There is no simulator, fitter, network access or random seed in the replay.

```sh
git merge-base --is-ancestor 2f6a734a456c23272b8a413afa880eb4234f83f5 HEAD
git diff --exit-code 2f6a734a456c23272b8a413afa880eb4234f83f5 HEAD -- \
  research/reality_framework/grut_f4_charge_selection_20261126/FROZEN_PROTOCOL.md \
  research/reality_framework/grut_f4_charge_selection_20261126/PARAMETERS.json
cd research/reality_framework/grut_f4_charge_selection_20261126
python3 verify.py
```

The verifier replays [run_f4.py](run_f4.py) and [independent_check.py](independent_check.py). It compares each byte stream with its archive, then checks agreement on all 100 common fixture values and the exact empirical comparison scale. Production applies the $B-L$ shift and evaluates anomaly, operator, triangle, charge, and RG checks. The independent program solves from the weak equation and Yukawa constraints, compares the covering-group centers, traces representations independently and evaluates the same physical outputs through the conjugate singlet charges. The details of the all-integer proof and the spin-bordism input require analytical review. Five fixtures are implementation controls.

| File | SHA-256 |
| --- | --- |
| FROZEN_PROTOCOL.md | 884684ba85668314e80420c7f58810b9f1b212d3637e052ae046bd7a8006e669 |
| PARAMETERS.json | 0c6f086c6ebe04d2e1c7db0bc5fa4b0da6237b2b96802340b872d9804e0b855b |
| RESULTS.json | c64e0564eef9738ddfd12812c5c92839335be38a3d5e968c4c605fe6eb414974 |
| INDEPENDENT_RESULTS.json | 37cb6f31bf22ae01b621104e115f2bb4c4f18bacb2c55d44dc7cc4bafd7d12ed |
| INITIAL_RUN_FAILURES.json | 2da9eab3c2b583cad89fe377c33228e8ad886b5aa96eac35effd3961ed2366da |

Expected status is 69/69 production checks and 101 exact values in agreement with byte-identical replays. The first failed beta-function evaluation is retained as an instrument-control record. In that version the Weyl contribution was missing the common factor of three. The historical output is archival and the replay verifies the corrected arithmetic against the original frozen inputs. [Negative results](NEGATIVE_RESULTS.md) · [Internal criticism](INTERNAL_CRITICAL_REVIEW.md)
