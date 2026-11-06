# Replay X4

Use the published X4 branch with parent X3 checkpoint
`2b0316237a8f8a85b3db22453eb2dfb57231d1bc`. Commands below run from
the TestingGRUT root. They need Python 3 with NumPy and SciPy, GCC with
OpenMP, and the input files already in the checkout.

```bash
export OPENBLAS_NUM_THREADS=1
X4=research/reality_framework/grut_working_theory_v2_20261101/tracks/x4_closure
X4_WORK=$(mktemp -d /tmp/grut-x4-replay.XXXXXX)
python -B "$X4/code/trace.py" "$X4_WORK/trace" --threads 6
X4_TAPE="$X4_WORK/trace/runs/h55_trace/moments.txt"
python -B "$X4/code/trace_audit.py" "$X4_TAPE" "$X4_WORK/source_audit.json"
python -B "$X4/code/kinetic_replay.py" "$X4_WORK/kinetic" "$X4_TAPE" "$X4_WORK/kinetic.json" --workers 5
python -B "$X4/code/restart_check.py" "$X4_WORK/kinetic" "$X4_WORK/kinetic.json" "$X4_WORK/restart.json"
python -B "$X4/code/mean_modes.py" "$X4_WORK/modes" "$X4_TAPE" "$X4_WORK/mean_modes.json" --threads 7
python -B "$X4/code/verify.py"
```

The trace command can resume a budget-limited run with the same arguments.
It recreates the full timestep tape and records its SHA-256. OpenMP
reduction order can affect low-order digits. Use the source recurrence and
archive comparison in the audit to verify the run. Large scratch tapes,
C input files and state dumps can be regenerated with these commands.

The kinetic runner reads the four probability tables directly from the
hash-pinned X3 input. Its `--collect` option verifies the generated inputs
against existing inputs and then collects already finished engine runs.
The first collection stopped when the source hash lookup dropped the
decimal part of the half-step engine filename. All eight numerical
cases had finished. The corrected collector checked their stored inputs
byte for byte and loaded those completed outputs.

The archived result files in `results/` describe the executed cases.
The last command checks those files and their inherited inputs. To verify
a fresh replay as well, give it `--replay "$X4_WORK"`.
