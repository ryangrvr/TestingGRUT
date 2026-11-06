# Replay X11

Run from the TestingGRUT root with Python, NumPy, SciPy and GCC with
OpenMP. One BLAS thread per process gives repeatable controls.
Use a scratch directory for the full time histories and native spectra.
The original wide-cutoff midpoint failure is an expected outcome in
the specified case. The refined cutoff case continues separately.

~~~bash
export OPENBLAS_NUM_THREADS=1
X11=research/reality_framework/grut_working_theory_v2_20261101/tracks/x11_h55
RUN=$(mktemp -d /tmp/grut-x11.XXXXXX)
mkdir -p "$RUN/results"
python -B "$X11/code/production.py" "$RUN/production" --workers 3 --threads 2
python -B "$X11/code/axial_law.py" "$RUN/results/axial_law.json"

python -B "$X11/code/kinetic.py" "$RUN/kinetic/legacy_at055" "$RUN/production/reference" --legacy-expansion
for pair in at035:early_boundary at055:reference late:late homogeneous_late:late_homogeneous radial055:radial spatial055:spatial; do
    label=${pair%%:*}
    source=${pair#*:}
    python -B "$X11/code/kinetic.py" "$RUN/kinetic/$label" "$RUN/production/$source" --step-factor 0.125
    python -B "$X11/code/abundance.py" "$RUN/abundance/$label" "$RUN/production/$source" "$RUN/kinetic/$label"
done
python -B "$X11/code/kinetic.py" "$RUN/kinetic/at055_half" "$RUN/production/reference" --step-factor 0.0625
python -B "$X11/code/abundance.py" "$RUN/abundance/at055_half" "$RUN/production/reference" "$RUN/kinetic/at055_half"

python -B "$X11/code/kinetic.py" "$RUN/kinetic/at055_laterN2" "$RUN/production/reference" --step-factor 0.125 --later-n2
python -B "$X11/code/baryon.py" "$RUN/baryon/clock_later" "$RUN/kinetic/at055_laterN2" "$RUN/results/axial_law.json" inherited_parity_later_widths
python -B "$X11/code/baryon.py" "$RUN/baryon/native" "$RUN/kinetic/legacy_at055" "$RUN/results/axial_law.json" inherited_parity_native_widths
python -B "$X11/code/kinetic.py" "$RUN/legacy_n2" "$RUN/production/reference" --legacy-expansion --later-n2
python -B "$X11/code/baryon.py" "$RUN/baryon/later" "$RUN/legacy_n2" "$RUN/results/axial_law.json" inherited_parity_later_widths

python -B "$X11/code/collect.py" "$RUN" "$RUN/results"
python -B "$X11/code/verify.py" --results "$RUN/results"
~~~

The published archive includes both early output-only work integrals
and the later integration-step work ledger. The additional ledger in
a fresh replay leaves the evolution equations unchanged. The direct
homogeneous replay is executed by the collector. The two kinetic
clock settings and the alternate washout values are diagnostic
controls. A short evidence replay regenerates only the axial law
JSON and passes it to the archived verifier with `--axial-replay`.
[Production](results/production.json) · [Continuation](results/continuation.json)
