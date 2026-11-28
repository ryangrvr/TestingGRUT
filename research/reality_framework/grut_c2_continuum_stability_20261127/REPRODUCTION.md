# Replay the half-line certificate

Use Python 3.12, NumPy 2.3.5, SciPy 1.17.0, and mpmath 1.3.0. The retained run uses extended binary arithmetic with 63 long-double fraction bits and outward interval operations. The mpmath interval factor comparison uses 55 decimal digits. The acceptance bounds assume the host's default round-to-nearest elementary long-double operations and directed nextafter, plus the outward arithmetic of mpmath. The floating numerical tangent solver provides only the interpolation data. The verifier independently bounds the residual of the resulting exact polynomial.

Use the C1 profile archive in the sibling folder. Run from this directory with fresh output paths.

~~~sh
timeout 180 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 -u certify.py --output /tmp/grut-c2-results.json \
  > /tmp/grut-c2-execution.log
cmp RESULTS.json /tmp/grut-c2-results.json
cmp EXECUTION.log /tmp/grut-c2-execution.log
~~~

The corrected production and internal replay agree byte-for-byte in the retained environment. The protocol and C1 profile archive hashes are in the JSON. The first run, prior to the residual-factor amendment, can be replayed at its commit. Its failure is the widening of a propagated determinant interval. [First run](INITIAL_ATTEMPT.json) · [Amendment](AMENDMENT_01.md)

The production retains all three profiles from the pinned archive, rebuilt into exact nodal Hermite polynomials. The enclosure integrates their defects directly and uses the first accepted radial shift in the frozen list. The output and logs include every attempted shift. [Protocol](PROTOCOL.json) · [Results](RESULTS.json)
