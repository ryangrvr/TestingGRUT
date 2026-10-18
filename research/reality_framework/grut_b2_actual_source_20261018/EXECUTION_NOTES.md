# What ran, and what the arithmetic proves

Charter/dependency freeze: 21ef4bd618937f404a7683fd8cdbb11c9adb5463.
Actual proof/code freeze: ce4664cd94cd40ef972c984fc4a060976a1e80d3.
All eight files in EXECUTED_SOURCES.json were hashed before primary
execution and unchanged through replay. No fallback or extra variants.

One Python 3.12.14 primary child and one unchanged replay executed
controls.py, standard library only. The Q class uses exact pairs of
Fractions; Gaussian elimination, block inversions, coordinate assembly,
site/fuel identities and all probability margins use exact arithmetic.
There is no floating point spectral or operating-probability calculation.

The inherited domains are Gamma(N,2), N=0,1,2, with original diagonal and
root counterterm only. At z=0 and z=i, the code independently assembles
the original coordinate matrix and the block matrix, checks equality,
checks finite-depth Weyl Schur inversions, and reconstructs every co-pair
site resolvent and fuel resolvent. It detects that a single stationary
tail substituted at finite depth is not exact. These compressions are
not ended physical paths or replacements for the infinite source.

The error-envelope check uses positive 64-term rational Taylor lower
bounds for exponentials at worst-case L=2016. Monotonicity supplies the
all-L scalar envelope; no length or time grid is used. It verifies only
the probability composition IF the missing source inequalities hold.

| Worker | Child wall seconds | Peak child RSS KiB | Result |
|---|---:|---:|---|
| Primary | 0.37693205300001864 | 11136 | 45 checks, six resolvent cases pass |
| Replay | 0.3620585740000024 | 11136 | Byte-identical output |

Combined child wall time: 0.738990627... seconds. This includes child
startup/output, not assistant derivation, delegated challenge or connector
publication. RSS is the operating system's child high-water mark.

CONTROLS.json SHA-256:
`47d6dc5d99a0b6121a7705956faf3c7568ba8e04435985a143d71fd8aa92fdf6`.
Logs record stdout/stderr, exit status, version, wall time and RSS.
MANIFEST.json covers every final file except itself; verify.py checks
its hashes, executed sources, resource bounds and output agreement.

A separate delegated Challenger read the actual operator and target,
then exchanged formulas with the Builder. Their agreement has shared
repository and derivation ancestry. No external reviewer was contacted.
No physical source dynamics, spectral-type classification, late-window
probability, measured experiment, paid compute, installation, persistent
service or background job ran. No independent external validation occurred.
