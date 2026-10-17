# WM2-S9 — a bounded methods test and an analytic operational regime

**Completed S9: an analytic operational regime is derived; the bounded finite recovery fails.**
The unchanged model admits a Builder-derived, nonshrinking-disorder formation
bound at sufficiently large size. The finite threshold is extremely conservative:

    N>=2^28, L=64N^2, 4N^2<=t<=8N^2
      => Pr(m>=ceil(N/2))>1/2.

This is a theorem derivation about a supplied-resource quantum path model,
not external mathematical validation or a distinctive prediction about nature.
It does not settle five/six-event calculations or the large intervening
finite-size range. It also does not establish permanent geometry or spacetime.

## What changed physically

No Hamiltonian, neutral preparation, disorder class or receiver resource
changes. S9 addresses the missing history dependence analytically. Exact
append identities and the receiver equation imply contact occupation O(1/N).
Both noncommuting creation and absorption terms enter a finite-time current
estimate. The error is explicitly bounded back to the complete closed law:

    K'(t)=(g/N)^2 sum_v sqrt(4-delta_v^2)<P_(v free)>+r_N(t),
    |r_N(t)|<=56g^2 N^(-5/4)+10g^2/(N sqrt(t))
                 +16g^2 N*2^(-120N^2).

The positive receiver spectral density covers every allowed detuning. A
Gronwall bound on the expected remaining-free count yields formation.
This is a controlled first-moment comparison, **not** an independently
emitting classical process, a reduced-state Lindblad theorem or a bath reset.
See FORMATION_ANALYSIS.md for the complete short proof and conservative constants.

## One integrated physical experiment

Combine that formation result with S4's receiver-uniform conditional theorem:
prepare, evolve, accept all macroscopic records, record the path, retain its
actual receiver, inject the neutral adjacent-root-pair probe, then propagate
with the **full** Hamiltonian. Rejects give zero signal. Recording/probe
injection and their energy/accounting remain supplied interventions.

The joint event of acceptance, consistency and retaining an original N/4
prefix has exponentially small conditional loss through an O(N) window.
Nearby phase response is nonzero. A separate complete-Dyson comparison to a
half-line probe proves a weak but nonvanishing injection-arrival probability
at recorded distance proportional to N, with no deletion of physical paths.
The full original-trial joint lower bound is 1/589824. It is a receiver-set
arrival event, not a phase-response lower bound at every distant site.
See INTEGRATED_PROTOCOL.md and BALLISTIC_TRANSPORT.md.

At the proof threshold the receiver has 2^62 modes and formation starts
at 2^58 hopping time units. Its protected use window is O(N), versus
O(N^2) assembly; their ratio shrinks. This is an asymptotic operational
relational-geometry protocol with stated resources, not an equilibrium phase.

## Bounded numerical construction

Freeze `85fc1bd` predates controls and primaries. The production sweep
certificate is implemented: nested ranges, direct in-range residual,
unknown-prefix-Gram correction, square sums inside eligible monotone
sweeps and linear sums across reversals/returns. Full physical controls
pass; the 333,505-amplitude N=3/L=576/T=36 error is .01359811 inside
.19647972. The arithmetic claim retains the declared IEEE/BLAS assumptions.

Exactly the two original five-event root components were run at T=100,
with rank 96 and a predeclared 1200 seconds each. Full original-weight
recovery fails. No fallback or N=6 trajectory is launched; this numerical
branch is closed. The worker limit differs from S8's 900 and S7's 600;
no equal-budget or equal-error superiority claim is made.

| Remaining counts | Status | Completed / target model time | Wall s | Peak MiB |
|---|---|---:|---:|---:|
| [2, 2] | UNRESOLVED | 65.032680 / 100 | 1202.430 | 68.305 |
| [3, 1] | UNRESOLVED | 76.746725 / 100 | 1208.667 | 59.840 |

All stop reasons and distinct error components are in RESULT_TABLES.md and
the raw JSON. The recorded stops are runtime limits, not exhausted
projection/error allowances. Neither is a negative formation observation.
The production certificate passes its controls but its evaluated cost
does not complete this finite recovery protocol. Another solver revision
is not opened. The analytic theorem does not apply at N=5 and cannot
substitute for the missing target-time numerical result.

## Validation and remaining architecture

FORMATION_CONTROLS.json checks full tensor Heisenberg identities, the
detuning-cancelled derivative, nonzero creation commutator, closed-walk
moments, spectral rates and exact rational final constants. TRANSPORT_CONTROLS
checks the finite support comparison and reference-walk identities. They
do not numerically establish all-N theorems or independently validate proofs.

Keep the three ledgers separate. Selected dimension, original charged-matter
universality, a common continuum causal structure, gravitational dynamics
and empirical discrimination remain open. Quantum structure, time, path
grammar, receiver resources and interventions are supplied. Standard quantum
mechanics contains the entire construction. OBSERVABLES.md states its earned
conditional constraints and rejection criteria without claiming otherwise.
The disorder theorem covers independent **storage** detunings and a common
g interval; arbitrary matter onsite disorder or receiver imperfections are
outside its scope. The limit is a sequence with supplied coupling g/N,
not a conventional fixed-local-coupling equilibrium thermodynamic phase.
