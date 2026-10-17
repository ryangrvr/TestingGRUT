# P1 charge, energy, switching and resource ledger

No energy from absolute nothing is claimed. This primary architecture uses
an explicitly time-dependent classical work source after recording. A finite
autonomous quantum implementation of that controller is not constructed.

## Every account

| Account | Hamiltonian and status |
|---|---|
| Geometry/storage | 16(Q-m)+64 sum(q-a)^2+sum delta(q-R)-4(Q-1); inherited supplied law |
| Receiver | dGamma(4+h_L), frequencies [2,6], n_R=Q-1<=N-1; actual conditional state retained |
| Electric | +sum E_e^2, fixed coefficient one, three levels/link |
| Matter mass | +N_++N_-, fixed positive unit mass energy per excitation |
| Fuel | sum (3+Delta_e)n_fe, |Delta_e|<=1/2; starts in vacuum |
| Conversion | V_C plus its exact adjoint, consumes fuel and creates two charges and one link flux |
| Operating interactions | Gauge-covariant hopping and kappa nearest-bond density energy |
| Controller | Explicit load drive, conversion envelope and operating switch; work below |

The storage+receiver rewrite is exactly the original relative Hamiltonian
on Q-1=n_R. Receiver energy already released during formation is not reused
as fuel. All fuel energy comes from the subsequent external drive. Geometry
and receiver energy need not be separately conserved; their exchange is kept.

For any smooth interval d<H(t)>/dt=<partial_t H(t)>; at a sudden switch the
work is <H_after-H_before>, evaluated in the unchanged state. Thus

    E_final-E_initial=W_load+W_convert_switch+W_operate_switch. (12)

Between switches every declared time-independent total Hamiltonian conserves
its full energy. The loading drive's explicit time dependence supplies work.
Equation (12), including switch impulses, is the ledger, not separate
conservation of bare matter/fuel energies while an interaction is on.

## Exact per-bond reference ledger (not only scaling)

The initial word recording is also a supplied intervention. An explicit
ideal instrument uses a degenerate pointer, initially blank, with controlled
unitary sum_w Pi_w tensor X_(code(w)); X writes the word code reversibly.
Reading/dephasing that pointer produces the declared word instrument. All
pre-record matter/fuel/fields are vacuum. Every diagonal and receiver-only
term commutes with the word projectors, while nonselective recording removes
the assisted word-changing expectation. Its mean system-energy change is
-<H_assisted>, of magnitude <=g_max=1.1. Conditional record energies can
differ; the per-branch pulse ledger below uses each actual conditional state.
The blank pointer, O(N log N) record memory and ideal measurement/control
capability are supplied. No cyclic memory erasure or zero apparatus cost is
claimed. The record can be kept; its reset work is outside this one-shot
protocol. Ideal intervention time is an assumption, not a hardware result.

Write p_L from (6), p=p_L p_C. At the end of loading, fuel has mean energy
Omega p_L. In its rotating frame the conserved energy is initially zero,
so the drive expectation just before turn-off is -Delta p_L. The continuous
drive has delivered 3p_L; its switch-off adds Delta p_L. Total loading work
is therefore Omega p_L, without retuning the drive frequency.

The continuous drive's positive work is at most three per addressed bond:
over this interval the two-level population crosses at most its first peak.
Its positive switch-off contribution is at most 1/2. Thus the 3.5K load
bound also controls positive supplied system work, not just its signed net.
These are ideal system-energy transfers, not a bound on an apparatus's
wall-plug consumption, controller dissipation or memory erasure.

Conversion switch-on has zero mean work because matter is vacuum. During
the conversion pulse,

    E_mass+E_field =3p,
    E_fuel =Omega(p_L-p),
    <V_C>=Delta p.

These sum to Omega p_L. Switch-off work is -Delta p, leaving

    E_diag_after =3p+Omega(p_L-p)
                 =Omega p_L-Delta p <=7/2.                  (13)

At exact resonance p_L=p_C=1, all three input work units become two mass
units plus one positive electric unit per pair. There is no residual fuel
energy and no mean conversion-switch work. Away from resonance, unconverted
fuel, coherent vacuum/pair branches and work of the switches are all retained.

Turning on hopping has zero reference expectation: it changes just one
species' position, whereas differing prepared configurations change an
opposite-charge pair. Initial minus occupations determine the plus
occupations, and vice versa. Pair/vacuum coherences can survive tracing fuel;
they are NOT replaced by a classical many-body mixture. Turning on density
interactions costs

    W_density=kappa[sum_e p_e+sum_(adjacent matching pairs)e p_e p_(e+1)]
               <=kappa(2K-1)<K/2.                          (14)

This is additional supplied controller work, not a free interaction energy.
The exact reference energy increase at operation start is <=4K.

For the actual moving-geometry experiment, loading factorizes exactly and
costs <=3.5K. Conversion switch-on expectation is zero and switch-off work
has magnitude <=lambda_max K<=1.1K. At operation start N_+=N_-<=K, and

    ||V_op|| <=[2j_++2j_-+2kappa]K=(53/96)K.                 (15)

The density norm uses at most two opposite-species neighbors per carrier;
the hopping norm uses at most two moves per carrier and partial-shift norm
one. Thus the net additional work and the sum of positive switch/load work
are each bounded by 6K, independent of N. Including nonselective recording,
mean supplied system work per original trial is bounded by 1.1+6K (rejected
trials receive no subsequent pulses). There is no missing extensive root
charge. Actual-reference differences in an extensive energy can cost O(N e_p),
not just e_p; (12) and the bounds above do not require suppressing that term.

## Field profile and separation

Immediately after reference preparation a successful matching bond has E=+1;
all other fields are zero. Therefore sum E^2=X<=K and the two opposite
charges are compensated locally, not at the root. Their total charge is
exactly zero in every branch, not just on average.

For a single separated + and - on an otherwise empty path, the field is
one unit on the intervening string and zero elsewhere. Its electric energy
equals their recorded distance: separation costs energy, rather than being
free as in D1. With many charges the unique cumulative-charge fields can
cancel. Every physical state actually accessed still has |E|<=1 by the exact
link capacity, hence

    0<=H_E<=m-1<=N-1                                      (16)

as an operator on the physical space, throughout formation, preparation and
operation. Partial shifts forbid overflow; the full three-level link law is
retained. This is a **supplied capacity restriction**, not a claim that energy
alone derives it or that an unbounded electric theory has been simulated.

In the reference operation, mass and fuel numbers are separately conserved,
the density interaction is positive, and hopping norm is at most
2(j_++j_-)K=5K/96. Energy conservation further gives

    <H_E(t)> <=K+K/2+5K/96=(149/96)K                       (17)

for all operation times. Equation (17) is an expectation bound, not a
restriction to basis states below that energy. Complete-H comparison controls
its transfer to the physical experiment over the declared use window.

Charge neutrality alone would NOT suffice. If a capacity were allowed to
grow, K positive charges followed by K negative charges have

    sum E^2=K^2+2 sum_(r=1)^(K-1)r^2=K(2K^2+1)/3,           (18)

despite zero net charge. P1 avoids this D1-like cubic price through explicit
local production, fixed capacity, positive energy and a restricted operating
law. It does not produce an unconstrained charged plasma or 3D electromagnetism.

## Resources and remaining prices

* K=N/32 fuel registers are addressed, with O(N) work and two constant-time
  ideal pulses. Potential-pair electric/fuel registers still number O(N^2).
  Electric levels are constant, rather than D1's O(N) levels/link. No root
  countercharge/counter registers are required.
* A classical record of m labels, O(N log N) bits, and relational feed-forward
  select a matching and switches. Parallel local addressing and an ideal work
  source are supplied capabilities, not derived implementation costs.
* The receiver still has 64N^2 modes and <=N-1 quanta; positive receiver
  energy is <=6(N-1). Storage energy has O(N) norm, with explicit offsets
  in (2). Geometry-changing norm is <=2.1. The inherited formation time is
  O(N^2); operation lasts O(N) in the fixed declared units.
* Observational units require an externally calibrated energy E_0, with
  t_SI=(hbar/E_0)t. No conversion to observed masses, lengths or gravity is
  earned. The path dimension and common relational support rule are supplied.
* Keeping conversion on would permit fuel return and matter depletion. Turning
  it off is a physical intervention whose work is (13), not an autonomous
  stabilization theorem. The controls must test this reverse channel explicitly.
* The optional weak bond-source pulse used for response is a further supplied
  local intervention. In the retained-prefix channel its work is O(|theta|),
  bounded using only terms touching its two cells, not the extensive total
  norm. It inserts no particles. Readout calibration and apparatus reliability
  remain separate implementation obligations; no experiment has been performed.
