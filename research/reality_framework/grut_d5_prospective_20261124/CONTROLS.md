# Matched D5 interventions and failure cases

The central intervention takes place after exactly m=N/2 displayed
cells have been built on the ordered path. In every main control
F=2 flips the very same cell, A_1, and the untouched suffix stays
blank. For each original h,g and coin value, all seven arms have
identical damaged A, B, C and P immediately after the central
update. The selectors, corresponding head fields and backup
contents are the declared differences. The subsequent write
and maintenance rules are the same in all main arms.

| Central intervention | Main original whole-word probability | Pilot-null probability |
|---|---:|---:|
| Intact | 1 | 1 |
| Invert G | 0 | 0 |
| Scramble G | 1/2 | 1/2 |
| Disable G | 0 | 0 |
| Invert M | 0 | 1 |
| Scramble M | 1/2 | 1 |
| Disable M | 0 | 1 |

All selector arms have exactly the same initial energy law, same
Hamiltonian, same dimensions, same blank backups, and same coin
mixture. F and K are chosen in the preparation and stored explicitly.
No result is conditioned on the coin. For inversion (without disable),
the new whole word is different from the original. Gene inversion
changes its alternating versus uniform rule. Memory inversion
changes the phase of every cell. Disabling either register in
the main model places its content in its backup and yields blank
deposits. Scrambling is one controlled flip with an independent
fair coin. The complete state preserves that coin and the pilot.
Neither disabling nor scrambling erases the information globally.

The null's local deposition reads the retained seed, so failure
after scrambling an accessible record cannot be predicted from
the appearance of an intact structure or the existence of that
record. Gene controls serve as positive controls for the read
argument shared by both models.

## Three fault stages

| Stage | Occupied displayed cells on the canonical band | Fault | Authentication at endpoint |
|---|---:|---|---|
| Early growth | N/4 | Flip A_0 | B_0 stores the flipped value |
| Prefix complete | N/2 | Flip A_1 | B_1 stores the flipped value |
| Full construction | N | Flip A_0 | C_0 stores the flipped value |

The internally selected fault runs implement these in time-independent
hops. Before the first two faults, no full word is present on their
ordered branch. The full rebuild archives the damaged prefix and makes
the rest. For the third, the reverse pass maintains the earlier
completed output. The fixed overwrite rules perform their full
passes regardless of damage, so this is a regeneration test.

The separate lab pulses apply exactly the same local tape permutations
without a stage read or program trigger. Each run uses a single probe.
Stage probability and joint end-arrival, target and archive probability
are scored separately. The output-phase amplitudes are coherent,
including the tails that miss the expected stage. Neither successful
clock arrival nor its expected stage is used as a filter.

The small independent verification deliberately also uses L=64 at N=8.
The late probe in that diagnostic has only about 10.95% mass in its
intended clock band and about 55.14% joint end, tape and archive
success. Short-route performance is kept in the packet and is
outside the frozen large-clock success gate. It isolates why the
physical padding matters. The internal inverted and disabled
cases retain their zero original-target scores. The supplied finite
clock also has an exact recurrence to the initial blank data in
unperturbed trajectories at 2T.

[RESULTS.json](RESULTS.json) · [INDEPENDENT_CHECK.json](INDEPENDENT_CHECK.json)
