# A measurable interface, with parameter discipline

This specifies tests of the pinned channel class, not measured evidence or
a prediction of the unknown large-size populations. The full Hamiltonian,
initial receiver, fuel preparation, output event and time must be declared
before collecting data. Hardware feasibility is not established here.

## Operational quantities without environmental decoding

Fix N,w,t and the same conditional receiver preparation for every trial.
Use the existing final success flag and the two decoded color ports A,B.
Count failed trials as zero for SUCCESS-WEIGHTED expectations; do not
renormalize each observable before composing them. The actual receiver
state is retained during evolution, not reset to a favorable packet.

Three input preparations suffice to identify the success instrument under
its proved identity/SWAP structure:

| Prepared input | Measured original-trial quantity | Model value |
|---|---|---|
| $\lvert00\rangle$, a triplet | Success frequency | u |
| $(\lvert01\rangle-\lvert10\rangle)/\sqrt2$, singlet | Success frequency | v |
| $\lvert01\rangle$ | Success-weighted (Z_A-Z_B)/2 | Re(c) |
| $\lvert01\rangle$ | Success-weighted (Y_A X_B-X_A Y_B)/2 | Im(c) |

The last row uses separately measured local Pauli correlations; it does
not read fuel or geometry records. Its sign follows from the convention
c=<chi_-,chi_+>. Statistical errors and apparatus errors would have to be
propagated jointly in a real experiment; no sample size or empirical
confidence interval is supplied by this packet. Finite input tomography
does not certify a size/time continuum.

For a direct sign check, use t=(01+10)/sqrt2 and s=(01-10)/sqrt2.
The successful output of input01 is [[u,c],[c*,v]]/2 in the t,s basis.
The displayed two observables are sigma_x and -sigma_y in that basis,
which give Re(c) and Im(c). This interface derivation is analytic; it is
not an additional physical run or part of the64 frozen execution checks.

Likewise retention frequencies on the triplet and singlet preparations
give u_T,v_T using the SAME retained-prefix/particle-count measurement.
Output cannot substitute for retention at its higher threshold.

## Held-out consequences and a rejection test

Once the four real instrument coefficients are calibrated at that fixed
setting, leave them fixed. For any other prepared input rho, predict
p_success=u Tr(P_+rho)+v Tr(P_-rho) and the complete unnormalized color
state from the B2-IS operation. Its reference extension is fixed too.
For any reference-entangled triplet qutrit input, predict success u and
exact conditional identity. These are falsifiable held-out-input relations.
One must compare measured departures with independently specified preparation,
measurement and symmetry-breaking errors, rather than refit each input.

An observed output outside span{identity,SWAP} or a reproducible triplet
logical error beyond that account rejects this color-blind candidate or
its claimed implementation. A measured small success probability rejects
an operating task even if successful logical states are perfect. Conversely
the finite diagnostic's nonzero perfect-recovery deficit does not alone
reject the1% task.

Instrument calibration is NOT a derivation of u,v,c from microscopic g,
delta,jp,jm,lambda,Omega. For a predictive MICROSCOPIC test, calibrate those
microscopic inputs independently, then compute these functions with controlled
errors and compare without further adjustment. The missing large-size moving
calculation still prevents that numerical prediction here. The exact finite
coefficient test supplies only the explicitly described short-time diagnostic.

## Comparison and assumption price

Ordinary quantum mechanics with a color-blind exchange-coupled environment
has the same instrument relations and noiseless code. This is standard
interaction-algebra physics, not a distinctive GRUT signature. Testing it
could validate a specified engineered channel, not the ontology of reality.

The quantum framework, supplied relational interaction rule, path grammar,
receiver resources, prepared fuel/input entanglement, record and final
decoding access remain inputs. MC-E1-T narrows the alphabet; an optimal
parity phase correction or environment feed-forward would further change
the decoder and must be separately authorized/frozen. No such correction
was inserted into the original operating theorem in this cycle.
