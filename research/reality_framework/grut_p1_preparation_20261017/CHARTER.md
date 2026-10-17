# GRUT-P1: controlled charge-balanced preparation

This is a separately named candidate and a bounded analytical construction,
not unchanged D1, an autonomous cosmology, or creation from nothing. The
design and elementary estimates were developed before this charter; the
charter freezes the primary protocol before executable controls. All new
derivations are author-derived, independent review pending.

## Primary architecture and exact question

Use **explicit post-formation control powered by a declared external work
source**. Start all new matter, electric and fuel registers in vacuum.
Record a macroscopic WM2 word, load fuel by a specified finite-duration
drive, convert fuel reversibly into opposite charges, then operate with the
conversion pulse off. Geometry and the actual receiver remain dynamical
throughout. No receiver reset, geometry clamp, favorable buffer measurement,
or one-way production channel is permitted.

The question is whether this protocol prepares extensive matter with fixed
positive electric energy, no root compensator, and a paid operational buffer.
Success establishes an externally powered preparation protocol, not an
autonomous matter-bearing phase or a distinctive prediction about nature.

## Frozen sequence and parameters

* N is a multiple of 64, N>=256. K=N/32 disjoint conversion bonds occupy
  recorded positions (0,1),(2,3),...,(2K-2,2K-1).
* Geometry: the WM2 rooted-word/storage/physical receiver law, g in [0.9,1.1],
  delta_i=bar_delta+epsilon_i, |bar_delta|<=1/2, |epsilon_i|<=1/4,
  L_N=64N^2, epsilon=1/4. New empty-head gates are explicit in MODEL.md.
* Two commuting hard-core species have charges +1 and -1, mass energy 1
  each. Every potential electric link has exactly three levels E=-1,0,1.
  Electric energy is +sum E^2 with coefficient one, not zero.
* Each potential unordered pair has a fuel qubit. Its energy is
  Omega_e=3+Delta_e, |Delta_e|<=1/2; reference Delta_e=0.
* Matter hopping/density interactions and conversion are off during formation.
  At fixed formation time T_f=4N^2, record the word. Accept m>=N/2; rejected
  trials count as zero. The original root weights are unchanged.
* Load the selected K fuel qubits using carrier frequency 3, amplitude
  eta_e in [0.9,1.1], for pi/2. Then convert for pi/2, with lambda_e in
  [0.9,1.1]. Both amplitudes are one in the reference instance. No pulse
  duration or frequency is retuned to detunings or outcomes.
* Then switch both directions of conversion off together, and switch on
  j_+=1/64, j_-=1/96 hopping and kappa=1/4 nearest-bond density interaction.
  Mass, positive electric energy, fuel energies and geometry remain on.
  Every switch and time-dependent drive has an energy/work account.
* Retain ell=N/4 positions for the comparison. Initial matter support is
  b=2K-1=N/16-1. Operation window: 0<=|t|<=N/65536.

## Gates, including dependencies

1. Specify registers, signs, all constraints, reverse terms and their
   invariance. No projected implementation may replace this calculation.
2. Derive a uniform conditional lower bound for P>=N/32 and full original-
   trial accounting. The original-trial target is >1/4 at N>=2^28, explicitly
   conditional on the archived, author-derived S9 formation theorem.
3. Bound electric energy and controller work by O(N), with fixed electric
   capacity. Explain the restriction this imposes on later charge separation.
4. Re-derive a complete-Hamiltonian finite-time comparison and information
   front with positive electric energy and matter backreaction. Do not reuse
   D1's matter-only tree isometry, which is not locally compatible with E^2.
5. Exhibit nonzero nearby transport for both dynamically created species.
6. Use bounded exact controls, not a new formation or compression campaign.

Stop if any structural gate fails. Do not change the primary parameters,
discard recombination, replace the receiver, relax an error gate after a
failure, or infer an all-size autonomous phase. Controls may use small
receivers only as explicitly identified algebra/implementation diagnostics.
No all-size claim may be inferred from those diagnostics. Ordinary floating
point checks are not outward-rounded arithmetic certificates.

## Preserved boundaries

D1 at 830746b remains the baseline. S9 formation and IG1 remain pending
independent mathematical validation; P1 cannot upgrade their status. The
path grammar, quantum framework, external time, relational generator rule,
electric cutoff and controller remain supplied. No dimension selection,
gravity, empirical distinction or autonomous preparation claim is authorized
by a successful P1 protocol.
