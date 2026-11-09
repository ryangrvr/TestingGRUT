# D1-01: prospective source-observable checks

The mathematical source family and physical calibration are fixed in DERIVATION.md before production. N22's completed numerical verdict, B2's separate failed N20 component, and P1's open parent obligations are historical and unchanged. D1 does not authorize a massive transfer.

## Physical experiment and new observables

Keep the historical L=5, m_R=0, Gaussian sigma=0.8, phase zero, coupling2 and B2 source functional F. Its corrected staggered flux is justified as the conjugate of the calibrated background source, using the exact rewritten Hamiltonian. The map is an operational matching prescription for generalized forces, not an exact operator equality between different Hilbert spaces.

Probe functions are g_m(z)=pi/(2L)sin(m pi z/L), m=1,2,3, and g_beta(z)=30(z/L)(1-z/L)^4/L. Their L1 norms are one; their signed integrals are 1,0,1/3,1. Fix v_j=integral_{I_j}g, with no quadratic reconstruction, and identity offset c_g=(1/4)integral_outer g. Match the lattice sum v_j y_j+c_g to the continuum sum (v_j/a)integral_{I_j}Phi+c_g. The outer term is calibration, not an assertion that the continuum field is constant in an outer cell.

Measure four sets of quantities, none selected using a D1 result:

1. Vacuum means, second moments and all cross-covariances of the four probes.
2. Full complex insertion matrices for G and G^2, and normalized moments in vacuum, W|Omega> and both K outcomes, observed at new times 0,0.43,1.31.
3. All four conditional branches of detector measurements at 0 and0.43, probed at1.31, retaining the 4x4 insertion matrices before normalization.
4. Mean injected energy in W|Omega> and both first-measurement branches, with an independent lattice hopping-cosine identity and the exact continuum Parseval comparator. Do not gate energy variance: it generally diverges for the continuum step-smearing kick.

Branch probabilities are retained as normalization data. No history-TV comparison is used as a physical D1 gate. Rare branches cannot be silently dropped. Require p>=0.01 for one measurement and p>=0.001 for two in both calculations.

## Frozen numerical decisions

Run N14,18,20,22 from the historical ground-state checkpoints named and hashed in INPUTS.json. Their reuse does not constitute a fresh independent ground-state solver. Refine N22 using its already archived tighter ground state and splitting each new propagation interval in two. Verify all input bases, hashes, normalization and ground residuals.

Use 8192 and16384 modes in an independently implemented Gaussian reference, retaining complex ordered covariances. A synthetic oscillator calculation at cutoffs24 and32 verifies both one- and two-backaction matrices and conditional moments, including a negative test for the kick-shift sign. An independent N6 operator construction checks mass-shift/source identities, ancilla calibration and commutators, with a missing-calibration negative control. No D1 production observable has been evaluated before the source freeze.

At both N22 solves, require every first/second moment and every complex insertion discrepancy, plus the fixed0.0001 numerical reserve, to be <=0.002. Require every vacuum covariance discrepancy plus that reserve <=0.001. Require relative kick and conditional-energy discrepancies, after adding1e-6 absolute numerical reserve, <=0.02. Require maximum N20-to-N22 moment/insertion change <=0.001, and N22 refinement change <=1e-8. These are prospective operational tolerances for these dimensionless probes, not universal operator-error bounds or fitted physical parameters. All gates must pass; no allowances are reassigned afterwards.

Require numerical residuals/norm checks<=1e-9, reference mode doubling<=1e-7, and the conservative conditional-moment tail estimates<=1e-4. The source-precheck threshold is1e-11. The machine-readable protocol.json is authoritative.

## Reference tail estimate

For each zero-extended step source, total variation B gives variance tail T=B^2L^2/(2 pi^4 M^2). Cross-covariance tails are at most sqrt(T_F T_G). For an ordered word with total absolute charge Q (4 for one kick and8 for two), let D=Q^2T_F/2, e_Z=D exp(D), e_L=Q sqrt(T_F T_G), L_0=|m_G|+Q sqrt(V_F V_G). Then conservative insertion errors are

    e_A=e_Z L_0+(1+e_Z)e_L,
    e_B=e_Z(L_0^2+V_G)+(1+e_Z)(2L_0 e_L+e_L^2+T_G).

The branch coefficient vectors have l1 norm1, so these bound unnormalized matrix contractions as well. For a normalized mean/second moment M with truncated branch probability p, use (e_A_or_B+|M|e_Z)/(p-e_Z), requiring p>e_Z. The means are integrated analytically. These analytic estimates are evaluated in floating point, not certified by interval arithmetic. Exact Parseval removes mode truncation from unitary mean kick energy. Conditional energy still uses p: bound its tail by E_W e_Z/[4p(p-e_Z)] and require this <=the1e-6 energy reserve. Also require same-size conditional/kick energy refinement <=1e-8.

## Budget and stop rule

One run at each of the four sizes, plus one N22 refinement. Cap each point at300 seconds and total production at900 seconds, one BLAS thread and4GiB committed memory per child. The inherited tested Windows Job Object runner enforces time/memory limits; job peak includes the numerical child behind the Python launcher. Preserve every exception or cap failure; no automatic retry, changed dictionary or larger size. Source and decisions are committed before production in a separate branch. No merge, release or Zenodo deposit.

## Reproduction and scope

Use Python3.12, numpy2.3.5 and scipy1.17.0 on Windows for the bounded runner, with input checkpoints matching INPUTS.json:

    python research/reality_framework/grut_source_dictionary_20261109/run.py --inputs INPUT_FOLDER --output NEW_EXTERNAL_FOLDER

The Gaussian reference is a separate implementation from the historical history code; the full lattice Hamiltonian is inherited. A pass supports this finite source dictionary on these observables and one/two-measurement states. It does not justify every possible detector history, a unique spatial interpolation, the endpoint extrapolation, a finite-duration switching approximation, or a rigorous continuum/ground-state error bound. If a derivation later changes those choices, it defines another experiment and leaves N22 unchanged.
