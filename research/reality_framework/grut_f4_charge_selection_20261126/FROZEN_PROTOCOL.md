# GRUT-F4 frozen analytical test

Freeze date 2026-11-26 (America/Chicago). Parent is [draft PR #26](https://github.com/ryangrvr/TestingGRUT/pull/26) at 95ad5c0f158b128d337da8fe862894015b5d8f8c. Accept F3 Verdict C. The freeze covers this document and [PARAMETERS.json](PARAMETERS.json). The published experimental input is a known benchmark at freeze. No production calculation, alternative principle or parameter search precedes or changes the test.

## Candidate

Would forbidding exact connected ordinary internal phase symmetries of the light local matter theory, modulo gauge transformations, uniquely fix the electric charge direction among the gauge representations allowed by the other assumptions below? The physical motivation is the absence of exact internal global symmetries in holographic quantum gravity under stated assumptions. We test its ordinary continuous phase consequence as a low-energy selection criterion. The extrapolation to all quantum gravity is conjectural. [Harlow and Ooguri](https://doi.org/10.1103/PhysRevLett.122.191601)

## Model class and assumptions

- Four-dimensional spin spacetime, Lorentzian local perturbative QFT and a Hermitian Wilsonian EFT about a stable electroweak vacuum. Gravity may serve as a background for the mixed anomaly test. A nonperturbative quantum-gravity completion is outside the test.
- Fix the compact group $(SU(3)\times SU(2)\times U(1)_y)/\mathbb Z_6$, with the generator $(e^{2\pi i/3}\mathbf 1_3,-\mathbf 1_2,e^{i\pi/3})$. Use integral weights $y$ on this specified circle. Supply one Higgs doublet of positive weight 3, with the usual neutral component taking a vev. The Higgs weight and global form are assumed.
- Supply three families with identical gauge representations, the usual quark and lepton nonabelian representations, and one additional left-handed singlet $n^c$ per family. The heavy or light singlet spectrum and neutrino mass mechanism of nature are unresolved. Use generic full-rank Yukawa matrices for all four Dirac couplings with common commutant a scalar on each left-handed flavor space. The exact phase test is on the specified EFT, including its symmetry-violating perturbations.
- Allow Lorentz- and gauge-invariant local polynomials of arbitrary finite engineering dimension, with covariant derivatives and field strengths, and no symmetry-breaking external spurions. Work order by order below each EFT cutoff with perturbatively small couplings. Test no exact ordinary phase symmetry and do not impose absence of discrete, higher-form or emergent symmetries. The perturbative local-operator class supplies the observable dimension grading.
- Cancel perturbative gauge and mixed gravitational anomalies and use the standard spin-bordism global gauge anomaly test for the specified group. Fix common naming, orientation and the normalization convention $Q(W^+)=1$. Compare charges as ratios to an electron or a charged weak boson. Test no equality of effective propagation speeds.

## Proposition and controls

Solve the Yukawa and anomaly constraints and the quotient condition with the stated matter content. Represent any remaining freedom by the integer $k=(y_Q-1)/6$ and the reference weights in the input table. Apply the candidate to this whole family. The decisive obstruction is at least two physically inequivalent allowed members with no exact ordinary continuous phase symmetry. To support the stronger all-integer result, require an explicit nonvanishing invariant that breaks $B-L$ at every $k$, an independent $B+L$-violating direction, and proof of all allowed continuous phases. Track the weak and mixed anomalies of the residual phase at the renormalizable level.

The witnesses to test are the gauge-invariant $QQQL$ contraction and $n^cn^c$ dressed with powers of $H^T i\sigma_2 D_\mu H$ and its conjugate. For the dimension bound, test the center-triality and Lorentz parity rules, and the triangle inequalities with the two slopes in the input file. Include the unbroken-phase engineering dimension, gauge-invariant Weinberg-operator control, and the no-$B-L$-violator control. Retain the hydrogen-atom neutrality degeneracy and the dimension-truncated ambiguity. Compare with conventional charge dequantization and an explicit Majorana assumption.

| Gate | Frozen criterion |
|---|---|
| Representation and anomaly | All required charge, center and anomaly equations equal zero exactly for the family. Perturbative anomaly cancellation and the spin global-anomaly input must be stated separately. |
| Symmetry | Renormalizable terms allow an anomaly-free continuous $B-L$ modulo gauge. The two explicit breaking directions have rank two modulo gauge. The proposed dressed witness is local and nonzero. |
| Dimension | Prove the declared analytical lower and upper bounds on the least engineering dimension for $B-L$ violation. Verify the controls on both signs of $k$ and the $k=0$ endpoint. A bound is distinct from an exact minimum. |
| Inequivalence and reality | With QCD-like masses and a long-lived $udd$ hadron, predict $q_n/|q_e|$. A nonzero branch is rejected if the absolute value implied by the experimental benchmark is smaller than its exact minimum. Compare to conservative $|{\rm center}|+5$ times the reported uncertainty as a comparison scale. Do not interpret far-tail extrapolation as a statistical significance. |
| Reproduction | A standard-library production implementation and separate internal check agree exactly on the registered integer fixtures and charge ratios. Prove the infinite-family statements algebraically. |
| Decision | A requires a distinctive independent empirical restriction. B requires a useful conventional restriction with the selection goal satisfied at least conditionally. C applies if the decisive inequivalent witnesses are established and the selection goal fails. D applies if the candidate cannot be defined defensibly or neither restriction nor counterexample is justified. |

The primary conventional comparator is the same ordinary QFT and measured electric-charge assignments, with no further law-selection principle. The effect target is an electric monopole charge of the free neutron. For the selected assignment to be unique, a pre-data derivation must select $k$ without supplying neutron, neutrino, or proton neutrality or equality premises. The measured neutron charge is held outside the derivation. An independent demonstrated unbroken-phase $B-L$ violation at sufficiently low operator dimension would restrict the class, and this observation will be treated as a conditional control.

## Budget and stopping

One analytic proof, the five values $k=-2,-1,0,1,2$, exact rational arithmetic in two implementations and a byte-for-byte replay. No numerical fitting, simulation, external dependencies or parameter scan. Stop when the symbolic all-integer proof and independent criticism resolve the scope, or record an inconclusive test. No additional principle or next stage is preauthorized by a negative result.

The reviewer will separately challenge global symmetry, faithful group action, the nonzero contraction, the operator-grade lower bound, EFT/RG and UV assumptions, charge-measurement independence and literature priority. External RAI receives a standalone packet for review.
