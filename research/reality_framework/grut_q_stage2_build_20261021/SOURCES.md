# Sources, credits and evidence scope

21 October 2026. This contribution is based on the published Claude
checkpoint 8ec108c and the preceding Codex response contribution a5b9f44.
No newer Claude stage-2 result was published when the builder handoff
was checked. No background agent completion is assumed.

## Primary papers used

| Source | Read scope in this continuation | Use |
|---|---|---|
| Klinkhamer--Volovik, [0806.2805](https://arxiv.org/abs/0806.2805), PRD 78, 063528 (2008) | Section II was read in the preceding response pass; the current work uses that pinned action and the stage-1 map | Credits q theory and its fixed-mu scalar-tensor representation |
| Klinkhamer--Volovik, [1601.00601](https://arxiv.org/abs/1601.00601), MPLA 31, 1650160 (2016) | Sections 5--7 source equations, their qualifications and late attractor; related stage-1 notes cross-checked | Credits curvature-produced receiver energy and the dissipative q example; our reduced-Planck-unit equations and finite-N calculation are explicit |
| Klinkhamer--Savelainen--Volovik, [1601.04676](https://arxiv.org/abs/1601.04676), JETP 125, 268 (2017) | Abstract and the stated reversible/dissipative distinction, with stage-1 source notes | Establishes why adding dissipation alone is not a blanket attractor claim; its full alternate system is not rerun |
| Lee et al., [2002.11761](https://arxiv.org/abs/2002.11761), PRL 124, 101101 (2020) | Abstract, Yukawa definition, Figure 6 viewed as a page image, and conclusion | Historical inverse-square-law constraint; approximate alpha=1/3 crossing is a visual bracket only |
| Froggatt--Nielsen, [hep-ph/9511371](https://arxiv.org/abs/hep-ph/9511371) | Opening assumptions and high-field lambda/beta criticality argument | Credits the multiple-point condition; its old mass estimates are not adopted as current predictions |
| Buttazzo et al., [1307.3536](https://arxiv.org/abs/1307.3536) | Sections 3.1--3.5 historical coupling/matching statements; relevant leading RGE conventions in appendix B | Supplies gauge inputs and conventions; this work runs one-loop equations only, not their precision code |

The source papers use differing curvature, Planck-mass, Higgs-quartic
and hypercharge conventions. Here R=12H² in de Sitter, M is reduced
Planck mass, V_H=lambda_H h^4/4 and gY is non-GUT-normalized. The
distinctions are needed for the numerical coefficients quoted.

All source observations are already exposed. This is not a current-global-
constraint compilation, a new laboratory confidence limit or a held-out
cosmological test. No unpublished raw experimental likelihood was used.

## Work completed here

The code implements the conserved scalar transition, a second formulation
in the Einstein frame, the full leading background source correction,
conditional multiple-point running, a radiation-receiver model for mu
formation, and an explicit positively biased metastable-potential version.
The scalar potential and phase preparation in QV1-M are **supplied**
choices. They are conventional EFT constructions, not a novelty claim.

The nondissipative q action does not derive the transport source of the
dissipative comparison. Its effective changing mu must not be identified
silently with the original bare action's conserved canonical flux. A
microscopic action/current and its quantum noise remain unconstructed.

## Verification meaning

One Codex author carried out this continuation. No subagents, separate
reviewers or Claude background results are claimed. The separate file
verify.py uses different variables and numerical methods without importing
build.py, but it is still a same-author check informed by the construction.
This is weaker than independent scientific reproduction or outside review.

The numerical scales are compressed hierarchy controls. For example the
dimensionless baseline m_q=20--80 is not a quantum-gravity EFT benchmark
with those masses above its Planck scale. The physical low-curvature
illustration uses separate stated inputs and analytic decoupling; it is
not a direct numerical replay of electroweak-to-Planck time scales.

Parameter values, state choices, numerical tolerances and reference
physical constants are in INPUTS.json. The code/result hashes establish
source identity, not empirical truth or an all-times error certificate.
There is no preregistration claim. No dependencies were installed.
