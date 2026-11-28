# Commit-pinned reproduction index

## Immutable source stack

Start from a clone of TestingGRUT with the listed commits and complete packet inputs. Fetch the stacked branches if the clone follows only the default branch. Published drafts link directly to the protocol, retained data, scripts and review at the revisions below. Use the result revision for a replay with archived outputs and the earlier pin when checking what was frozen before a run. [Repository](https://github.com/ryangrvr/TestingGRUT)

| Packet | Draft | Protocol or source pin | Produced result pin |
| --- | --- | --- | --- |
| S1 integrated synthesis | [#29](https://github.com/ryangrvr/TestingGRUT/pull/29) | [d3dd471a0ad69a95cd6f0dc1eba23fc0695a39cd](https://github.com/ryangrvr/TestingGRUT/commit/d3dd471a0ad69a95cd6f0dc1eba23fc0695a39cd) | [a2b6d3a2a806e52b088ad3fa10f8ce789d6659a9](https://github.com/ryangrvr/TestingGRUT/commit/a2b6d3a2a806e52b088ad3fa10f8ce789d6659a9) |
| C1 fixed charge and one loop | [#30](https://github.com/ryangrvr/TestingGRUT/pull/30) | [a3ff9f40a1588501a24afc87de996668ef237374](https://github.com/ryangrvr/TestingGRUT/commit/a3ff9f40a1588501a24afc87de996668ef237374) | [65f13b0b30b162e3cf86e5b8fe41931545b8f94d](https://github.com/ryangrvr/TestingGRUT/commit/65f13b0b30b162e3cf86e5b8fe41931545b8f94d) |
| C2 continuum stability | [#31](https://github.com/ryangrvr/TestingGRUT/pull/31) | [fc17fb4a08e5bb81b0260010c33439ca74767f11](https://github.com/ryangrvr/TestingGRUT/commit/fc17fb4a08e5bb81b0260010c33439ca74767f11), with [amendment da9b5dcf12bf677854d2c4560bdf37a0461dbac0](https://github.com/ryangrvr/TestingGRUT/commit/da9b5dcf12bf677854d2c4560bdf37a0461dbac0) | [414e5c2a0aa601851e27808a3deaea1db15af819](https://github.com/ryangrvr/TestingGRUT/commit/414e5c2a0aa601851e27808a3deaea1db15af819) |
| C3 reporter feasibility | [#32](https://github.com/ryangrvr/TestingGRUT/pull/32) | [35a70c726abf8b5ea684c428d1f50a248647fc15](https://github.com/ryangrvr/TestingGRUT/commit/35a70c726abf8b5ea684c428d1f50a248647fc15) | [0d1ff1e428fe49ffbdb55e5083da5f9e1df0d13c](https://github.com/ryangrvr/TestingGRUT/commit/0d1ff1e428fe49ffbdb55e5083da5f9e1df0d13c) |
| C4 gravity and closure | [#33](https://github.com/ryangrvr/TestingGRUT/pull/33) | [54d4de2872fb1d9b372de79d25e9849186bc6ecc](https://github.com/ryangrvr/TestingGRUT/commit/54d4de2872fb1d9b372de79d25e9849186bc6ecc) | [6c96cd26b223097fb93f15408da1590b47acd63d](https://github.com/ryangrvr/TestingGRUT/commit/6c96cd26b223097fb93f15408da1590b47acd63d) |

C4's tree SHA is 76e804d1e22ec20c93dd430ea03e5df9ea2ede15. The recursive NUL-delimited tree inventory SHA-256 is 762465abda4f9ecfbdc243f269fbc40a74914682d49c7dc58167f19795fee74b. The corresponding S1 digest is 21eef8faec1355321591afdeb077f0aad5ac66f2c8de744aea3808e0cf2284f1. Reproduce the C4 tree using:

~~~sh
c4=6c96cd26b223097fb93f15408da1590b47acd63d
git rev-parse "$c4^{tree}"
git ls-tree -rz "$c4" | sha256sum
git ls-tree -rz --name-only "$c4" | python3 -c 'import sys; p=sys.stdin.buffer.read().split(b"\0")[:-1]; print(len(p), sum(x.endswith(b".md") for x in p))'
git diff --name-status a2b6d3a2a806e52b088ad3fa10f8ce789d6659a9 "$c4"
~~~

The C4 counts are 8,122 and 1,722. The S2 diff against C4 should add the 12 requested documents and append one coordination entry. The C2 failed-run result is pinned at [b7a913bcabd304ab543c06f2fb59f00cc31f0681](https://github.com/ryangrvr/TestingGRUT/commit/b7a913bcabd304ab543c06f2fb59f00cc31f0681). [S1 corpus and replay](../grut_s1_integrated_synthesis_20261127/REPRODUCTION_INDEX.md)

## Scalar matter replays

Use a detached worktree for each desired return point, and create fresh result paths for full replays. A complete C4 checkout also contains the unchanged C1–C3 replay sources. A sparse checkout needs the four C-series packets as well as the original Candidate 1 inputs for the C1 full replay. Record the actual numeric environment. The original C1 host used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and C++17 with OpenMP. C2 required 63 long-double fraction bits, round-to-nearest with outward nextafter and mpmath 1.3.0 with the documented interval precision. C4 inherits the relevant C2 interval conditions.

| Packet | Run from its tracked packet directory | Verification target |
| --- | --- | --- |
| C1 archive only | python3 verify_packet.py | Frozen protocol/source hashes, retained raw log integrity, result flags, signs, absolute allowances and workload. [C1 instructions](../grut_c1_quantum_channels_20261127/REPRODUCTION.md) |
| C1 full replay | Use the profile regeneration and OpenMP production commands in the frozen packet and a fresh output location. | Paired and separate one-carrier intervals, equal-fission interval, signed bounds and raw outcomes. The archive integrity check is exact. New wall times and parallel log ordering can vary. [C1 instructions](../grut_c1_quantum_channels_20261127/REPRODUCTION.md) · [Protocol](../grut_c1_quantum_channels_20261127/PROTOCOL.json) |
| C2 | One-thread, 180 s bounded certify.py command. | Three exact-charge results, initial pivot failure and corrected proof, archived RESULTS.json and EXECUTION.log. [C2 instructions](../grut_c2_continuum_stability_20261127/REPRODUCTION.md) |
| C3 | One-thread, 90 s bounded scatter.py command. | Original reference phase, parent phase, retained angular sum and empirical refinements, RESULTS.json and EXECUTION.log. [C3 instructions](../grut_c3_empirical_feasibility_20261127/REPRODUCTION.md) |
| C4 | One-thread, 45 s bounded verify_gate.py command with C1–C3 files present. | Pinned cross-input hashes, compactness and source scale, nonvacuum tail and isolated material tide bounds, RESULTS.json and EXECUTION.log. [C4 instructions](../grut_c4_gravitational_accessibility_20261127/REPRODUCTION.md) |

The archived C1 profile SHA-256 is 8ecb6242f74a6e0a2e180d617ac244db1c8abaebdb3f8a7a0f0bcfd025ac9f0a. C2 and C4 store this input hash. The C4 result also stores the SHA-256 of all three earlier results. When byte comparison fails across arithmetic environments, retain both outputs, identify the discrepancy and separately check the registered inequalities. [C2 hash and bounds](../grut_c2_continuum_stability_20261127/RESULTS.json) · [C4 hashes and bounds](../grut_c4_gravitational_accessibility_20261127/RESULTS.json)

## Developmental and foundational results

Use the immutable S1 result pin above for D1–D7, I1–I3 and F1–F5 source versions accumulated at the end of S1. The earlier S1 index gives the per-packet frozen SHAs, PRs, proof files and program invocations. Its full source boundary includes separate legacy decisions and archive replay scripts. For the consolidated review start with these original instructions and use the frozen thresholds unchanged. [S1 per-campaign pins and commands at a2b6d3a](https://github.com/ryangrvr/TestingGRUT/blob/a2b6d3a2a806e52b088ad3fa10f8ce789d6659a9/research/reality_framework/grut_s1_integrated_synthesis_20261127/REPRODUCTION_INDEX.md)

| Scope | Primary replay entries |
| --- | --- |
| Developmental models D1–D7 | [D1/I1 evidence](../grut_i1_integrated_20261125/EVIDENCE_MAP.md) · [D2 with corrected control](../grut_d2_development_20261123/README.md) · [D3](../grut_d3_record_repair_20261124/README.md) · [D4](../grut_d4_internal_recognition_20261124/README.md) · [D5](../grut_d5_prospective_20261124/README.md) · [D6](../grut_d6_adaptive_20261125/README.md) · [D7 with both verify commands](../grut_d7_functional_recovery_20261125/README.md) |
| Integrated models I1–I3 | [I1](../grut_i1_integrated_20261125/REPRODUCTION.md) · [I2](../grut_i2_persistence_20261125/REPRODUCTION.md) · [I3 and 2,624-value independent comparison](../grut_i3_predictive_memory_20261126/REPRODUCTION.md) |
| Obstruction and reconstruction F1–F5 | [F1](../grut_f1_operational_reconstruction_20261126/REPRODUCTION.md) · [F2](../grut_f2_complementary_interactions_20261126/REPRODUCTION.md) · [F3](../grut_f3_coupled_consistency_20261126/REPRODUCTION.md) · [F4 and correction](../grut_f4_charge_selection_20261126/REPRODUCTION.md) · [F5](../grut_f5_closure_consistency_20261127/REPRODUCTION.md) |

## Older and parallel archives

The F5 archive index includes the historical hierarchy, Level-0 and successor closures, T3 absorptive calculation and amendment, CP-12, the original Candidate 1, separate matter and vacuum protocols, and older cosmology. The inherited pin table identifies exact parent and side-branch commits. The separate QED validation line includes finite massless passes and the live-source failure. Those data are separate from the counted F5 corpus. [Legacy sources at the S1 pin](https://github.com/ryangrvr/TestingGRUT/blob/a2b6d3a2a806e52b088ad3fa10f8ce789d6659a9/research/reality_framework/grut_s1_integrated_synthesis_20261127/REPRODUCTION_INDEX.md#older-and-parallel-archives) · [QED branch #14](https://github.com/ryangrvr/TestingGRUT/pull/14)

| QED record | Exact result pin |
| --- | --- |
| B1-01 and interacting parent | [7d57fc7c8a36c8f349165afc946d9323f174e630](https://github.com/ryangrvr/TestingGRUT/commit/7d57fc7c8a36c8f349165afc946d9323f174e630) |
| B2 N20 allocation miss | [1fa9f643e01f2b1a012e82cb8061592cba2a40ef](https://github.com/ryangrvr/TestingGRUT/commit/1fa9f643e01f2b1a012e82cb8061592cba2a40ef) |
| N22 finite massless pass | [748c74ff65206aff27ac550e404b02fab0989f9a](https://github.com/ryangrvr/TestingGRUT/commit/748c74ff65206aff27ac550e404b02fab0989f9a) |
| Source D1-01 failed registered moments and finest step | [7ba04fd35e47870cbdb41268b2d7ffd33f567f51](https://github.com/ryangrvr/TestingGRUT/commit/7ba04fd35e47870cbdb41268b2d7ffd33f567f51) |

## Verification status

S2 inspected the pinned source stack, checked the archived tree metadata, four protocol digests and all C4 cross-input digests, and ran the existing C1 archive-only verifier. The verifier confirmed 1,932 Jost points, 34 radial batches and the separate-allowance emission interval $[-0.253448528261,-0.146689876395]\,g_Xv$. The new work is evidence reconciliation, derivation exposition and literature comparison. Scientific production outputs and protocols are inherited from the named campaigns. Archived internal replays, internal critical checks and external replication retain their own status. [External RAI packet](EXTERNAL_RAI_REVIEW_PACKET.md)
