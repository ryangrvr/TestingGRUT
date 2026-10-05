# CONJECTURE MODE — IMPORT MANIFEST 01

**Source:** `ryangrvr/TestingGRUT`, branch `grut-conjecture-mode-1`.
- Files in §1–§2 are pinned to the source commit `a8dbf3f45b4502b65913dcbb709c171580ef5e34` ("CARD-01 v1R2 TERMINAL: COMPUTATIONALLY UNRESOLVED").
- `CONJECTURE_MODE_FINAL_SYNTHESIS_01.md` and this manifest are added in the closing commit.

**Hashes.** `sha256` is computed over the exact git blob content at the source commit (`git show <sha>:<path> | sha256sum`).

**Rule for importers.** Import byte-exact. **Never edit an imported file.** Put status in the destination capsule README and the import manifest.

**Roles:**
- ACTIVE GOVERNANCE REFERENCE: rules a future campaign may cite.
- ARCHIVED SCIENTIFIC RECORD: frozen record; status must not change.
- REPRODUCIBILITY SUPPORT: code, configurations, logs and inputs.

**Counts:** 54 files at the source commit = 3 `program_governance/` + 51 `conjecture_mode/`. Adding the synthesis file gives **55 files to import**.

## 1. Closing document

| Path | Role | sha256 |
|---|---|---|
| `CONJECTURE_MODE_FINAL_SYNTHESIS_01.md` | ARCHIVED SCIENTIFIC RECORD (closing synthesis) | `fe1d1304cced67c3d380c0510c78d40e36ee39ada39e0301762a0dddb2743374` |

## 2. Files to import (pinned to the source commit)

| Path | Role | sha256 |
|---|---|---|
| `conjecture_mode/CARD_LEDGER.md` | ARCHIVED SCIENTIFIC RECORD | `536df2aa74ca63446a2148ed05fde49fb9c3b9fb50efaa1dd8472bc81d1cbcb8` |
| `conjecture_mode/CONJECTURE_GENERAL_REGISTER.md` | ARCHIVED SCIENTIFIC RECORD | `0a562e8ee69e38d92b918c6f8114a79ee868d949f76277a6026861950e15cfe8` |
| `conjecture_mode/EQUIVALENCE_REGISTRY.md` | ARCHIVED SCIENTIFIC RECORD | `b73e4d03716ad7661cfc4b9eb79a0b55cf972a278262fbb3f169ee4da4148a0a` |
| `conjecture_mode/cards/CARD_01_D1_OFFICIAL_CHAINS.sha256sum` | ARCHIVED SCIENTIFIC RECORD | `48bf1e172f4ea3e2d754bdc9b139e688d4709d83727c3a23fdb6d61e880b675a` |
| `conjecture_mode/cards/CARD_01_DATA_PROVENANCE.md` | ARCHIVED SCIENTIFIC RECORD | `2ea29cd53657cf43407276ee92f98dece7e0dc96af78abc883293d16c358fe1c` |
| `conjecture_mode/cards/CARD_01_EXECUTION_REPAIRS.md` | ARCHIVED SCIENTIFIC RECORD | `b93ae8d4cc575bd783329e9e01efe7fbf0b0e9db992062033e92381958b885a4` |
| `conjecture_mode/cards/CARD_01_HORIZON_RELAXOR_SPEC.md` | ARCHIVED SCIENTIFIC RECORD | `bc6b3e8bdf4aa7fbe9c2068e4c370651b289daae4a4e25ce54c0a8a2c814c038` |
| `conjecture_mode/cards/CARD_01_NETWORK_PREFLIGHT.log` | ARCHIVED SCIENTIFIC RECORD | `a993900b5090a3b06a76e1701b23312d0fd0958ff598a14fc0ab6a0abda261c7` |
| `conjecture_mode/cards/CARD_01_OWNER_CONVERGENCE_RULING.md` | ARCHIVED SCIENTIFIC RECORD | `bca421e87b00b98c3747bef829b87e3fac4b9a34fe1d2be83054819684be0f7c` |
| `conjecture_mode/cards/CARD_01_OWNER_THRESHOLD_RULING.md` | ARCHIVED SCIENTIFIC RECORD | `f3a05b3cb3d628980dac9f1f053dd99931e680a6a244a9c12c9f7b3232b14a64` |
| `conjecture_mode/cards/CARD_01_PRODUCT_SHA256.txt` | ARCHIVED SCIENTIFIC RECORD | `a0d3fa845f6a958ad81ce797e94dfc60b7ff02c23e19b80a3954be369ac079f7` |
| `conjecture_mode/cards/CARD_01_RESULT.md` | ARCHIVED SCIENTIFIC RECORD | `12150a7c39f18be745b586452076e71c26f9ed104fb42bb2de9a4b9ebe2971af` |
| `conjecture_mode/cards/CARD_01_RUN_CORRECTIONS.md` | ARCHIVED SCIENTIFIC RECORD | `bcf14aa40f581d96e6aaae310736868d2f011f614bd205efef1ff0b3be15c64d` |
| `conjecture_mode/cards/CARD_01_V1R2_CHARTER.md` | ARCHIVED SCIENTIFIC RECORD | `20598f741d30138cb146147821385b0c8aa3d64a09656b3b07660044aee08b2c` |
| `conjecture_mode/cards/CARD_01_V1R2_RESULT.md` | ARCHIVED SCIENTIFIC RECORD | `d51e43c8340507d55dfa46482ac749348bd0fc2c69ee4ddd1706780be18aac0e` |
| `conjecture_mode/cards/CARD_01_V1R_CHARTER.md` | ARCHIVED SCIENTIFIC RECORD | `0c5cc7a1bff21377f4d305c95bec313f426fb3b430297d7a4429b232067acee5` |
| `conjecture_mode/cards/CARD_01_V1R_OFFICIAL_CHAINS.sha256sum` | ARCHIVED SCIENTIFIC RECORD | `86015fea7ee56eb173ea1b7cb5aec031c99777da073f8cdfad0123dca4f63802` |
| `conjecture_mode/cards/CARD_01_V1R_RESULT.md` | ARCHIVED SCIENTIFIC RECORD | `d10aff30015f4c32c99d116e593a05ae783e52629583fa2c1cc83ecef3a55281` |
| `conjecture_mode/cards/code/card01_cobaya_plumbing_test.log` | REPRODUCIBILITY SUPPORT | `ba85fb2b79aae94fce1888e98e19e77c89bc5692371f7bd086010878c47ba73d` |
| `conjecture_mode/cards/code/card01_cobaya_plumbing_test.py` | REPRODUCIBILITY SUPPORT | `8f182af223f856ba10d6249a8b93752471521d187a9971bd75eab3e19146be5a` |
| `conjecture_mode/cards/code/card01_cobaya_theory.py` | REPRODUCIBILITY SUPPORT | `d69ecddb75e1ff1c9cedbe97ead1b624262ec0aa9e3eb4a8a411ce7aa2c2880d` |
| `conjecture_mode/cards/code/card01_cpl_projection.log` | REPRODUCIBILITY SUPPORT | `0d8460d766d3e9e65d00dd5664124ff3776ff2a466974ddbaaf9982d11ea3b5d` |
| `conjecture_mode/cards/code/card01_cpl_projection.py` | REPRODUCIBILITY SUPPORT | `633134d038df052ee874e5fdb8e5383cc983d0304db1d6c1e181710076542f9e` |
| `conjecture_mode/cards/code/card01_d1d2_diagnostic.json` | ARCHIVED SCIENTIFIC RECORD (result summary; located in code/) | `7e5493d6037973dbae097f108bc19c09add42d15d06c8b178216d7bbf5dee7eb` |
| `conjecture_mode/cards/code/card01_d1d2_diagnostic.log` | REPRODUCIBILITY SUPPORT | `c4bdbc003a2d3b21cc7204b693e10cd1229a5767f5a096242760837005fafd3a` |
| `conjecture_mode/cards/code/card01_d1d2_diagnostic.py` | REPRODUCIBILITY SUPPORT | `1fe7a6b3e64e31519e559b32a90b4596e5ec8c8f1deccb3ee067809198a59e75` |
| `conjecture_mode/cards/code/card01_d3_analyze.py` | REPRODUCIBILITY SUPPORT | `0f489630dccf51b530ed9ab1e09b8908490457cedb6aa080742fbb8214bc9460` |
| `conjecture_mode/cards/code/card01_d3_driver.py` | REPRODUCIBILITY SUPPORT | `29545e4ae8997752956fca3ad879f2e28e77b92cf0800ae40591462b85cbf158` |
| `conjecture_mode/cards/code/card01_d3_feasibility.log` | REPRODUCIBILITY SUPPORT | `03fdcc420090222d3e48b6c92b66c6a537d50930cabbb3171432fea779760e17` |
| `conjecture_mode/cards/code/card01_d3_feasibility.py` | REPRODUCIBILITY SUPPORT | `c0e345e66a4bd54e8473a7a6cb9c401234205e41bfc5aa1f5b414e7bf5344f57` |
| `conjecture_mode/cards/code/card01_d3_primary_summary.json` | ARCHIVED SCIENTIFIC RECORD (result summary; located in code/) | `9f3c2163403bfc1d388ce297fe6aa7084597a3a801873543a94540698d5c2928` |
| `conjecture_mode/cards/code/card01_d3_run.py` | REPRODUCIBILITY SUPPORT | `78c50644acd855bee1e6d6c8a9e7428ba19a95a6c64b68a01153237138ccf156` |
| `conjecture_mode/cards/code/card01_run_config.py` | REPRODUCIBILITY SUPPORT | `62f1fb4f974b85072cee29d04f98e5125b97b04c95aaf6df0c117f56e1460bda` |
| `conjecture_mode/cards/code/card01_v1r2_analyze.py` | REPRODUCIBILITY SUPPORT | `5e247efcbd0d9e5c13b90328146d63aa42b9908d7d8dc09d7b8ff6971147009c` |
| `conjecture_mode/cards/code/card01_v1r2_driver.log` | REPRODUCIBILITY SUPPORT | `7da78fb0f57b190a3573ac41d5df6bf1731656f303f1d94ec995815042402c76` |
| `conjecture_mode/cards/code/card01_v1r2_driver.py` | REPRODUCIBILITY SUPPORT | `1dda6d80724c7184659d21df300db4786b47c8cc2b2fc468367982a0e5cd22ff` |
| `conjecture_mode/cards/code/card01_v1r2_lib.py` | REPRODUCIBILITY SUPPORT | `3520821b616bf2767ac8cf198877fb6390a56ad9f80edbf68ef5e2a1e916d4d0` |
| `conjecture_mode/cards/code/card01_v1r2_run.py` | REPRODUCIBILITY SUPPORT | `5d87e419d49ad8a55518f53edc03e6ba673e07b5f990936799ac8311c4bcf4e1` |
| `conjecture_mode/cards/code/card01_v1r2_starts.py` | REPRODUCIBILITY SUPPORT | `ff1f56c2d288f8544a7cecbc03065663b3d383e8450081ffceee96be788f0dab` |
| `conjecture_mode/cards/code/card01_v1r2_summary.json` | ARCHIVED SCIENTIFIC RECORD (result summary; located in code/) | `f6954e7f38500344ad405d62a3f5621324c8ca73406440f80bba3c4abbaef4ca` |
| `conjecture_mode/cards/code/card01_v1r_analyze.py` | REPRODUCIBILITY SUPPORT | `bf8dc3c632f2e0c117e630df92b9762dc3bd99ec326bbf749f9ce701099b87b5` |
| `conjecture_mode/cards/code/card01_v1r_config.py` | REPRODUCIBILITY SUPPORT | `1d29c75e0eedf62f3e4b3ee71244707ac119db08860fa6dca77f6253ed0bc4d5` |
| `conjecture_mode/cards/code/card01_v1r_driver.py` | REPRODUCIBILITY SUPPORT | `fd2577746941b7f410840ca48b2ab94362f2164d50c5ffecbdd6c2208296e643` |
| `conjecture_mode/cards/code/card01_v1r_inputs.py` | REPRODUCIBILITY SUPPORT | `6562c0677bb4ce5095e3b36c735445cd12c49a6eaafd1ee9d9febe8ea4ca98fe` |
| `conjecture_mode/cards/code/card01_v1r_primary_summary.json` | ARCHIVED SCIENTIFIC RECORD (result summary; located in code/) | `3ed4330ac725a3803bf2d28cb66b892248b5bf0d1c888e7ee37968b22d5ff64d` |
| `conjecture_mode/cards/code/card01_v1r_refine.py` | REPRODUCIBILITY SUPPORT | `560cbd53c62459c3fd44ca34674b9f6f0360946ea8dac90647098fce6754c210` |
| `conjecture_mode/cards/code/v1r2_inputs/stageA_starts.json` | REPRODUCIBILITY SUPPORT | `083bc1169bbf77f0d6161f7fca407d68695c9c0a6e674b7185986945c4f169a7` |
| `conjecture_mode/cards/code/v1r_inputs/v1r_covmat_lcdm.txt` | REPRODUCIBILITY SUPPORT | `9a968281b6726dcef2c76bfe7c7a39dc1df2fd0b91371082c8be2590aa26be0a` |
| `conjecture_mode/cards/code/v1r_inputs/v1r_covmat_w0wa.txt` | REPRODUCIBILITY SUPPORT | `a49ac04db7fbabaac697dcea51c357c3876148090b6a9d8bdc8bd40ab3fb2667` |
| `conjecture_mode/cards/code/v1r_inputs/v1r_cpl_startB.json` | REPRODUCIBILITY SUPPORT | `6a232928de2defce1f9addd1a65112a1b16d5456111645df9f16f25a84110898` |
| `conjecture_mode/cards/code/v1r_inputs/v1r_nuisance_param_defaults.json` | REPRODUCIBILITY SUPPORT | `5dbf20f168407ea836d28c9fa6b894ccb8f293325f1bbdc5a44a59153ddc1ffe` |
| `program_governance/CONJECTURE_MODE_CHARTER_01.md` | ACTIVE GOVERNANCE REFERENCE | `a1fcdd3a14a02eff0a3c0282cd4f681f71e98c6f8c1af107cb0a1b898f796d8a` |
| `program_governance/PROGRAM_CAMPAIGN_GATE_01.md` | ACTIVE GOVERNANCE REFERENCE | `067ea7d8935f51f0bf6dd6e2c82ed3ed73d114ac5c6fd4a00af6cb660ef7fc9f` |
| `program_governance/PROGRAM_GOVERNANCE_OWNER_RULING_01.md` | ACTIVE GOVERNANCE REFERENCE | `9f5c893c0311dfcba9978887f2dd0dd1fd96d9c539789dd1947d0cb1c4bcba13` |

## 3. Omitted from this import (present on the branch, belonging to other campaigns)

`grut-conjecture-mode-1` descends from `grut-program-governance-1-frozen`, which descends from `grut-selector-screen-1-frozen` and the earlier frozen campaign branches. These files belong to those campaigns. They should be imported, if at all, as their own capsules from their own frozen branches, not as part of Conjecture Mode.

| Top-level path | Files | Reason |
|---|---|---|
| `(repository root files)` | 4 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `bridge` | 33 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `charter` | 1 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `gravity_scout_1` | 19 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `ledgers` | 8 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `probes` | 43 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `provenance` | 3 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `qft_scout_1` | 22 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `reconnaissance` | 3 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `residual_synthesis` | 7 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `results` | 17 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `review` | 31 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `selector_screen_1` | 15 | inherited from an earlier campaign branch; not Conjecture Mode material |
| `zoomouts` | 8 | inherited from an earlier campaign branch; not Conjecture Mode material |

## 4. NOT IN THE REPOSITORY (owner decision needed)

The per-fit raw outputs and run logs were written **only to the executing container's ephemeral scratch space**. They were **never committed**. If the container is reclaimed, they are lost.

The committed record keeps the analyzer summary JSONs (§2) and every driver, configuration and result document. The raw outputs behind those summaries are not otherwise reproducible without rerunning the fits.

| Scratch item | Content | Size |
|---|---|---|
| `d3_out/` | v1 primary: 66 per-fit JSONs | 268 KB |
| `d3_out_archive_misseeded_free/` | v1: 2 archived, excluded mis-seeded free-ε fits (owner: do not delete) | 12 KB |
| `v1r_out/` | v1R: 29 per-target JSONs | 120 KB |
| `v1r2_out/` | v1R2: 10 per-start JSONs plus start files | 88 KB |
| `d3_logs/`, `v1r_logs/`, `v1r2_logs/` | per-fit logs | ≈ 1.2 MB |
| Driver logs (`d3_driver*.log`, `v1r_driver*.log`) | phase logs, including the OOM / restart attempts | small |

The large likelihood data products (≈ 4 GB) and the official DESI chains are **deliberately not committed**. They are re-acquirable from pinned sources, and their hashes are committed in `CARD_01_PRODUCT_SHA256.txt`, `CARD_01_D1_OFFICIAL_CHAINS.sha256sum` and `CARD_01_V1R_OFFICIAL_CHAINS.sha256sum`.


## 5. Notes for the destination importer

- **Registries.** `conjecture_mode/EQUIVALENCE_REGISTRY.md` (EQ-01) and `conjecture_mode/CONJECTURE_GENERAL_REGISTER.md` (no registered entries) are archived as-is. Any program-wide registry should cross-reference them, not copy or edit them.
- **Governance dependency.** `PROGRAM_GOVERNANCE_OWNER_RULING_01.md` is what freezes the other two governance files. All three must travel together.
- **Path dependencies.** Code under `conjecture_mode/cards/code/` expects its sibling layout (`v1r_inputs/`, `v1r2_inputs/`) and environment variables for data paths. It is reproducibility support, not an executable guarantee: re-execution requires re-acquiring the pinned data products.
