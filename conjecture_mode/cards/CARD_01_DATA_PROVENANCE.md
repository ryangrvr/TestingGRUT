# CARD #1 — DATA PROVENANCE (preflight stage)

**Preflight:** 2026-10-04T12:04Z.
**Raw log:** `CARD_01_NETWORK_PREFLIGHT.log`.
**Method:** HEAD requests (`curl -sSI -L`), `git ls-remote`, and one blobless no-checkout clone of the ACT lensing **software** repository, used only to list file names.
**Products acquired: NONE.** No likelihood data file, chain, measurement vector or covariance was downloaded or opened. No SHA-256 product hashes exist, because no product was acquired.

**Order of commits:**
1. owner ruling `86ccbb2e`;
2. run-configuration freeze `c9899abc`;
3. then this preflight.

## Pinned software (installed from PyPI, which is reachable)

| Package | Version | Upstream tag / commit (from `git ls-remote`) |
|---|---|---|
| cobaya | 3.6.2 | `CobayaSampler/cobaya` tag v3.6.2 → `899f30a49f85de610dac321e91a1af50018e56aa` |
| camb | 2.0.4 | PyPI |
| act_dr6_lenslike | 1.2.1 | `ACTCollaboration/act_dr6_lenslike` tag v1.2.1 → `b386ddbb5821c1216c709f051c9289292f174d30` |
| Py-BOBYQA | 1.5.0 | PyPI |

## D3 components: access and provenance status

| Component (S6) | Distribution | Pinned ref | Access | O3 status |
|---|---|---|---|---|
| DESI DR2 BAO (`bao.desi_dr2.desi_bao_all`) | `CobayaSampler/bao_data` release v2.6 | tag v2.6 → `b7b8a36e9bccb063081f811f323cada21ab5fbdd` | GitHub archive tarball: **403**. Git transport (`ls-remote`): reachable | Not acquired. Traceability to the DESI DR2 release could not be cross-checked, because the official DESI hosts are blocked (below) |
| Planck PR4 NPIPE CamSpec TTTEEE | `CobayaSampler/planck_native_data` release v1, asset `CamSpec_NPIPE.zip` | tag v1 → `abcbf96ae4a88f966f572b4209ef45d35098c4b4` | **200** | Not acquired (preflight did not pass as a whole) |
| Planck low-ℓ TT (Commander, native) | same release, `planck_2018_lowT.zip` | same | **200** | Not acquired |
| Planck low-ℓ EE (SimAll, native) | same release, `planck_2018_lowE.zip` | same | **200** | Not acquired |
| **Planck + ACT DR6 lensing, likelihood v1.2** | data tarball `ACT_dr6_likelihood_v1.2.tgz`. The only documented sources are **NASA LAMBDA** (`lambda.gsfc.nasa.gov`, per the package code and README) and NERSC (`portal.nersc.gov`, chains). The software repository at v1.2.1 contains **no data** (file listing: code, tests, `get-act-data.sh`) | v1.2 | **BLOCKED**: `CONNECT tunnel failed, response 403` for both `lambda.gsfc.nasa.gov` and `portal.nersc.gov` | **FAILS ACCESS.** No collaboration-maintained distribution is reachable. A substitute is forbidden (O3; Charter §4B) |
| Pantheon+ / Union3 / DESY5 (extensions) | `CobayaSampler/sn_data` release v1.8 | tag v1.8 → `61d96434cafc2770928322c38e5a750e686368ae` | archive tarball **403**; git transport reachable | Not acquired |

## D1/D2 (official DESI w0wa chains)

| Host | Result |
|---|---|
| `data.desi.lbl.gov/public/` | **403** (blocked) |
| `data.desi.lbl.gov/public/papers/y3/` | **403** (blocked) |
| `www.desi.lbl.gov` | **403** (blocked) |

**D1/D2: UNAVAILABLE.** Under O3, no mirror or reconstructed chains were sought or used.

---

# RESUMPTION — owner enabled full network access (2026-10-04T12:08Z)

The ruling (`86ccbb2e`) is unchanged, and so is the run configuration (`c9899abc`, apart from the pre-evaluation corrections CR-1 to CR-4 in `CARD_01_RUN_CORRECTIONS.md`). The re-preflight is appended to `CARD_01_NETWORK_PREFLIGHT.log`.

## Acquisition routes

All acquired files are hashed in `CARD_01_PRODUCT_SHA256.txt`.

| Component | Route | Pin | O3 |
|---|---|---|---|
| DESI DR2 BAO (`bao.desi_dr2.desi_bao_all`) | git clone of `CobayaSampler/bao_data`. The GitHub archive and web routes still return 403 for this repository; git transport works | tag v2.6 = commit `b7b8a36e…`, verified MATCH | **PASS.** Pinned commit. All files hashed. Standard `bao.desi_dr2` implementation. **Traceability:** the official DESI DR2 README (`data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.md`) states that "the DESI DR2 BAO likelihoods used to produce these results are publicly available at github.com/CobayaSampler/bao_data/tree/master/desi_bao_dr2". The `desi_bao_dr2` git tree at v2.6 (`0c850733…`) is **identical** to the tree at `master` (`bb0c1c90…`) |
| Planck low-ℓ TT Commander / EE SimAll (clik) | official Planck Legacy Archive `COM_Likelihood_Data-baseline_R3.00`, downloaded by `cobaya-install`, read by clipy-like 0.15 (`benabed/clipy` tag `clipy_0.15` = `ad1aff3b…`, git clone; the clipy PyPI package is the same version) | PR3 baseline | **PASS** (collaboration distribution) |
| Planck PR4 NPIPE CamSpec TTTEEE | `CobayaSampler/planck_native_data` release v1 asset `CamSpec_NPIPE.zip` (HTTP 200), via `cobaya-install` | tag v1 = `abcbf96a…` | **PASS** (the standard public distribution, used by DESI through `planck_NPIPE_highl_CamSpec.TTTEEE`) |
| Planck + ACT DR6 lensing v1.2 | NASA LAMBDA `ACT_dr6_likelihood_v1.2.tgz` (HTTP 200), via `cobaya-install` / `act_dr6_lenslike` 1.2.1 | v1.2 | **PASS** (collaboration distribution; variant `actplanck_baseline`, as in the DESI input) |
| Pantheon+ / Union3 / DESY5 | git clone of `CobayaSampler/sn_data` (archive route 403) | tag v1.8 = `61d96434…`, verified MATCH | **PASS** (pinned commit, hashed; standard `sn.*` implementations) |
| `planck_supp_data_and_covmats` (an install dependency only; not a likelihood input) | git clone | tag v2.1 = `621f4c5d…` | — |

## D1/D2

The official DESI chains are now reachable at `data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya/base_w_wa/`. Only the DESI-released chains will be used.
