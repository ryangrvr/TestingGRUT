# CONJECTURE MODE — IMPORT MANIFEST 01 — ADDENDUM A (raw fit outputs and logs)

**Scope.** This addendum is additive to `CONJECTURE_MODE_IMPORT_MANIFEST_01.md` (frozen at `c7fdefd`). That manifest is **not edited**.

Manifest 01 §4 recorded these files as not yet in the repository. They are now committed, byte-exact, under `conjecture_mode/cards/raw/`. This preservation was approved in the owner's TestingGRUT final-closeout instruction.

**Source.** The executing container's scratch outputs for Card #1 v1, v1R and v1R2. The copy was verified byte-exact: 232/232 sha256 match between source and copy before commit.

**Hashes.** `sha256` is taken over the git blob content, i.e. the exact file bytes.

**Rule.** Import byte-exact. Never edit an imported file. These raw files change **no** frozen scientific status.

## 1. Layout

| Directory | Content |
|---|---|
| `raw/v1/fits/` | v1 PRIMARY: 66 per-fit JSONs (20 ε × 3 starts, CPL × 3, free-ε × 3) — source of `code/card01_d3_primary_summary.json` |
| `raw/v1/ARCHIVED_MISSEEDED_FREE_EXCLUDED/` | the 2 free-ε fits seeded from ε = 0 before the grid was complete (ER-B). **Archived and excluded; used for no statistic** |
| `raw/v1/logs/` | 69 per-fit logs, including the smoke test and the ER-B fits |
| `raw/v1/driver_logs/` | 4 driver logs: OOM attempt, two restart-interrupted attempts, final |
| `raw/v1r/fits/` | v1R: 29 per-target JSONs (22 arm A + 7 arm B) — source of `code/card01_v1r_primary_summary.json` |
| `raw/v1r/logs/`, `raw/v1r/driver_logs/` | 29 per-target logs; 2 driver logs (OOM attempt, final) |
| `raw/v1r2/fits/` | v1R2: 10 per-start JSONs — with `raw/v1r/fits/` and the v1R2 driver log, source of `code/card01_v1r2_summary.json` |
| `raw/v1r2/fits/starts/` | the 10 start points used by those runs |
| `raw/v1r2/logs/`, `raw/v1r2/driver_logs/` | 10 per-start logs; 1 driver log |

## 2. IMPORT/REPRODUCIBILITY CHECK — DOES NOT MODIFY FROZEN SCIENTIFIC STATUS

**Method.** The committed analyzers were re-run, unchanged, against the committed raw copies, writing to a scratch directory. Their output was then compared with the committed summary JSONs.

| Committed summary | Regenerated from | Result |
|---|---|---|
| `code/card01_d3_primary_summary.json` | `card01_d3_analyze.py raw/v1/fits` | **byte-identical** (sha256 `9f3c2163403bfc1d…`) |
| `code/card01_v1r_primary_summary.json` | `card01_v1r_analyze.py raw/v1r/fits` | **byte-identical** (sha256 `3ed4330ac725a380…`) |
| `code/card01_v1r2_summary.json` | `card01_v1r2_analyze.py raw/v1r/fits raw/v1r2/fits raw/v1r2/driver_logs/v1r2_driver.log` | **byte-identical** (sha256 `f6954e7f38500344…`) |

**3/3 summaries reproduce exactly. Nothing was adjusted.**

## 3. Files (232)

| Path | Role | sha256 |
|---|---|---|
| `conjecture_mode/cards/raw/v1/ARCHIVED_MISSEEDED_FREE_EXCLUDED/PRIMARY_DESI_CMB__free__+0.000__s0.json` | ARCHIVED, EXCLUDED: mis-seeded free-ε fit (v1 ER-B); used for no statistic | `99fe4883c7ae8706bff7840c9badf98ac0ae81c166b4bdf7bb8d9bf1efd57b16` |
| `conjecture_mode/cards/raw/v1/ARCHIVED_MISSEEDED_FREE_EXCLUDED/PRIMARY_DESI_CMB__free__+0.000__s2.json` | ARCHIVED, EXCLUDED: mis-seeded free-ε fit (v1 ER-B); used for no statistic | `20b4ec50273edc27f935a35a2c3e7e66e65b6663dc16873c66c0fc8a50e9aefd` |
| `conjecture_mode/cards/raw/v1/driver_logs/d3_driver.log` | REPRODUCIBILITY SUPPORT: v1 driver log, final attempt | `b47d7dde84e7dcf1d9585e7d4bc41ccd4cf950ae3ce4c91e5117d7941d7150ec` |
| `conjecture_mode/cards/raw/v1/driver_logs/d3_driver_attempt1_oom.log` | REPRODUCIBILITY SUPPORT: v1 driver log, OOM attempt | `d22ae78e4be328c5ef6aaf66438ea1f49bfa35c351c28e8499528faeca39cec3` |
| `conjecture_mode/cards/raw/v1/driver_logs/d3_driver_attempt2_restart.log` | REPRODUCIBILITY SUPPORT: v1 driver log, attempt interrupted by container restart | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `conjecture_mode/cards/raw/v1/driver_logs/d3_driver_attempt3_restart2.log` | REPRODUCIBILITY SUPPORT: v1 driver log, attempt interrupted by container restart | `7c0686daa18c5ba9c24496030cfe2705de14be76000ec55c1d4a93dedd19e390` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__+0.000__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `5278387f3f4d5cedf98544f67b48efd7da7dc3ed05bde036b9558c319488cbef` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__+0.000__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `ddeca97db6a4ea3e00c7108126b63c7af7102455dd7da8c4ca9ce7597b41d2d5` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__+0.000__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `c40792a3e6d92ba3b1c8294fd838a95f36e685431b888f1b37b2602780d76bed` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.005__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `27bdb5c38f03c5f6fe07a9d1f5f4819cd6ebd17e9f972245656456c2b5bc467f` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.005__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `e89a7c1698c6f7f36976d0ee8c67f4cd1c86662bf452541fe5b162f94d87613c` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.005__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `db0c8f260c0d4615737453c74388951739b5fd370f1fcb2c381b3b06787d60e1` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.010__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `b13de6092c62bd8d323ebc6d3211982a50c0e8c64017a684ce5d8bc8ee212496` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.010__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `add2b860f9880341984f126f944f42bbccbfebd6273f40296bd68ca066dfafdb` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.010__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `cec24dcfcd585dadbb04e5bf26c87613b9ad3449498345e4a3ad25b91c593761` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.015__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `4e07a8ed2be9a37b09b069057de30897de52adbef4d72927f1579d43a6e8077d` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.015__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `cabba8097e8ca2b6a0b11f8a9f7e43d1fb06769a6fd77c1ccaed7806eaf7b2ae` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.015__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `8b1c403c0e4d3821324004bdb3b87b1df7fb6c440f21a3b34f54599d3730b530` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.020__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `0fa55aa62f7fe0c03d1668cf2c26391bda719063c474a37132c02154e4000baf` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.020__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `5eafcc37d2954422b71b1c20e7859662d69a96a1bf70ece4d4c069d4026ec8fa` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.020__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `52f93c52664d85e92fdc1629f934eaa2d26e79d257ed0774ea1a900833863e76` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.030__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `9c4c9ea0f3a90f211bb37a7aa7700023d32c4843ffe0daf0f9ec5246d69e640a` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.030__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `72591ad46cf3ddd144e97a33a4b4b9048b931f6714001e3b68b14201eed4235e` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.030__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `f350a2a192eca093f807b5251dd5c3489da280988109ce4c376ceae6cde4c783` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.040__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `097523362962135cb4c7fb21211444f65f829c81d1d17cc7c35865f51328fe27` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.040__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `04c2685a1a2bde85f9a6c2cb8092ea7901fb47dea211b534c7a33ebd70133bac` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.040__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `833ee66ab5e36ddd0951a42ce57a2cf550fe45199a7f1707a805799466adae8d` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.050__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `4215a625df5c405c2b4e3f942b2d65a3ec3046438b6b7cf2786e7dcc53fbc616` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.050__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `b5b9a62626450558260b8b25d78c2bd13dcc1ea147d64543d82edd73951ff519` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.050__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `8f32b5192bcc5b92b422e21ef6325bd6501f8dfe7ed6fc8233f126981656638c` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.060__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `3ef2f65f5af1d42f98ab2d36f0382b7e13abb2bfe73df6a798d6ad0e82b90e38` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.060__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `7c720b360a374528dfc11b5b998d1cd673b560271963dc5c78b04570f9ae77c2` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.060__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `490228c2c3af130ef2f74d51e4e9334a52ad503655ffcaee6c1152fe3a1007e5` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.080__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `b825fd05ac9d99384dd669c34edced323519466b99fb8ae3d45bfb721b28dcb4` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.080__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `91996525c83325b2702c0aa1877b28b10bfc409f41708b6b7e34817cbc813f12` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.080__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `d1a32346699687fa16cb8b3aee36a9743c2caab4316947df001e5d86f28bc975` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.100__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `9758d2942a7e28d8e4aca9effd73080b11d5b15e043fba506a7c4b74ea57e0e7` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.100__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `b392797c3b217bfdd9199b1ef77229e9e1e4417121adffc67c0516ebe0e23bdc` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.100__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `cb146e8293475fdc59a113ac47047ad3b1865a56586015f2d72598d73c2800b4` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.120__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `78fb20bad8e631cfdde2bb48a0ec7c027f40f17b61fbe1409de8177c6757ab29` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.120__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `301e0a22731de82f5580f90840d97d77fd210ef1d50b4656113a9a26b9abcdca` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.120__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `63bb6385cddea2e878fafb000890f05d6b8c2f08395c9edaa623bb91812b0a3b` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.150__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `10184f719cb3c5b421528c707e0623c4499f361815e4fba1e28b096d3f2cb612` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.150__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `34d350cef643b66dde56de12e844d173fdce300e8e0525e2242815f69d798981` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.150__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `415e36cf63c453b1544e0e5e75f2924021f6fb67f8359dd1b8fa4b54ca533325` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.200__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `8c5646616586fd6f407b2d6565d4ac876f8b738c91c2555e94e92879acd53298` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.200__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `61e99b8cbdb7d40c511c30643b60aefb35a8e548a32512cc7467739f9fc131f2` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.200__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `62d30bc12a923ecd5eb07adcf3ac8c85f5b5e7311f440f154f15655094cd5a14` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.250__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `5d13892bb2ef5ec53cd34d744dfb64111b54a5c0ab8fea252ad6b1c2b97858ee` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.250__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `80bc61f7f04acb2fa6b64b2114441f914e7039cb127e8b2db37e826160f7bf89` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.250__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `1fddf229101d6e14babf9f2e84939b78518eb95d74997310d5776712a332cbbd` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.300__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `08d733ee9752f491142161765966a3aea3d99f449a52522e4b8932bcb74bd0f9` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.300__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `8459f88e5e5ccc89665f043c6b5b9f5d65ec0f70f2dab3a439f47cae9c13e5bf` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.300__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `470975ce604f9bc0e9341a3e084783ab35ba0753f9d3cdd81dd624be8f95f9b6` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.400__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `f0ce2d188da92d6f62069ea0b7bf1ec582ff194026d9e1606cbd40692859d203` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.400__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `03431c8f571ba50b55ee33f2d3b600ede940320056ab3c697bdad1d7fa765b66` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.400__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `6d5f5614e3ecfa38626f58cf3d129c98af7793774bbc3a91230e2d203434b53b` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.500__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `5a1ef0c2813bcc7e82df69eb7b3dbe7796f78ec03cf8394678e849fa06c03565` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.500__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `5adebaff04cb75fbe90f71c73ebff7e25e6cd5423b582cfd52ae2bbfdf90ad46` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.500__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `fbb0c2ad96ca346e45cf0beef42b637b763859fde08bba2b954dc8ce738e34c1` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.700__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `33e14a7052aeefa9e4e61936b369f646f8ce804645cba95d701e088b678bf314` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.700__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `cfb01d4f5ecb5fcb2d2025d0511e768bfb68af274158274930a53fcf400d6d27` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-0.700__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `95e98866b8fb557e1668f3b982ee596ec5bcd185d864b308c7559aa0465a250e` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-1.000__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `076801d9305ee480eb6bd125968293752bbdf1e1fc0f772400cc01fed7f7125d` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-1.000__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `7a30d899327422f7588118d288a1e7e522a0338646c6f6a232d73d7bdbfcbb15` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__card__-1.000__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `2265133c6a97154553e45af5455a4ef55b9ee6ff5bd9e935b8712fc9e09e53d2` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__cpl__+0.000__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `a4b90fcb9e78d6f4de453061b4d1691aaf81cd558685bd3ba1b84b5a5dc5fb08` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__cpl__+0.000__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `17de90af2bc96decc4c68529a42a3a80bdabd8502f6dfd82df4f15e505473f8c` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__cpl__+0.000__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `4c43ac3644c9cdf21933602b89d35587b75cc6fad069111d3b909db950b1695f` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__free__-0.080__s0.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `e7b9f6d7d58599e5ff6f1f7beec6601e3af0ec695fce4822b124e8cec6b69f2c` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__free__-0.080__s1.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `4da9d45b9ddc7fa696148570bc2ab167e42f4151968b9c8f80d42663d9dd74e0` |
| `conjecture_mode/cards/raw/v1/fits/PRIMARY_DESI_CMB__free__-0.080__s2.json` | ARCHIVED SCIENTIFIC RECORD: v1 raw per-fit output | `57369e41db917df764ffbca9ac57b435890b208839f7b8484efed115119e79e8` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__+0.000__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9b8207866c9d3288b7a89366e048c21c6b1ace4e429257662e359577bc847d2b` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__+0.000__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c1ab6ca1991441c669eadb0b53c90054aea0de2768aace00795164fd62689f84` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.005__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `2442922dc341877b692472e9a3e059ae9e6658b258789099e832864b3ebad619` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.005__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `d170a1f86cd036514922cb516aabc39c64a9eb3da50899a384d983b40b98fb31` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.005__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `845cd5db3e0fbffdf8d520bf9420e29d6ba16b7d35163fe382f319ee1fe0aeb8` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.010__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `2a3fff3994251c53ec27e5ce8b1cf381185abf2dc699d0755d59ac8dd58478ff` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.010__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c514edd9c66be9b8222de1e9361667472c0461d9b52f2a296c888f0a2df13550` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.010__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `00f86b87711285f6ab2a007e22a300aac7c0b833aa4618abd17b770b128664f4` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.015__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c6893c20c73f7ded96b6bc8705dffb925c0f09c08465d3705e32e8ce2611bd7b` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.015__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `bb1c9c855ddd7595dfaf7222f1d41a5aa63d043baa4d651e3f050c4ac5e57ccf` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.015__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `127411c80d1cc837570b91d54293acb8136f4a0b1f4b4a8b38de925cf832ad22` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.020__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `5389293ed972b69a1bf0b1c31c08edd3a72fc45af3e8b1f7c7464a014887bf70` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.020__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `e1e15f37ad6c9c16f02f6a7c21ad74d62184d9d58e2143c47e36f556a4cf8db3` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.020__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9e9621538973709618818b8ee67d2623defdc087b826330fec59edf5ffab3e68` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.030__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `01e2d7e61a905c91a455817378f53bf79bbd367016d91037da21e0dbbdef897b` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.030__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c64a0ac06f4a11d72ad68a2a30670aac1b3529cdb2f84584f573e46556717d42` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.030__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `d5fa153f79cd554dfec96a4a206f1633355dcb382ed4c2e81a852bc2bf95dc1c` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.040__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9ebde46e713be92b49bcd4ddde65cccfa66c58c4442d80859347844251b73e51` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.040__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `d96b1c0e7bae3d538a490630f98bff7031db3ae783c354a5695f91e6df2c749e` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.040__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `60bd9506961b6378f36cf603ff1acc893202f2b748ee49ac624d176661fac755` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.050__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `56ffec2f900d6370d13906a0d13dd26388553f03f5597d60c01085296b37ea1b` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.050__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `35b87c51a86c262ca7cca462eeaff019b55fa97ad98f89be0bd1818aa6b97e8e` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.050__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `1ab2d6ff433af1d1fd27b19c407447deabe351fbbeca042bf0140d1d319339fa` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.060__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `5f9f34c2fefa205ca816ecc12e6c3866c2ace09ee3e0f3fe074efbad4224780f` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.060__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `b70d922bd32925bcab615cfeff2daeb031d200567e9b2e0ef1d3af40a846464e` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.060__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `12e364fcc4766215f0ca1855e28e5f04c1300208eb04ea230ace68fcd3de1e01` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.080__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `0173714aa1fc4c531bbe6b2853742e0cd0653e930c456033f9e07a5c9553fe51` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.080__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `687753fcc7cdd4a08249b0e9cbf2584c825c76edef66616875eea681c327f50f` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.080__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `f1a91abe7609f2e3dea951a26fcf828c7df9d305a47d3e0c0b4f7a00fdd3cb68` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.100__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `683c39ba33f61d79d9fc78f2ee460dc860e5a048ba53031af15f58335b5a6577` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.100__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `b39ede94769aa4d95fb8055306bb3d5586bb71c5d14023be21f4e93568a8c958` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.100__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `1d5445408593ab4b6bafcf3f0df37c0e60a304467867dc39a63b305304477855` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.120__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `48854f2075796c8538fc06af81b1d291672195705dbf041d0c0214e4039cf84a` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.120__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `f30daafe2f3614d6f5e8badf1a77459bc7607fc16b942feb0da62900b877609f` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.120__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c2084ea0744288a7a2f2ea664f5e4fa956bae9655f47a9a3d1ab7e46bff11d57` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.150__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `3e47e0120fae4d8b316793e5620827436ba5e81e5b15b2b4a36f9e35f2ddb226` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.150__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `18f6285863a31bfa23989f81174767ff6d52ffd8c4d5d20a0025b8c0068a5f00` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.150__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9b826473c0398c69a2b3bcae75248599ea56ac4250e22959d7b5ca4db3078e66` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.200__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `581c8e5e69d7537d1dea55db190e1477bda6090a67bb862ada708a9de17a27c5` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.200__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `a73f22a8f660faeafac44dcfd4f356ea0f0acab4f45689a3778be6351dd13969` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.200__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9e61e9ce3b6deb2c09f0801e418a34c80018453ea1387fe15740c73d30e96e5a` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.250__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `460ac35b2091898f072ee51231e0e2e2e140778651c1d1800166dc7faef3c9fe` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.250__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `ea8fb4f19fb8c630375bc8a1d673119a235ab3e997bd7c894b46ae85a710b1ee` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.250__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `996147b24ad016b0ad4afacb788b92d2586e41b2909f4383889b0270377203bf` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.300__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `db5715d3ddfb55322b9521ab7e5f93e53c7728b61bb732512f965c78d9dbcea1` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.300__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c2518d729907f2a8613ca0d0c625ae3e76717d74ffcb60f790e978511bac0ced` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.300__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `10800a066380774b2f622445b9e550b6daa9f9257a5c4209d334316acb323364` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.400__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `a1e55e39f0cd5cefafb0216185088d99eb2586601a1c859082dd474ba12088f6` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.400__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `a3468f953c10d0a736d989e0068f55138362381f083da79e8d94dfbd1024fd93` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.400__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `306f6e85f93c5ec588725292a760f89456b479a7ab976733ae2ff5b35d804c36` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.500__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `9956f154a40428c85d1fd86572a310993553cbc58fa4da1da98d9a2e662691b3` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.500__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `f53af333cf779e7bb959cddee90df2e6dd3fd8feb2716c7510ede17e0598edbe` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.500__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `97a869a13b60f033e65d60485abed4f40d626bdfb17a0e0f98f80eba35c01f38` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.700__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `140adb8643a976d0f0cbfd8949b00ab897c53330887f06964c46c6542beffba6` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.700__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `0aa3ab1f960908ca0666cbe09f505792b7bd13e7a09e6eff8072fc85cdeba3e6` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-0.700__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `3a5e8986ab143f94cd526ef3f8f2185bc31b9400a2454d20a7c4f296fe615131` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-1.000__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `64c2f884fb6dcbf820d137da7ca96b1ae96417d103d041ebaffb061b1bc58872` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-1.000__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `c224afcdc721c88cceb4719303be9a4a87c8f07dad00c6557775aec857cd4f96` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__card__-1.000__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `eb408b0c7871c3a6a947b7989fde5ce55fb4b36b8c5a1a1e35ce980d4b183db3` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__cpl__+0.000__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `e5223b61d3358f1d82c2f9df3f57c0e90882f12d21088e7dc63042dfdf615644` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__cpl__+0.000__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `b231d5491b7667b413e22c4bc22a575928a62158469cf8c96aa2f1a0a018aa34` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__cpl__+0.000__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `d9cbd9b4b461ef75af9e1d794d3bf42a4169ef122c8346dae833966350104cc6` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__+0.000__s0.log` | REPRODUCIBILITY SUPPORT: log of an ARCHIVED, EXCLUDED mis-seeded free-ε fit (ER-B) | `522492251a9cb23b62b4c997e4de3e2a3a10ba77eadcf1c010723d6b2cad63ae` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__+0.000__s1.log` | REPRODUCIBILITY SUPPORT: log of an ARCHIVED, EXCLUDED mis-seeded free-ε fit (ER-B) | `9d29648d6486015d0cb1f6e5b82bdc2d8cdacea743ea32a86b7d4b1bab75e44e` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__+0.000__s2.log` | REPRODUCIBILITY SUPPORT: log of an ARCHIVED, EXCLUDED mis-seeded free-ε fit (ER-B) | `d41bddbf092806bdf4eef6b022903c806a3fedcc617fbb98f6feedc9d721a125` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__-0.080__s0.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `4fbcc1c767a8dfd6c1fb77378544a06cff27a2d27944bfb87c37f30ce8b64df0` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__-0.080__s1.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `840bde54f1e31647ea86a77273c10534e5fb0ba773533bac296cad6ff3e91042` |
| `conjecture_mode/cards/raw/v1/logs/PRIMARY_DESI_CMB__free__-0.080__s2.log` | REPRODUCIBILITY SUPPORT: v1 per-fit log | `5dc3b1b22d6006f17ae3037e65a94c6eff74cd4f59cf78eaf250169b32e82111` |
| `conjecture_mode/cards/raw/v1/logs/smoke_primary_card0_s0.log` | REPRODUCIBILITY SUPPORT: v1 single-start smoke-test log (ε=0, start 0) | `e25cc37005a0c9f8fd940cc54b7c772800f115553f53f3f483330b604fd611b0` |
| `conjecture_mode/cards/raw/v1r/driver_logs/v1r_driver.log` | REPRODUCIBILITY SUPPORT: v1R driver log, final attempt | `c5d92c2d8a5bbfcaf3a3bf3d598e581d240d1e590e27e57f5904383624f77272` |
| `conjecture_mode/cards/raw/v1r/driver_logs/v1r_driver_attempt1_oom.log` | REPRODUCIBILITY SUPPORT: v1R driver log, OOM attempt | `9d81b290a709b09101e94312cc539f69572b42f0f169ba448b1326abe5b617db` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__+0.000__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `66ee900644f590d3fa7ba87419aa6ad27e99f4cfa7d517d2527731bd5658ad44` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__+0.000__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `946ebcc79729c52b899fb10810498e68887212d7727025d66e69936b393d9fbe` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.005__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `ffa1f11f9e3176cddebce82865487e29e8aef4a8a85e3b7322a5e52081291e0c` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.005__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `643108d26a9e4234f027ab306bd2d066e58d1c5db87c55a3b9e6bf0e6aa56b16` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.010__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `a2a9894ec677f17dd28f71f01515c6f749f5ed9a4e3c9653a43bc6d0440d0b4e` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.015__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `aecbf7b75de761730405eb0d7c88c6314c90c5a6e67ef8b714a45e49b29f6580` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.020__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `64451f33c8ef5861928e7d158a598b2516fb8cc41ca9e9611d2cf8c4c35e43bb` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.030__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `f94028285a365cd1f8eca8868ec3bb504f15cc762b95b290913cd7de2945a587` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.040__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `d2ef80748e6b13d8c08c7bd67bdd5746ef018e2e073b1bd139f712af4fc3f2aa` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.050__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `51fe1c28281d7f826b1607892222a1518c059f0f370ddee2fdd32529e2725a04` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.060__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `1639715c4c581ffe8aab14b2d80ae53bb96657f3b37c78944ee99574913f7d83` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.080__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `9c980f1a5f78fcd2956a6f0332dfc404ccbf03f6b3a634c880001ad6d4b135ed` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.080__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `461825819be9c3a2a4f2f286ca506d3c3ab454e21ecd00716771542d10119510` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.100__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `c48790b7afad988e04ddef57b6a9466a96fa702e6568b2db5406a800a9f00ee4` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.120__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `d47f1ce2f529d7fdf38c3ae2de4218da516303c4fce6cc305ae3fb31e7494338` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.150__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `150435d50b63ecf4e3782861cbebe20a793dff54643aa9a10e9f6123b4e8ecdc` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.150__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `9a233b6f6e3f7b0259138d0bf98a3196478062e5cce89e8fe7cb300be35cf4ae` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.200__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `becac40a2bde805b10ab16a5c69f06b5270f269d00b99f50deefe4cb43a2bfdd` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.200__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `208db6ce00b77cba64c7e114d3e95e1eb0638ce6e8eaf3b51131da5da7deae91` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.250__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `55ea244877a315d7705118bd28dc4e47d99790b9059e27fa020c177998cf80a6` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.300__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `989feb4daf50a63d88883868039b83e7b48dd33c9846f4afb0fd806a04ce7efa` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.400__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `d14f333a9a4cf4f7fb1e7007e02bcce5ef6908c39d517943ba14144f013f5ac1` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.500__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `1009814e6fa1596c4f44412d74a4683bdaf3d1840cf37bf059b9ea1da6be2493` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-0.700__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `6a15ede9ce780b21a982228f3e84decfa2be50e54ef4c46f9743727cbe215a60` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__card__-1.000__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `4868368fa33bed685054fe8ba4a2db67e74c6110066f9616d9937c7d1e3f766d` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__cpl__+0.000__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `3f0d711be8e79d8d0b1ef18987356c27a0af6fdd96fa58db23146cc51beb53e8` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__cpl__+0.000__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `fe37e2de45255fe78e8805eb809dd3e11779cb24081fbe50c4f89cc02e321770` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__free__+0.000__A.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `2e35d3ebf731caef6944fcec375604806bc75912e59e3a7a14ada61f96fb5e72` |
| `conjecture_mode/cards/raw/v1r/fits/V1R__free__+0.000__B.json` | ARCHIVED SCIENTIFIC RECORD: v1R raw per-fit output | `cfc76e8e293efe46069b4e66d39c089c231d4747fa20532f3c9275dbe18697b1` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.005_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `4c2407dda9a59923c12d5a717c75e5525c5fef0d20fcd4fe684b7ca77bedc943` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.005_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `ea3fbff28738c046228341c5c5922d84a403f86d51f15acb8ed632eff6fa395f` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.015_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `ce3730411f736f89f4bccafdc478b5f9f4c2cf870ff65330005a8ec399eff92f` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.01_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `458cbc01e12e5e7e3d4f40f1e02142c04a6a5b5d3f54170091d2f6edae95bc3f` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.02_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `12472f5ad84840c97ab99170d7a5644e8270170e95aab87963b1cac5d9ddfd2b` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.03_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `c7f1d90e04131f599b4af5d12906f0e9bb7418166f70087f2c8eed9345e1df27` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.04_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `b3c7a1343417022b787ed4a7e24b1792ec21a98a5280725c2c71fbf2dcc11e59` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.05_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `f6a7399308a5db5ab60bc3ec50506a58f275474a09ea8f60e612d860f9fbcabd` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.06_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `85ce1454d84f7f23d7ef2c50989746fd72ddfe78e2703787f06bdbf21de6d512` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.08_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `0c067ced002513a63919df767be96289177df70dfdbb026a48b9c5d9ce666bc2` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.08_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `1355409d38727bdf5941156932f84dee1b3b981d09f37eec496c4b9e31161f01` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.12_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `69a51ad327ab3dde5a2de9e51101060a72aaa6942dcd1deb218d8d0ac226a231` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.15_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `70c0f9d4025fb94982b1c73e24ef1f520d10ebc77557a96f8e96d2d47700a78c` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.15_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `2af4b94b49b930d360877f371217d82248d8b3bdae31215e2909ea94e296f821` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.1_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `3cf09424c65f9cce653dce7fc08048ceb1d7da3d35ccfaf7009a0645a7a2cc6f` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.25_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `a97cffd2fcb4ae005486b2e7259f36b0f6d7cbb732d3261506a27b1bb9014c66` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.2_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `59150345a52627d9becc9eb90451a786de1135969733e5a85accf6f1d510fe1c` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.2_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `d5b076dfdd69228b28148094d7308cd82c2f54256b1444e543d7a71794ba3b7d` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.3_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `d8b2aa127c0498a4425a576aaf2d31c9709fbfe42ec50e28e88b12d73f3710d2` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.4_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `636bde96a9e967defebfaa53db77fca53218fc33ec14f08d44f547c1313d2d28` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.5_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `a3cb6896f23f4164a9ccb3531e859b15abe00063ec9fe2a4e5de2bb96ec40ec6` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-0.7_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `ed94e9f46d5381f29a5798fcedac819325c9430d6202b82bc8480072eb8afce1` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_-1.0_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `65a2287d1cc45dc5bd95e7ee79dc888cf7210a43d6f6f8e48cfe53b99cc59ffe` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_0.0_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `d634443d2780f59508115e3ef5acb8fd9d9d207eddd2f35aad4c11b851cd90c0` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_card_0.0_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `537c839f9b486100c117d78797ef9f690c8f23d46e29d2806ca8c5f6162226da` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_cpl_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `95e2b0ffcff6840f28597a4fa8d349e22dee0e4ffed62e0dc7e6fca3e43e4438` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_cpl_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `16d6628b7a41306f5984b01baa987ad4c65f9887a36b20cb36d3cdea3733c7af` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_free_A.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `32b5f8f011642ec2cad508146a2f13491e0e96d19a298c412625df1b891e11f9` |
| `conjecture_mode/cards/raw/v1r/logs/v1r_free_B.log` | REPRODUCIBILITY SUPPORT: v1R per-fit log | `d913172e8942eb3e231dd50f4873c37ce377a02d0e775ab304a4119ec6eb59b7` |
| `conjecture_mode/cards/raw/v1r2/driver_logs/v1r2_driver.log` | REPRODUCIBILITY SUPPORT: v1R2 driver log, final attempt | `7da78fb0f57b190a3573ac41d5df6bf1731656f303f1d94ec995815042402c76` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_-0.08__r2.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `9ba1c122f2dc569ba77a57699221925a192c2bc7ac52e294ae1d5b12ea893f8a` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_-0.08__r3.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `22a457b9cfbefe2977f425fdaeca84cb6a3485d123af67e0ba2ff05cb9fb0124` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_-0.08__r4.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `46aa1bdafabc746869f459fc708a40139e65a461f5dc63f76fea5c7231ef9ea5` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_0.0__r2.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `4035de5d322713b66f08b7afca1ed5e9ecbb14a2abf397976460638baa3df07e` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_0.0__r3.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `7bacb7f75e12a1aa0bf709be63821d60aa92da85e13c605f70ec7431975d135f` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__card_0.0__r4.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `b2225af2c8ec1eba1349977027a069451e0b6ae703f12e0cad02398e3aaa0e18` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__cpl__r2.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `31d03af95ec843eb9bf181580158aa60a67da9fdf071531f6f7230d098f67ef0` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__cpl__r3.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `d788a0b130a2943c5f1fb09769c1da331a5bb5bf879d3327c72f7cbbc3ee59b1` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__free__r2.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `42caae4e9748ea1f4a17ef06d5fc90a932b9237ffdbda4957037e13bc454fc0e` |
| `conjecture_mode/cards/raw/v1r2/fits/V1R2A__free__r3.json` | ARCHIVED SCIENTIFIC RECORD: v1R2 raw per-fit output | `9129ce54f985f37140461dfeb5445f74ba017a9d719533a986345679dee42fa3` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_-0.08__r2.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `3af4dac4701e7b983b91d9f74a1dd9c9e50c76f4dc01e2b73996802126667e7a` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_-0.08__r3.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `88ae92f038f16973dbf6ee3e37d9e87a691cece8ba484559b88787df874ea95e` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_-0.08__r4.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `f1f65cb55035c362a01ba102c570690dc03c3cba0b73b276c05daa3f6e8b87c4` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_0.0__r2.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `a88595643592d56242e5c61f1f7788af08e90277de36dc7b4ee87eead2733232` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_0.0__r3.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `7e5c2daa0d065ae4bece693efebbe91078ca2b3fe5c4baff2ea14819ebd0cb0e` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__card_0.0__r4.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `fee07acd76784909e0745b555a6c686baf6ead5efcb49e4a2f2f765e859e3750` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__cpl__r2.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `2616c5fe1a2f6eb0725137c1c0c56e35b1e14ec8b6bfa819e551b0a7e1cf7703` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__cpl__r3.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `c7fbb5d44c953ecacd36fbe1ec8f360b04dc3aa33eea439ae84b1fa4485468c1` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__free__r2.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `4f4804ff6b49601389534eaa1d40e938b5c99476f04cde646413adc328946bca` |
| `conjecture_mode/cards/raw/v1r2/fits/starts/V1R2A__free__r3.start.json` | REPRODUCIBILITY SUPPORT: v1R2 start point | `957ace544920a7bf8696354c072ff771f93501136eb7d12c9220adac6191db46` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_-0.08__r2.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `2360bd8d8d251f5e79566abfed44533ceffc9e991a612cc22997e7cc2079a74c` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_-0.08__r3.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `a5e4baeaad8442ce0f071ee630b2d8667a16e10f7b17d3d35162f7de401b353a` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_-0.08__r4.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `cac49a1fd76e277545e43e6ac4c4c5060e31c6684503b7ba43f79c8d9396f5bf` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_0.0__r2.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `c333bba347e07909ccc7de35a780c6323883cbeaebe8ae79712486ffd8e5c121` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_0.0__r3.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `048c3480af571fbeda04316f763953ef3bdfb78836320557eb356226347cc441` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__card_0.0__r4.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `cdad06ef5cb2d02a0232f9b683978601a71ceacb42cbe3a27a47b973235e67e2` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__cpl__r2.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `86b1e5fd8e58843f694ab4ccb47e10b9bb9a9b0dedd7eca67b30b6c183901fd0` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__cpl__r3.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `c06f144f274904ef866b2e1650e83b87a33621a327bd66c059d570701f660c8a` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__free__r2.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `881dffb412d1bd61230aa831114a46784f213fda5ab9a58e817cfa1a9c612d2f` |
| `conjecture_mode/cards/raw/v1r2/logs/V1R2A__free__r3.log` | REPRODUCIBILITY SUPPORT: v1R2 per-fit log | `13f5031a9379875f9ece570bf30e19b4a13b47d8c7fbb962bfff5ebbdac69206` |
