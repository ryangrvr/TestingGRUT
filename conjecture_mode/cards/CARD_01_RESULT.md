# CARD #1 — RESULT (v1)

## **CARD-01-ACCESS-BLOCKED**

| | |
|---|---|
| Pre-data spec | `8c634145178eaad29153934b0d2846c6628bdd07` |
| Owner threshold ruling | `86ccbb2ec873e9f22a894df035c2249a7c92d96d` |
| Run-configuration freeze | `c9899abc80ba607f431545453b7bba7781e1ad8a` |
| Preflight | 2026-10-04T12:04Z (`CARD_01_NETWORK_PREFLIGHT.log`, `CARD_01_DATA_PROVENANCE.md`) |

### Why

The owner-approved primary, verdict-bearing combination (S6) requires the **updated Planck + ACT DR6 lensing likelihood v1.2**.

- Its data are distributed only through NASA LAMBDA (`lambda.gsfc.nasa.gov`), with NERSC (`portal.nersc.gov`) as the other collaboration host.
- Both are **denied by this environment's egress proxy** (`CONNECT tunnel failed, response 403`).
- The collaboration software repository (v1.2.1) contains no data.

Under the owner ruling (step 3) and Charter §4B, a required D3 component that fails access ⇒ **ACCESS-BLOCKED; stop rather than improvise.** Every SN extension also uses the primary CMB, so all three extensions are blocked by the same component.

### What was and was not done

- **Not done:**
  - no product acquired;
  - no likelihood evaluated;
  - no D3 profile, CPL comparator or control run;
  - no D1/D2.

  No Card #1 likelihood value, data value or chain exists or was inspected.
- **Unchanged:**
  - the postulate, projection, primary branch;
  - thresholds S1–S8 and the S5 decision logic;
  - the S6 datasets and the run configuration.

  These remain exactly as committed. The card can be resumed **unchanged** once access exists.
- **No substitute used.** In particular, none of the following was used:
  - Planck-only lensing in place of Planck + ACT;
  - a compressed CMB;
  - any third-party mirror.

  Any such change would need a prior owner ruling and would make this a different card version or a C2-R downgrade.

### Other access findings (non-decisive)

- **D1/D2 diagnostics: UNAVAILABLE.** The official DESI hosts (`data.desi.lbl.gov`, `www.desi.lbl.gov`) are blocked. Under O3 no mirror chains may be used. This does not by itself block D3.
- **O3 #4 traceability:** the DESI DR2 BAO files could not be cross-checked against the DESI release, because the DESI hosts are blocked. The cross-check could be done against the DR2 paper metadata once acquired.
- **GitHub archive tarballs** for `bao_data` v2.6 and `sn_data` v1.8 return 403, while git transport to the same repositories works. Retrieving the pinned tag commit by git is the same content at the same commit. Whether this qualifies under O3 #2 is an **owner question** for the resumed run; it was not exercised.
- Planck native-data release assets (CamSpec NPIPE, low-ℓ TT, low-ℓ EE) are reachable (HTTP 200).

### What would unblock it (owner decision; environment setting)

Allow the host **`lambda.gsfc.nasa.gov`** in the environment's network policy. This is required for D3. Optionally also allow:
- `data.desi.lbl.gov` and `www.desi.lbl.gov` (official DESI chains for D1/D2, and DESI provenance cross-checks);
- `portal.nersc.gov`.

After that, the same frozen pipeline can be rerun from the preflight step, with the owner ruling unchanged.

### Firewall

- Card #1 remains POSTULATED. Its v1 test has **not** been performed.
- ACCESS-BLOCKED is **not** evidence for or against the card.
- The historical C2-R reconnaissance grade is unchanged.
