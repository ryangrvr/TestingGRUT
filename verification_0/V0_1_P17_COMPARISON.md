# V0-1 — COMPARISON PHASE (P-17; the original was unsealed only after the reproduction was committed at `64b0182`)

**Original:** `playground/SCOUT_0/probes/P17_RESULT.md` (blob `5fc0401b…`) on frozen `scout-0 @ ab2da47`.

**Reproduction:** `V0_1_P17_REPRODUCTION.md` and `code/v0_1_*`. These were produced by a context-isolated sub-agent from
the spec alone. The independent derivation is **not** edited below.

## 1. Line-by-line scientific comparison

| item | original | independent reproduction | classification |
|---|---|---|---|
| Class 𝓗 / assumptions | harmonic modes; coupling linear in x_j and q, with counterterm; any smooth confining V | the same, derived from the Hamiltonian. It also makes explicit that the counterterm removes the potential renormalisation | **exact agreement** |
| Bath elimination | variation of constants plus one integration by parts | the same steps, re-derived and symbolically checked for arbitrary paths | **equivalent derivation** (standard; KNOWN: Ford–Kac–Mazur / Zwanzig) |
| Memory kernel | γ(t) = Σ c_j c_jᵀ/ω_j² cos ω_j t | identical; γ(0) equals the counterterm matrix | **exact agreement** |
| Free force | F(t) = Σ c_j{[x_j(0) − c_jᵀq₀/ω_j²] cos ω_j t + p_j(0)/ω_j sin ω_j t} | F(t) = Σ c_j[ξ_j cos ω_j t + (η_j/ω_j) sin ω_j t], with ξ_j = x_j(0) − c_jᵀq₀/ω_j² and η_j = p_j(0) | **exact agreement** |
| Preparation (P) | the free-force law is independent of (q₀, p₀); shifted-Gibbs satisfies it | the same, plus: with distinct frequencies this is equivalent to a q₀-independent law of (ξ, η); conditioning joint Gibbs gives shifted-Gibbs | **agreement** (the reproduction is slightly more explicit) |
| Product preparation | (P) fails; slip −γ(t)q₀; "noise vs memory filing is a bookkeeping convention" | slip −γ(t)q₀ confirmed exactly. It notes that a re-integration to a q-memory form gives a q₀-independent forcing law | **agreement**: the same point as the original's bookkeeping sentence |
| Equality of path laws | q = Φ(q₀, p₀, F) is the same measurable map in A and B; Volterra well-posedness; energy bound | the same structure, with regularity stated: Picard, plus an energy / Grönwall argument | **equivalent derivation** |
| Thermal law of F | Gaussian, mean 0, covariance Tγ(t − s) | identical. Also the random-phase fourth cumulant per mode is −(3/2)T²c⁴/ω⁴ (non-Gaussian) | **agreement** (the reproduction adds the explicit cumulant) |
| "Matched covariance suffices" | only for Gaussian preparations (H3) | only among Gaussian laws | **agreement** |
| Realisability converse | stationary Gaussian F with S_F supported where J > 0 has a harmonic parent with mode-dependent T(ω) | B → A holds only for laws supported on the bath's cos / sin span (finite bath) | **consistent**: the original's converse is restricted to stationary Gaussians with support condition, matching the reproduction's restriction |
| Ontology conclusion | reduced data identify the forcing **law** and the memory kernel, never the ontology, within 𝓗 | the same non-identifiability conclusion. **Qualification CL-1** on the positive side-claim: identifying γ needs M and V known (the reproducer's sufficient condition adds data across initial momenta) | **agreement on the theorem**; CL-1 is a **scope clarification** of an unqualified table entry |
| C-B interpretation (v) | the discriminator reads forcing-law cumulants; exact C-B only at WB + OD | **not checkable from the spec** (no C-B model definition) | **unresolved in V0-1**, deferred with C3 / C4 |
| Exact Taylor checks (C2), T1 – T8 | the logged values (2 modes n ≤ 8; 1 mode n ≤ 12) | **every value reproduced exactly**, with a distinct implementation (independent power-series recursions for A and for the GLE, with Wick pairing) | **exact agreement** |

**Implementation independence.**
- The original computes A by repeated Lie derivatives of the Hamiltonian, followed by a monomial expectation. The
  reproduction uses coefficient recursions in sparse polynomial rings.
- **Shared by mathematical necessity:** the GLE comparator's covariance identity Cov(F⁽ᵏ⁾(0), F⁽ˡ⁾(0)) = T(−1)ˡγ⁽ᵏ⁺ˡ⁾(0).
- No code was copied. The orchestrator re-ran the reproduction scripts, and all logs reproduced byte-identically.

## 2. Differences recorded

- **CL-1 (scope clarification, owner may reclassify as a correction).** The original's table states, without
  qualification, that "what *is* identifiable is the forcing law … and the memory kernel". Identifying γ from reduced
  data presupposes the system part (M, V). The reproducer gave a sufficient condition (variation of initial momenta);
  necessity was not shown. The **non-identifiability theorem itself is unaffected.**
- **No numerical discrepancy.** No sign, coefficient or assumption error was found.

## 3. Grade

**V0-1-A — P-17 INDEPENDENT DERIVATION REPRODUCED — VER-I1 (orchestrator-exposed).**
- The theorem, its assumptions and scope, the kernel, the free force, the preparation condition, the product-preparation
  slip, and **all exact C2 Taylor checks** agree.
- **Recorded alongside the grade:**
  - scope clarification CL-1;
  - C-B interpretation (v), plus the C3 / C4 checks, deferred (they need the external C-B model definition).
- **Qualifiers:**
  - **VER-I1, not VER-I2.**
  - **Orchestrator-exposed:** the orchestrating agent had read the P-17 proof earlier in the session. The reproduction
    was done by a context-isolated sub-agent from the spec only, and its file-access report is in the ledger.

**Owner option:** if CL-1 is judged a substantive scope correction rather than a clarification, the grade becomes V0-1-C
(reproduced with correction). The theorem is unaffected either way.
