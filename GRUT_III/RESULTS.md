# GRUT III — verification results

Generated 2026-10-11 13:51 UTC by `GRUT_III/calc/run_all.py` (Python 3.13.16; numpy, sympy, mpmath).

| script | status |
|---|---|
| `v1_ds_thermal_kernel.py` | **PASS** |
| `v2_jacobson_eos.py` | **PASS** |
| `v3_horizon_fluid.py` | **PASS** |
| `v4_medium_absorptive_kernel.py` | **PASS** |
| `v5_anomaly_selfconsistency.py` | **PASS** |
| `v6_low_acceleration_scale.py` | **INFO** |

---

## `v1_ds_thermal_kernel.py`

```text
== 1.1 stationarity identity (sympy) ==
  a a' (eta-eta')^2 - 4 sinh^2(H dt/2)/H^2 = 0  -> PASS
== 1.2 worldline Wightman spectrum: contour-shifted quadrature vs Planck formula ==
  [PASS] G+(w=0.1): got 0.0341159443205  want 0.0341159443205  rel.err 0.0
  [PASS] G+(w=0.5): got 0.0831716441794  want 0.0831716441794  rel.err 0.0
  [PASS] G+(w=1): got 0.1594527119  want 0.1594527119  rel.err 1.55e-31
  [PASS] G+(w=2): got 0.318310996243  want 0.318310996243  rel.err 0.0
  [PASS] G+(w=5): got 0.795774715459  want 0.795774715459  rel.err 1.17e-17
  [PASS] G+(w=-0.5): got 0.00359417263343  want 0.00359417263343  rel.err 1.07e-31
  [PASS] G+(w=-2): got 1.11005941965e-6  want 1.11005941965e-6  rel.err 0.0
  symmetrised S(0.1) [record normalisation, flat line w/2pi] = 0.052316394  (= (w/2pi) coth(pi w/H));  record's quadrature value 0.03393; exact floor S(0)=H/(2pi^2)=0.050660592
== 1.3 image-sum identity ==
  [PASS] sum_n 1/(z+i n beta)^2 vs (pi/beta)^2/sinh^2(pi z/beta): got 1.95948563096  want 1.95948563096  rel.err 1.01e-31
== 1.4 Wick contraction: <T00 T00>_c and tracelessness for the improved conformal scalar ==
  <T00(t) T00(0)>_c = 1/(4*pi**4*t**8)   expected (3/4)C_T/t^8 = 1/(4 pi^4 t^8): PASS
  sum_mu <T_mumu(t) T00(0)>_c = 0  -> PASS
  tracelessness at a generic point: 0  -> PASS
== 1.5 stress-tensor Wightman spectrum on the dS worldline: quadrature vs closed form ==
  [PASS] G+_T(w=0.05): got 2.14349834831e-5  want 2.14349834831e-5  rel.err 4.21e-31
  [PASS] G+_T(w=0.5): got 8.21647769999e-5  want 8.21647769999e-5  rel.err 1.46e-31
  [PASS] G+_T(w=1): got 0.000320554315569  want 0.000320554315569  rel.err 1.5e-31
  [PASS] G+_T(w=3): got 0.022460890198  want 0.022460890198  rel.err 1.37e-31
  [PASS] G+_T(w=8): got 8.25895246049  want 8.25895246049  rel.err 7.64e-30
  [PASS] G+_T(w=-1): got 5.98616826727e-7  want 5.98616826727e-7  rel.err 0.0
  [PASS] G+_T(w=-3): got 1.46274573912e-10  want 1.46274573912e-10  rel.err 1.57e-31
  [PASS] KMS ratio G+_T(1)/G+_T(-1) = e^(beta w): got 535.491655525  want 535.491655525  rel.err 1.89e-31
  [PASS] KMS ratio G+_T(3)/G+_T(-3) = e^(beta w): got 153552935.395  want 153552935.395  rel.err 0.0
  [PASS] flat limit G+_T(50) / (w^7/10080 pi^3) -> 1: got 1.0056078423  want 1  rel.err 0.00561
  [PASS] flat limit G+_T(200) / (w^7/10080 pi^3) -> 1: got 1.00035003063  want 1  rel.err 0.00035
  [PASS] classical floor G+_T(0) = H^7/(560 pi^4): got 1.83321111691e-5  want 1.83321111691e-5  rel.err 1.64e-31
  [PASS] noise-kernel decay rate (tau=10..12) -> 4H: got 4.00015702692  want 4.0  rel.err 3.93e-5
  dissipative part Im chi_T(w) propto G+(w)-G+(-w) = w(w^2+H^2)(w^2+4H^2)(w^2+9H^2)/(10080 pi^3): an odd polynomial
  -> the retarded response at a point is a finite sum of delta-function derivatives (no memory: Huygens);
  -> all memory sits in the NOISE, with derived correlation time 1/(4H) and Planckian occupation at T=H/2pi.

V1 OVERALL: PASS
```

## `v2_jacobson_eos.py`

```text
== 2.1 Clausius bookkeeping ==
  Clausius dQ = T dS  =>  T_kk = R_kk*eta*hbar/(2*pi)
   matches (hbar eta/2pi) R_kk: PASS
  with eta = 1/(4 G hbar):  T_kk = R_kk/(8*pi*G)  -> PASS
== 2.2 tensor lemma: null-cone data determine X up to the metric ==
  singular values: [3.992e+01 3.079e+01 2.319e+01 1.863e+01 1.149e+01 1.018e+01 8.444e+00
 5.947e+00 4.935e+00 6.192e-15]
  one-dimensional nullspace, spanned by:
 [[-1.  0.  0. -0.]
 [ 0.  1. -0. -0.]
 [ 0. -0.  1. -0.]
 [-0. -0. -0.  1.]]
   = eta_ab: PASS
== 2.3 Bianchi closure on a spatially flat FRW example ==
  d/dt [ f - R/2 ] = 0  -> f = R/2 + const: PASS
  => G_ab + Lambda g_ab = 8 pi G T_ab, Lambda = integration constant (not a response of the medium).

V2 OVERALL: PASS
```

## `v3_horizon_fluid.py`

```text
== 3.1 horizon fluid transport ==
  eta  = c^3/(16 pi G) = 8.0313e+33 Pa s
  s    = k_B c^3/(4 G hbar) = 1.3213e+46 J K^-1 m^-2
  eta/s = 6.0783e-13 K s   vs hbar/(4 pi k_B) = 6.0783e-13 K s  -> PASS
== 3.2 Planckian memory time ==
  T_dS(H0) = 2.655e-30 K ;  tau = hbar/(2 pi k_B T) = 4.5782e+17 s = 1/H0 = 4.5782e+17 s -> PASS
  BH M=1 Msun: T_H=6.170e-08 K, memory time 1.9703e-05 s = 4GM/c^3 = 1.9703e-05 s
  BH M=10 Msun: T_H=6.170e-09 K, memory time 1.9703e-04 s = 4GM/c^3 = 1.9703e-04 s
  BH M=1e+06 Msun: T_H=6.170e-14 K, memory time 1.9703e+01 s = 4GM/c^3 = 1.9703e+01 s
== 3.3 de Sitter static-patch identities (hbar=c=k_B=1, exact) ==
  E = rho_Lambda V = 1/(2*G*H),  T S = 1/(2*G*H),  F = E - TS = 0 -> PASS
  |E_Komar| = 1/(G*H),  (1/2) N T = 1/(G*H) -> PASS

V3 OVERALL: PASS
```

## `v4_medium_absorptive_kernel.py`

```text
== 4.1 real scalar ==
  m=0.0: Sum_pol|eps.T|^2 = 0.1666666667  vs 8p^4/3 = 0.1666666667  [PASS]
  m=0.2: Sum_pol|eps.T|^2 = 0.1176000000  vs 8p^4/3 = 0.1176000000  [PASS]
  m=0.4: Sum_pol|eps.T|^2 = 0.0216000000  vs 8p^4/3 = 0.0216000000  [PASS]
  => threshold law beta^5 (beta = 2p/M); Im Pi_S = (1/2)(1/2)(kappa^2/4)(beta/8pi)(8p^4/3)/5 = kappa^2 beta^5 s^2/(3840 pi) = G beta^5 s^2/120
== 4.2 Dirac fermion ==
  massless: Sum|M|^2_F / (1/2 Sum|M|^2_S) = 6.0000000000  (expected 6 = C_T^Dirac/C_T^scalar)  [PASS]
  threshold shape: F_F(r) := beta * Sum|M|^2_F(m)/Sum|M|^2_F(0)  vs  beta^3 (1 + 8r/3), r = m^2/s
   r=0.0025: numeric 0.9916044798  formula 0.9916044798
   r=0.0100: numeric 0.9656868362  formula 0.9656868362
   r=0.0400: numeric 0.8519924732  formula 0.8519924732
   r=0.0900: numeric 0.6348800000  formula 0.6348800000
   r=0.1600: numeric 0.3081600000  formula 0.3081600000
   r=0.2025: numeric 0.1275413831  formula 0.1275413831
   r=0.2401: numeric 0.0129257942  formula 0.0129257942
  max |numeric - beta^3(1+8r/3)| = 2.22e-16  [PASS]  (HLZ 1999 form reproduced from scratch)
== 4.3 massless vector (photon pair) ==
  (1/2)Sum|M|^2_gamma / (1/2)Sum|M|^2_S = 12.0000000000  (expected 12)  [PASS]
== 4.4 Im Pi(s) = c_eff(s) G s^2 : Standard-Model staircase ==
  sqrt(s)           c_eff      120 c_eff   (c_SM(45 Weyl) = 283/120 = 2.358 at high energy)
    1.0e-12 GeV     0.1000      12.00
    1.0e-10 GeV     0.1493      17.92
    1.0e-09 GeV     0.1748      20.97
    1.0e-06 GeV     0.1750      21.00
    1.0e-03 GeV     0.1750      21.00
    2.0e-03 GeV     0.2123      25.47
    1.0e-02 GeV     0.2246      26.95
    3.0e-01 GeV     0.2490      29.88
    1.0e+00 GeV     0.2935      35.22
    3.0e+00 GeV     1.2399     148.79
    1.0e+01 GeV     1.7219     206.63
    2.0e+02 GeV     2.0480     245.77
    1.0e+03 GeV     2.3363     280.35
    1.0e+04 GeV     2.3581     282.97
  at LIGO frequencies (sqrt s ~ 4e-22 GeV): c_eff = 0.1000 (photon only; neutrinos closed by their masses)
  on-shell flat-space absorption of a GW: Im Pi(s=0) = 0 exactly.  FRW estimate Gamma ~ c_eff t_P^2 H0^4/(2 omega) = 5.2e-162 s^-1 at 100 Hz
  => the medium's absorptive response is derived, finite-memory-with-thresholds, and unobservably weak for GWs: consistent with GW170817 trivially.

V4 OVERALL: PASS
```

## `v5_anomaly_selfconsistency.py`

```text
  SM (45 Weyl)                : a = 2.76528  c = 2.35833  (120c = 283.0, 360a = 995.5)  C_T = 40c/pi^4 = 0.968424
  SM + 3 nu_R (48 Weyl)       : a = 2.81111  c = 2.43333  (120c = 292.0, 360a = 1012.0)  C_T = 40c/pi^4 = 0.999222
  one real conformal scalar   : a = 0.00278  c = 0.00833  (120c = 1.0, 360a = 1.0)  C_T = 40c/pi^4 = 0.003422
  scalar checks: -3aH^4/2pi^2 = -H^4/240pi^2 and C_T = 1/(3 pi^4): [PASS]
== 5.2 self-consistent de Sitter from the anomaly alone ==
  SM (45 Weyl): H_sc = sqrt(pi/a) M_P = 1.066 M_P = 1.301e+19 GeV
  SM + 3 nu_R: H_sc = sqrt(pi/a) M_P = 1.057 M_P = 1.291e+19 GeV
  observed bound H_inf < 4.7e+13 GeV (r<0.036)  =>  anomaly coefficient needed a > 2.12e+11  (~1.2e+12 vector species)
  => anomaly-driven self-consistent dS is Planckian at SM content: GRUT III inherits Starobinsky's conclusion that the R^2 scale is supplied.

V5 OVERALL: PASS
```

## `v6_low_acceleration_scale.py`

```text
H0 = 67.4 km/s/Mpc: cH0 = 6.548e-10 m/s^2;  a0_obs/cH0 = 0.183
   cH0/(2 pi)                       = 1.042e-10  ratio to a0_obs = 0.868
   cH0/6 (Verlinde)                 = 1.091e-10  ratio to a0_obs = 0.909
   c sqrt(Lambda/3)/(2 pi)          = 8.626e-11  ratio to a0_obs = 0.719
   cH0/(2 pi) * sqrt(Omega_L)^-1    = 1.259e-10  ratio to a0_obs = 1.049
H0 = 73.0 km/s/Mpc: cH0 = 7.092e-10 m/s^2;  a0_obs/cH0 = 0.169
   cH0/(2 pi)                       = 1.129e-10  ratio to a0_obs = 0.941
   cH0/6 (Verlinde)                 = 1.182e-10  ratio to a0_obs = 0.985
   c sqrt(Lambda/3)/(2 pi)          = 9.342e-11  ratio to a0_obs = 0.779
   cH0/(2 pi) * sqrt(Omega_L)^-1    = 1.364e-10  ratio to a0_obs = 1.137
Discriminating prediction if a0 is set by the apparent-horizon temperature: a0(z)/a0(0) = H(z)/H0
   z=0.5: H(z)/H0 = 1.322   (Lambda-horizon alternative: 1.000)
   z=1.0: H(z)/H0 = 1.790   (Lambda-horizon alternative: 1.000)
   z=2.0: H(z)/H0 = 3.032   (Lambda-horizon alternative: 1.000)
   z=3.0: H(z)/H0 = 4.566   (Lambda-horizon alternative: 1.000)
Status: CONJECTURE. Must survive: CMB acoustic peaks, cluster masses, lensing, BBN (unaffected), the RAR scatter. Not derived here.

V6 OVERALL: INFORMATIONAL (no pass/fail)
```
