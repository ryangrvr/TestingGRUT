# GRUT III — verification results

Generated 2026-10-11 19:26 UTC by `GRUT_III/calc/run_all.py` (Python 3.13.16; numpy, sympy, mpmath).

| script | status |
|---|---|
| `v1_ds_thermal_kernel.py` | **PASS** |
| `v2_jacobson_eos.py` | **PASS** |
| `v3_horizon_fluid.py` | **PASS** |
| `v4_medium_absorptive_kernel.py` | **PASS** |
| `v5_anomaly_selfconsistency.py` | **PASS** |
| `v6_low_acceleration_scale.py` | **INFO** |
| `v7_access_geometry_memory.py` | **PASS** |
| `v8_residual_framework.py` | **PASS** |
| `v9_induced_gravity_matching.py` | **PASS** |
| `v10_graviton_sector_memory.py` | **PASS** |
| `v11_quasi_ds_kms_defect.py` | **PASS** |
| `v12_two_temperatures.py` | **PASS** |

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
  symmetrised S(0.1) [flat line w/2pi] = 0.052316394  (= (w/2pi) coth(pi w/H));  exact floor S(0)=H/(2pi^2)=0.050660592
  one-sided G+(0.1) = 0.034115944;  the record's calc/worldline_reduction.py 'S(w)' is 2 Re Int_0^inf e^(iw tau) g = G+(w),
  so its quadrature value 0.03393 is to be compared with G+(0.1): relative gap (G+ - record)/G+ = 0.545%.
  The identification 'record value = exact Planck G+' holds only if the record's eps-extrapolated quadrature precision
  is coarser than 0.545%; the record itself quotes <3% validation against the flat half-line result, so it is.
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
  [PASS] noise-kernel asymptotic decay rate (tau=10..12) -> 4H: got 4.00015702692  want 4.0  rel.err 3.93e-5
  dissipative part Im chi_T(w) propto G+(w)-G+(-w) = w(w^2+H^2)(w^2+4H^2)(w^2+9H^2)/(10080 pi^3): an odd polynomial
  -> the retarded response at a point is a finite sum of delta-function derivatives (no memory: Huygens);
  -> all memory sits in the NOISE, with derived ASYMPTOTIC DECAY RATE 4H (large-|tau| tail e^{-4H|tau|}) and Planckian occupation at T=H/2pi.
  NOTE: 4H is a decay rate, not a 'memory time'. The pointwise <T00 T00> correlator diverges as tau^-8 as tau -> 0,
  so no finite normalised autocorrelation time Int K dtau / K(0) exists without smearing the operator over a finite region.

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
== 3.1 horizon (membrane) fluid transport ==
  eta  = c^3/(16 pi G) = 8.0313e+33 kg s^-1  (= Pa s m: membrane-paradigm SURFACE shear viscosity of the horizon, a 2D fluid;
         not a bulk Pa s viscosity of empty space)
  s    = k_B c^3/(4 G hbar) = 1.3213e+46 J K^-1 m^-2  (surface entropy density)
  eta/s = 6.0783e-13 K s   vs hbar/(4 pi k_B) = 6.0783e-13 K s  -> PASS
== 3.2 thermal (beta/2pi) time scale tau_thermal = hbar/(2 pi k_B T) = 1/kappa ==
  T_dS(H0) = 2.655e-30 K ;  tau_thermal = hbar/(2 pi k_B T) = 4.5782e+17 s = 1/H0 = 4.5782e+17 s -> PASS
  BH M=1 Msun: T_H=6.170e-08 K, tau_thermal 1.9703e-05 s = 4GM/c^3 = 1.9703e-05 s
  BH M=10 Msun: T_H=6.170e-09 K, tau_thermal 1.9703e-04 s = 4GM/c^3 = 1.9703e-04 s
  BH M=1e+06 Msun: T_H=6.170e-14 K, tau_thermal 1.9703e+01 s = 4GM/c^3 = 1.9703e+01 s
  NOTE: tau_thermal = 1/kappa is the beta/2pi scale of the Planckian occupation. It is a DIFFERENT quantity from the
  1/(4H) of V1, which is the asymptotic (large-|t|) decay rate of the pointwise <T00(t) T00(0)> correlator. That pointwise
  correlator diverges as t^-8 as t -> 0, so no finite normalised autocorrelation time exists without smearing;
  neither number is a 'memory time' of the medium in the usual sense.
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
== 4.4 Im Pi(s) = c_eff(s) G s^2 : Standard-Model staircase in the timelike invariant mass sqrt(s) ==
  ILLUSTRATIVE ONLY: massive-vector (W, Z) threshold polynomial NOT computed (leading beta only);
  hadronic region 0.3-2 GeV schematic (3 real pions as proxy, perturbative quarks+gluons above 2 GeV).
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
  for a timelike graviton of invariant mass sqrt(s) ~ 4e-22 GeV: c_eff = 0.1000 (photon only; neutrinos closed by their masses)
  a free on-shell GW has s = 0, so no species is open regardless of frequency: the above is NOT a statement about LIGO-band GWs.
  on-shell flat-space absorption of a GW: Im Pi(s=0) = 0 exactly.  FRW order-of-magnitude estimate (s ~ H0^2 off-shellness,
  photon channel open) Gamma ~ c_eff t_P^2 H0^4/(2 omega) = 5.2e-162 s^-1 at omega = 2 pi x 100 Hz
  => the medium's absorptive response in sqrt(s) is derived (thresholds fixed by masses/spins), and unobservably weak for GWs: consistent with GW170817 trivially.

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
Status: numerical coincidence (Milgrom 1999, Verlinde 2016); no mechanism within GRUT III.
== dimensional-coincidence table ==
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
== discriminator: a0 ∝ H(z) (apparent horizon) vs a0 = const (Lambda-horizon) ==
  (this test has been proposed before; it is not a GRUT III prediction)
   z=0.5: H(z)/H0 = 1.322   (Lambda-horizon alternative: 1.000)
   z=1.0: H(z)/H0 = 1.790   (Lambda-horizon alternative: 1.000)
   z=2.0: H(z)/H0 = 3.032   (Lambda-horizon alternative: 1.000)
   z=3.0: H(z)/H0 = 4.566   (Lambda-horizon alternative: 1.000)
Status: numerical coincidence (Milgrom 1999, Verlinde 2016); no mechanism within GRUT III; the a0 ∝ H(z) test has been proposed before.

V6 OVERALL: INFORMATIONAL (no pass/fail)
```

## `v7_access_geometry_memory.py`

```text
== A. free scalar: smeared retarded kernel is fixed by the access geometry ==
-- A.1 convolution of Gaussians and radial integral against delta(t-r) (sympy) --
  [PASS] f_L * f_L = f_{sqrt2 L} (per Cartesian factor; 3D Gaussian factorises): (f_L*f_L)(d) = exp(-d**2/(4*L**2))/(2*sqrt(pi)*L)
  [PASS] chi_L(t) = int d^3r f_{sqrt2 L}(r) delta(t-r)/(4 pi r) = t exp(-t^2/4L^2)/(8 pi^{3/2} L^3): chi_L = t*exp(-t**2/(4*L**2))/(8*pi**(3/2)*L**3)
  [PASS] direct 3D quadrature (f_L*f_L)(1.3 L) vs f_{sqrt2 L}(1.3 L): got 0.0147128153797  want 0.0147128153797  rel.err 1.18e-16
-- A.2 normalisation, static limit, mean memory time (sympy) --
  [PASS] int_0^oo chi_L dt = 1/(4 pi^{3/2} L): = 1/(4*pi**(3/2)*L) = (1/(4 pi L^2)) * (L/sqrt(pi))
  [PASS] = smeared static Coulomb potential int d^3r f_{sqrt2 L}/(4 pi r) (independent integral): 1/(4*pi**(3/2)*L)
  [PASS] mean memory time <t> = sqrt(pi) L: <t> = sqrt(pi)*L
  [PASS] rho_L(w) = FT of <[F_L(t),F_L(0)]> = (w/2pi) exp(-w^2 L^2)  [= pointwise w/2pi x form factor]: rho_L = w*exp(-L**2*w**2)/(2*pi)
-- A.3 compact support: top hat of radius R, chi vanishes for t > 2R --
  [PASS] int d^3x (f_R*f_R) = 1 (lens-volume formula correctly normalised)
  [PASS] int_0^{2R} chi_R dt = int_0^{2R} t (f_R*f_R)(t) dt = 3/(10 pi R): = 3/(10*pi*R)
     d/R     numeric (f*f)(d)      analytic lens formula     chi_R(t=d) = d (f*f)(d)
      0.5       0.151072856138       0.151072856138       0.075536428069
  [PASS] top hat (f*f)(d=0.5R) vs lens volume: got 0.151072856138  want 0.151072856138  rel.err 0.0
      1.0      0.0746038795743      0.0746038795743      0.0746038795743
  [PASS] top hat (f*f)(d=1.0R) vs lens volume: got 0.0746038795743  want 0.0746038795743  rel.err 1.86e-16
      1.5      0.0205160668829      0.0205160668829      0.0307741003244
  [PASS] top hat (f*f)(d=1.5R) vs lens volume: got 0.0205160668829  want 0.0205160668829  rel.err 1.69e-16
      1.9    0.000880325778977    0.000880325778977     0.00167261898006
  [PASS] top hat (f*f)(d=1.9R) vs lens volume: got 0.000880325778977  want 0.000880325778977  rel.err 2.46e-16
     1.99       8.937544773e-6       8.937544773e-6     1.77857140983e-5
  [PASS] top hat (f*f)(d=1.99R) vs lens volume: got 8.937544773e-6  want 8.937544773e-6  rel.err 1.9e-16
      2.0                  0.0                  0.0                  0.0
  [PASS] top hat chi_R(t=2.0R) = 0 (t > 2R): got 0.0  want 0.0  abs.err 0.0
      2.1                  0.0                  0.0                  0.0
  [PASS] top hat chi_R(t=2.1R) = 0 (t > 2R): got 0.0  want 0.0  abs.err 0.0
      3.0                  0.0                  0.0                  0.0
  [PASS] top hat chi_R(t=3.0R) = 0 (t > 2R): got 0.0  want 0.0  abs.err 0.0

== B. conformal scalar stress tensor: smeared response kernel ==
-- B.0 conformal structure of <T00 T00> at a generic point (exact Wick contraction, xi=1/6) --
  [PASS] <T00(x)T00(0)>_Wick = C_T (I_00^2 - 1/4)/(x^2)^4, C_T = 1/(3 pi^4), I_00 = 1 - 2 tau^2/x^2, all x
  [PASS] time axis: I_00 = 1, I_{00,00} = 3/4, <T00 T00> = 1/(4 pi^4 t^8) (v1)
-- B.1 closed form of chi_L(t) two ways (sympy) --
  chi_L(t) = -(-120*L**6 + 180*L**4*t**2 - 30*L**2*t**4 + t**6)*exp(-t**2/(4*L**2))/(184320*pi**(7/2)*L**13*t)
  [PASS] position-space (distributional) closed form == momentum-space (Ward/CFT) closed form
  [PASS] small-t edge: t chi_L(t) -> 1/(1536 pi^{7/2} L^7) = rho_L(oo)/pi: series: 1/(1536*pi**(7/2)*L**7*t) - 7*t/(6144*pi**(7/2)*L**9)
  [PASS] rho_L(oo) = (1/(1440 pi^3)) int_0^oo K^6 e^{-K^2 L^2} dK = sqrt(pi)/(1536 pi^3 L^7)
-- B.2 quadrature on the rotated contour r = s e^{i pi/8} (eps -> 0 exactly) vs closed form --
  [PASS] contour rotation valid at finite eps=0.2L, t=0.8L: Re W (rotated vs real axis): got -4.27643215421e-6  want -4.27643215421e-6  rel.err 4.29e-29
  [PASS] contour rotation valid at finite eps=0.2L, t=0.8L: Im W (rotated vs real axis): got -1.87033777207e-6  want -1.87033777207e-6  rel.err 7.31e-29
     t/L      L^8 chi_L (quadrature)      L^8 chi_L (closed form)    L^8 Re W_L     chi_L(t)/chi_L(L)
     0.05      0.0002358912365475854      0.0002358912365475854   -1.067592e-5
  [PASS] chi_L(t=0.05 L) quadrature vs closed form: got 0.000235891236548  want 0.000235891236548  rel.err 1.14e-26
     0.25       4.232267927982421e-5       4.232267927982421e-5   -1.025506e-5
  [PASS] chi_L(t=0.25 L) quadrature vs closed form: got 4.23226792798e-5  want 4.23226792798e-5  rel.err 1.51e-29
      0.5        1.42556302650117e-5        1.42556302650117e-5   -9.0189526e-6
  [PASS] chi_L(t=0.5 L) quadrature vs closed form: got 1.4255630265e-5  want 1.4255630265e-5  rel.err 1.37e-30
     0.75       3.209385341224302e-6       3.209385341224302e-6   -7.2045739e-6
  [PASS] chi_L(t=0.75 L) quadrature vs closed form: got 3.20938534122e-6  want 3.20938534122e-6  rel.err 1.64e-30
        1      -2.383367742835294e-6      -2.383367742835294e-6   -5.1153059e-6
  [PASS] chi_L(t=1 L) quadrature vs closed form: got -2.38336774284e-6  want -2.38336774284e-6  rel.err 1.58e-31
      1.5      -5.419208795405391e-6      -5.419208795405391e-6   -1.3185309e-6
  [PASS] chi_L(t=1.5 L) quadrature vs closed form: got -5.41920879541e-6  want -5.41920879541e-6  rel.err 0.0
        2       -3.34115273622407e-6       -3.34115273622407e-6   7.5208761e-7
  [PASS] chi_L(t=2 L) quadrature vs closed form: got -3.34115273622e-6  want -3.34115273622e-6  rel.err 1.13e-31
        3       6.971322952178255e-7       6.971322952178255e-7    5.740814e-7
  [PASS] chi_L(t=3 L) quadrature vs closed form: got 6.97132295218e-7  want 6.97132295218e-7  rel.err 0.0
        4       3.724708384718247e-7       3.724708384718247e-7   -1.1563266e-7
  [PASS] chi_L(t=4 L) quadrature vs closed form: got 3.72470838472e-7  want 3.72470838472e-7  rel.err 1.26e-31
    chi_L(2L)/chi_L(L) = 1.40186;  Gaussian tail exp(-t^2/4L^2) x polynomial
    chi_L(3L)/chi_L(L) = -0.292499;  Gaussian tail exp(-t^2/4L^2) x polynomial
    chi_L(4L)/chi_L(L) = -0.156279;  Gaussian tail exp(-t^2/4L^2) x polynomial
  [PASS] edge: t chi_L(t) at t = 0.05 L vs rho_L(oo)/pi (1/t law, 1 - 7 t^2/(4 L^2) correction): got 1.17945618274e-5  want 1.17945132611e-5  rel.err 4.12e-6
-- B.3 scaling chi_L(t) = L^{-8} g(t/L) (numerical, L = 1/2, 1, 2 at t/L = 0.5 and 2) --
  [PASS] L^8 chi_L at t/L=0.5: L=1/2 vs L=1: got 1.4255630265e-5  want 1.4255630265e-5  rel.err 0.0
  [PASS] L^8 chi_L at t/L=0.5: L=2 vs L=1: got 1.4255630265e-5  want 1.4255630265e-5  rel.err 7.86e-28
  [PASS] L^8 chi_L at t/L=2: L=1/2 vs L=1: got -3.34115273622e-6  want -3.34115273622e-6  rel.err 0.0
  [PASS] L^8 chi_L at t/L=2: L=2 vs L=1: got -3.34115273622e-6  want -3.34115273622e-6  rel.err 0.0
  (F_L = int d^3x f_L T00 has mass dimension 4 => <F_L F_L> ~ L^-8; the task text's L^-9 is corrected to L^-8)
-- B.4 Fourier transform of the smeared Wightman function at T=0 (shifted line Im t = -L, numpy) --
      w L      W~_L(w) [L^7]        rho_L^{(b)}(w)        |W~_L(-w)|/W~_L(w)    W~_L/(w^7/10080 pi^3)
     0.25   1.8603091588e-10   1.8603091588e-10   2.01e-12            0.95261017
  [PASS] W~_L(w=0.25/L) = rho_L^{(b)}(w) (Ward/CFT form): got 1.86030915881e-10  want 1.86030915879e-10  rel.err 1.04e-11
  [PASS] positive frequency: |W~_L(-0.25/L)|/W~_L(0.25/L) < 1e-9: 2.01e-12
      0.5   2.0599957730e-08   2.059995773e-8   2.64e-14            0.8241124
  [PASS] W~_L(w=0.5/L) = rho_L^{(b)}(w) (Ward/CFT form): got 2.05999577304e-8  want 2.05999577304e-8  rel.err 1.49e-13
  [PASS] positive frequency: |W~_L(-0.5/L)|/W~_L(0.5/L) < 1e-9: 2.64e-14
      1.0   1.4945960123e-06   1.4945960123e-6   1.33e-17            0.46712592
  [PASS] W~_L(w=1.0/L) = rho_L^{(b)}(w) (Ward/CFT form): got 1.49459601234e-6  want 1.49459601234e-6  rel.err 3.02e-13
  [PASS] positive frequency: |W~_L(-1.0/L)|/W~_L(1.0/L) < 1e-9: 1.33e-17
      2.0   2.4838437512e-05   2.4838437512e-5   1.24e-18            0.060649113
  [PASS] W~_L(w=2.0/L) = rho_L^{(b)}(w) (Ward/CFT form): got 2.48384375116e-5  want 2.48384375116e-5  rel.err 2.55e-13
  [PASS] positive frequency: |W~_L(-2.0/L)|/W~_L(2.0/L) < 1e-9: 1.24e-18
      3.0   3.6770897736e-05   3.6770897736e-5   1.69e-18            0.0052549138
  [PASS] W~_L(w=3.0/L) = rho_L^{(b)}(w) (Ward/CFT form): got 3.67708977363e-5  want 3.67708977363e-5  rel.err 5.22e-14
  [PASS] positive frequency: |W~_L(-3.0/L)|/W~_L(3.0/L) < 1e-9: 1.69e-18
      4.0   3.7214866290e-05   3.721486629e-5   3.82e-18            0.00070991553
  [PASS] W~_L(w=4.0/L) = rho_L^{(b)}(w) (Ward/CFT form): got 3.72148662901e-5  want 3.72148662901e-5  rel.err 3.18e-13
  [PASS] positive frequency: |W~_L(-4.0/L)|/W~_L(4.0/L) < 1e-9: 3.82e-18
      6.0   3.7216377867e-05   3.7216377867e-5   4.45e-18            4.1551385e-5
  [PASS] W~_L(w=6.0/L) = rho_L^{(b)}(w) (Ward/CFT form): got 3.72163778669e-5  want 3.7216377867e-5  rel.err 1.38e-12
  [PASS] positive frequency: |W~_L(-6.0/L)|/W~_L(6.0/L) < 1e-9: 4.45e-18
  [PASS] small-w limit W~_L(0.25/L)/(w^7/10080 pi^3) -> 1 - (7/9)(wL)^2 (pointwise v1 spectrum recovered): got 0.952610170026  want 0.951388888889  rel.err 0.00128
  [PASS] large-w plateau W~_L(6/L) -> sqrt(pi)/(1536 pi^3 L^7): got 3.72163778669e-5  want 3.7216377867e-5  rel.err 1.96e-12

== C. KMS for smeared observables: state x commutator = noise, for any smearing ==
-- C.1 image sum closed form --
  [PASS] sum_n 1/(4pi^2 (r^2-(t+i n beta)^2)) = [coth(pi(r-t)/beta)+coth(pi(r+t)/beta)]/(8 pi r beta) at a complex point: got 6.3526500706e-33  want 0.0  abs.err 6.35e-33
-- C.2 smeared thermal Wightman spectrum, shifted line Im t = -beta/4 (generic: not the KMS-symmetric beta/2) --
  beta = 6.283185, L = 1.0:
      w        G+_L(w)            G+_L(-w)          G+/G- e^-bw - 1    rho_L/rho_A - 1    N_L/(coth rho/2) - 1
     0.25   4.7187343013e-02   9.809285e-03   +1.33e-15          +2.29e-13          +2.29e-13
  [PASS] KMS  G+_L(0.25)/G+_L(-0.25) = e^(beta w)  [beta=6.283, L=1.0]: 1.33e-15
  [PASS] rho_L(0.25) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(0.25) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(0.25) = 0: 6.4e-19
      0.5   6.4774141616e-02   2.799144e-03   +2.66e-15          +2.29e-13          +2.28e-13
  [PASS] KMS  G+_L(0.5)/G+_L(-0.5) = e^(beta w)  [beta=6.283, L=1.0]: 2.66e-15
  [PASS] rho_L(0.5) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(0.5) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.28e-13
  [PASS] Im G+_L(0.5) = 0: 0.0e+00
      1.0   5.8659374547e-02   1.095430e-04   +1.49e-14          +2.28e-13          +2.28e-13
  [PASS] KMS  G+_L(1.0)/G+_L(-1.0) = e^(beta w)  [beta=6.283, L=1.0]: 1.49e-14
  [PASS] rho_L(1.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.28e-13
  [PASS] N_L(1.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.28e-13
  [PASS] Im G+_L(1.0) = 0: 0.0e+00
      2.0   5.8300692615e-03   2.033145e-08   +5.14e-11          +2.11e-13          +2.11e-13
  [PASS] KMS  G+_L(2.0)/G+_L(-2.0) = e^(beta w)  [beta=6.283, L=1.0]: 5.14e-11
  [PASS] rho_L(2.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.11e-13
  [PASS] N_L(2.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.11e-13
  [PASS] Im G+_L(2.0) = 0: 2.0e-17
  beta = 6.283185, L = 0.5:
      w        G+_L(w)            G+_L(-w)          G+/G- e^-bw - 1    rho_L/rho_A - 1    N_L/(coth rho/2) - 1
     0.25   4.9451910886e-02   1.028004e-02   +1.11e-15          +2.29e-13          +2.28e-13
  [PASS] KMS  G+_L(0.25)/G+_L(-0.25) = e^(beta w)  [beta=6.283, L=0.5]: 1.11e-15
  [PASS] rho_L(0.25) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(0.25) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.28e-13
  [PASS] Im G+_L(0.25) = 0: 2.6e-18
      0.5   7.8132528998e-02   3.376413e-03   +2.66e-15          +2.29e-13          +2.28e-13
  [PASS] KMS  G+_L(0.5)/G+_L(-0.5) = e^(beta w)  [beta=6.283, L=0.5]: 2.66e-15
  [PASS] rho_L(0.5) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(0.5) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.28e-13
  [PASS] Im G+_L(0.5) = 0: 9.5e-19
      1.0   1.2418189689e-01   2.319026e-04   +1.33e-14          +2.28e-13          +2.28e-13
  [PASS] KMS  G+_L(1.0)/G+_L(-1.0) = e^(beta w)  [beta=6.283, L=0.5]: 1.33e-14
  [PASS] rho_L(1.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.28e-13
  [PASS] N_L(1.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.28e-13
  [PASS] Im G+_L(1.0) = 0: 4.2e-18
      2.0   1.1710007142e-01   4.083680e-07   +2.66e-12          +2.27e-13          +2.27e-13
  [PASS] KMS  G+_L(2.0)/G+_L(-2.0) = e^(beta w)  [beta=6.283, L=0.5]: 2.66e-12
  [PASS] rho_L(2.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.27e-13
  [PASS] N_L(2.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.27e-13
  [PASS] Im G+_L(2.0) = 0: 2.0e-17
  beta = 1.000000, L = 0.3:
      w        G+_L(w)            G+_L(-w)          G+/G- e^-bw - 1    rho_L/rho_A - 1    N_L/(coth rho/2) - 1
     0.25   1.7886841091e-01   1.393029e-01   +0.00e+00          +2.28e-13          +2.29e-13
  [PASS] KMS  G+_L(0.25)/G+_L(-0.25) = e^(beta w)  [beta=1, L=0.3]: 0.00e+00
  [PASS] rho_L(0.25) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.28e-13
  [PASS] N_L(0.25) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(0.25) = 0: 3.7e-18
      0.5   1.9774595703e-01   1.199390e-01   +4.44e-16          +2.30e-13          +2.29e-13
  [PASS] KMS  G+_L(0.5)/G+_L(-0.5) = e^(beta w)  [beta=1, L=0.3]: 4.44e-16
  [PASS] rho_L(0.5) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.30e-13
  [PASS] N_L(0.5) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(0.5) = 0: 0.0e+00
      1.0   2.3010905713e-01   8.465239e-02   +2.22e-16          +2.29e-13          +2.29e-13
  [PASS] KMS  G+_L(1.0)/G+_L(-1.0) = e^(beta w)  [beta=1, L=0.3]: 2.22e-16
  [PASS] rho_L(1.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(1.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(1.0) = 0: 2.8e-19
      2.0   2.5683628306e-01   3.475901e-02   +2.00e-15          +2.29e-13          +2.29e-13
  [PASS] KMS  G+_L(2.0)/G+_L(-2.0) = e^(beta w)  [beta=1, L=0.3]: 2.00e-15
  [PASS] rho_L(2.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(2.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(2.0) = 0: 0.0e+00
      3.0   2.2353314038e-01   1.112906e-02   +4.66e-15          +2.29e-13          +2.29e-13
  [PASS] KMS  G+_L(3.0)/G+_L(-3.0) = e^(beta w)  [beta=1, L=0.3]: 4.66e-15
  [PASS] rho_L(3.0) = G+ - G- = (w/2pi) e^(-w^2 L^2) (commutator from A): 2.29e-13
  [PASS] N_L(3.0) = (G+ + G-)/2 = (1/2) coth(beta w/2) rho_L: 2.29e-13
  [PASS] Im G+_L(3.0) = 0: 0.0e+00

  (runtime 26.5 s)
V7 OVERALL: PASS
```

## `v8_residual_framework.py`

```text
================================================================================================
V8  residual framework:  Delta rho = rho_GRUT - rho_conventional is the decisive object   (H = 1, beta = 2 pi/H)
================================================================================================

== (A) helpers, tested on exact objects ==
  [PASS] Planckian G+: rho(w) = w/2pi exactly at w in {0.1, 1, 3}  (rho(1) = 0.159154943091895 vs 1/2pi = 0.159154943091895)
  [PASS] Planckian G+: N(w) = (w/4pi) coth(beta w/2)
  [PASS] KMS defect of the Planckian spectrum = 0  (max defect 9.86e-32)
  [PASS] mutation: G+ at temperature T/2 is flagged (defect > 1e-3)  (defect at w=0.5: 0.00344529  -- tester live)
  [PASS] residual of a model against itself vanishes
  [PASS] shape residual of w^3 against w^7 = -4  (got -4.0)

== (B) FDT cancellation for G+_T (v1 closed form): rho_T is an odd polynomial ==
  [PASS] sympy: G+_T(w) - G+_T(-w) == P7(w)/(10080 pi^3)  (P7 = 36*H**6*omega + 49*H**4*omega**3 + 14*H**2*omega**5 + omega**7)
  [PASS] P7 coefficients (1, 0, 14H^2, 0, 49H^4, 0, 36H^6, 0): odd polynomial of degree 7  (coeffs = [1, 0, 14*H**2, 0, 49*H**4, 0, 36*H**6, 0])
  [PASS] thermal denominators cancel: 1/(1-e^-x) + 1/(1-e^x) = 1
  [PASS] sympy: N_T(w) == (1/2) coth(pi w/H) rho_T(w)
  [PASS] N_T(0) = G+_T(0) = H^7/(560 pi^4)  (classical Nyquist floor)  (limit = H**7/(560*pi**4))
  [PASS] numeric: rho_from_Gplus(G+_T) == P7/(10080 pi^3) at the table frequencies
  [PASS] numeric: KMS defect of G+_T = 0  (max defect 2.41e-35)
  [PASS] rho_T local exponent: 1 at w = 1e-3 H (IR, from 36 H^6 w) and 7 at w = 1e3 H (flat limit)  (slope(1e-3) = 1.000002722, slope(1e3) = 6.999972)
  -> rho_T(t) = sum_k c_k delta^(k)(t): the dissipative kernel at a point is LOCAL (Huygens); all memory sits in
     N_T = (1/2) coth(pi w/H) rho_T, which is fixed by rho_T and the temperature T = H/2pi -- nothing is left to supply.

== (C) old GRUT kernels against the conformal medium  (normalised to agree with rho_T at w = 1) ==
  Conformal medium reference: rho_T(w) = (w^7 + 14 w^5 + 49 w^3 + 36 w)/(10080 pi^3), slope 1 (w<<H) -> 7 (w>>H).

  --- C1  rung-3 noise model  rho = J = w^3 exp(-(w/omega_c)^2), omega_c = 1  [calc/finite_T_exponent.py] ---   normalisation constant c = 0.00086972976 (the record left this amplitude free)
       w |  rho_old/rho_T | Delta rho / rho_T |  slope_old |  slope_T | shape resid
     0.1 |    0.073749826 |       -0.92625017 |       2.98 |  1.02701 |     1.95299
     0.3 |     0.55174296 |       -0.44825704 |       2.82 |  1.22895 |     1.59105
     1.0 |            1.0 |               0.0 |        1.0 |      2.6 |        -1.6
     3.0 |  0.00012902409 |       -0.99987098 |      -15.0 |  5.18462 |    -20.1846
    10.0 |  8.8321224e-46 |              -1.0 |     -197.0 |  6.73814 |    -203.738
  IR exponent (w = 1e-3): old = 2.999998  vs conformal medium = 1.0000027;  UV exponent of the medium (w = 1e3) = 6.999972
  [PASS] C1  rung-3 noise model  rho = J = w^3 exp(-(w/omega_c)^2), omega_c = 1  [calc/finite_T_exponent.py]: IR exponent = 3  (got 2.999998)
  [PASS] C1  rung-3 noise model  rho = J = w^3 exp(-(w/omega_c)^2), omega_c = 1  [calc/finite_T_exponent.py]: Delta rho != 0 -- |Delta rho|/rho_T >= 10% at every tabulated w != 1  (min |Delta rho|/rho_T off w=1: 0.448257)
  [PASS] C1  rung-3 noise model  rho = J = w^3 exp(-(w/omega_c)^2), omega_c = 1  [calc/finite_T_exponent.py]: shape residual |slope_old - slope_T| >= 0.5 at w = 0.1 and 10.0  (shape residual: at 0.1: 1.95299, at 10.0: -203.738)
  (the record's S(w,T) = J coth(w/2T) is 2 N(w) for rho = J: KMS-locked by construction, so N adds no information)
  UV: the Gaussian cutoff gives slope 3 - 2 (w/omega_c)^2 = -15.0 already at w = 3: no power law, vs 7 for the medium; the cutoff scale omega_c/H is a second free choice

  --- C2  S5-1 Lorentzian, odd extension  rho = sgn(w) (w^2+kappa^2)^-1, kappa = 1/(2 sqrt 2.3)  [P06/P06B, p = 1] ---   normalisation constant c = 0.00035473349 (the record left this amplitude free)
       w |  rho_old/rho_T | Delta rho / rho_T |  slope_old |  slope_T | shape resid
     0.1 |      255.96878 |         254.96878 |  -0.168498 |  1.02701 |    -1.19551
     0.3 |      45.897501 |         44.897501 |  -0.905908 |  1.22895 |    -2.13486
     1.0 |            1.0 |               0.0 |   -1.80392 |      2.6 |    -4.40392
     3.0 |     0.00173388 |       -0.99826612 |   -1.97613 |  5.18462 |    -7.16075
    10.0 |   9.6729586e-8 |        -0.9999999 |   -1.99783 |  6.73814 |    -8.73597
  IR exponent (w = 1e-3): old = -1.8399831e-5  vs conformal medium = 1.0000027;  UV exponent of the medium (w = 1e3) = 6.999972
  [PASS] C2  S5-1 Lorentzian, odd extension  rho = sgn(w) (w^2+kappa^2)^-1, kappa = 1/(2 sqrt 2.3)  [P06/P06B, p = 1]: IR exponent = 0  (got -1.8399831e-5)
  [PASS] C2  S5-1 Lorentzian, odd extension  rho = sgn(w) (w^2+kappa^2)^-1, kappa = 1/(2 sqrt 2.3)  [P06/P06B, p = 1]: Delta rho != 0 -- |Delta rho|/rho_T >= 10% at every tabulated w != 1  (min |Delta rho|/rho_T off w=1: 0.998266)
  [PASS] C2  S5-1 Lorentzian, odd extension  rho = sgn(w) (w^2+kappa^2)^-1, kappa = 1/(2 sqrt 2.3)  [P06/P06B, p = 1]: shape residual |slope_old - slope_T| >= 0.5 at w = 0.1 and 10.0  (shape residual: at 0.1: -1.19551, at 10.0: -8.73597)
  UV exponent of (w^2+kappa^2)^-p is -2p = -2 (p = 1): the family value;  kappa = 0.32969024
  family member p = 0.5: IR exponent -9.19992e-6 (0 for every p), UV exponent -1.0 (= -2p)
  family member p = 2.0: IR exponent -3.67997e-5 (0 for every p), UV exponent -4.0 (= -2p)

  --- C2' S5-1 Lorentzian, Drude-friction reading  rho = w (w^2+kappa^2)^-1  (Ohmic-Drude) ---   normalisation constant c = 0.00035473349 (the record left this amplitude free)
       w |  rho_old/rho_T | Delta rho / rho_T |  slope_old |  slope_T | shape resid
     0.1 |      25.596878 |         24.596878 |   0.831502 |  1.02701 |   -0.195507
     0.3 |       13.76925 |          12.76925 |  0.0940919 |  1.22895 |    -1.13486
     1.0 |            1.0 |               0.0 |  -0.803922 |      2.6 |    -3.40392
     3.0 |     0.00520164 |       -0.99479836 |  -0.976134 |  5.18462 |    -6.16075
    10.0 |   9.6729586e-7 |       -0.99999903 |  -0.997828 |  6.73814 |    -7.73597
  IR exponent (w = 1e-3): old = 0.9999816  vs conformal medium = 1.0000027;  UV exponent of the medium (w = 1e3) = 6.999972
  [PASS] C2' S5-1 Lorentzian, Drude-friction reading  rho = w (w^2+kappa^2)^-1  (Ohmic-Drude): IR exponent = 1  (got 0.9999816)
  [PASS] C2' S5-1 Lorentzian, Drude-friction reading  rho = w (w^2+kappa^2)^-1  (Ohmic-Drude): Delta rho != 0 -- |Delta rho|/rho_T >= 10% at every tabulated w != 1  (min |Delta rho|/rho_T off w=1: 0.994798)
  [PASS] C2' S5-1 Lorentzian, Drude-friction reading  rho = w (w^2+kappa^2)^-1  (Ohmic-Drude): shape residual |slope_old - slope_T| >= 0.5 at w = 1.0 and 10.0  (shape residual: at 1.0: -3.40392, at 10.0: -7.73597)
  NOTE: in the Drude reading the IR exponent 1 coincides with the medium's IR exponent 1 (w << H), so the shape residual
        at w = 0.1 is only -0.1955 (tested at w = 1 and 10 instead); the amplitude residual there is still
        24.6 x rho_T, and for w >= 1 Drude falls as 1/w while the medium rises as w^7: Delta rho != 0 regardless.

  CONCLUSION (C): every registered GRUT kernel asserts Delta rho != 0 against the conformal vacuum medium of the
  same (observable, state), with exponents 3 / 0 (or 1) in place of the medium's 1 -> 7, and gives no operator or state
  that would produce the difference.  That is exactly what SCOUT_0 recorded as 'supplied' (noise profile, kernel class).

== (D) residual_report template, demonstrated ==
  residual_report: observable = T^00 at a point on a comoving worldline (one conformal scalar); state = Bunch-Davies / geodesic observer, KMS at T = H/2pi
    max |Delta rho|/rho_ref = 0; max |shape residual| = 0; IR exponent model/ref = 1.0000/1.0000; KMS defect of implied N = 0.0e+00
    VERDICT: EMPTY (Delta rho = 0)
  [PASS] template on the conformal medium itself returns EMPTY (Delta rho = 0)
  residual_report: observable = unspecified ('fast modes', rung 3); state = unspecified thermal state, T free
    max |Delta rho|/rho_ref = 1; max |shape residual| = 203.7; IR exponent model/ref = 3.0000/1.0000; KMS defect of implied N = 0.0e+00
    VERDICT: SUPPLIED (Delta rho != 0, no mechanism named)
  [PASS] template on C1 returns SUPPLIED  (IR exponent 3.0000 vs ref 1.0000)
  residual_report: observable = unspecified (S5-1 scaled-limit resolvent density); state = unspecified, kappa = 1/(2 sqrt 2.3)
    max |Delta rho|/rho_ref = 255; max |shape residual| = 8.736; IR exponent model/ref = -0.0000/1.0000; KMS defect of implied N = 0.0e+00
    VERDICT: SUPPLIED (Delta rho != 0, no mechanism named)
  [PASS] template on C2 returns SUPPLIED  (IR exponent -0.0000 vs ref 1.0000)
  residual_report: observable = unspecified; state = unspecified
    max |Delta rho|/rho_ref = 1; max |shape residual| = 203.7; IR exponent model/ref = 3.0000/1.0000; KMS defect of implied N = 0.0e+00
    VERDICT: TESTABLE (Delta rho != 0; mechanism: (hypothetical) named operator content)
  [PASS] template on C1 with a named mechanism returns TESTABLE
  [PASS] template: implied noise of every model is KMS-consistent by construction (defect < 1e-25)

V8 OVERALL: PASS
```

## `v9_induced_gravity_matching.py`

```text
== A1. Flat-plane entanglement entropy from the cone heat kernel (Solodukhin eqs. 20, 26-28), d = 4 ==
  W_area(alpha) = A*(alpha**2 - 1)/(96*pi*alpha*epsilon**2)
  S = (alpha d/dalpha - 1) W |_1 = A/(48*pi*epsilon**2)   ->   S/A = 1/(48*pi) / eps^2
  [PASS] S_scalar/A = 1/(48 pi eps^2)  [Solodukhin eq. 28 with d = 4]
  [PASS] eq. (28) general-d formula 1/(6(d-2)(4 pi)^((d-2)/2)) reduces to 1/(48 pi) at d = 4
  <phi^2>_eps = int_(eps^2)^inf ds (4 pi s)^-2 = 1/(16*pi**2*epsilon**2)
  [PASS] <phi^2>_eps = 1/(16 pi^2 eps^2)

== A2. One-loop delta(1/G) from a_1 per species (Solodukhin eqs. 112-115, 131-148; Kabat c_1 table) ==
  W_div (int R part) = -IntR*k_1/(32*pi**2*epsilon**2)   =>   delta(1/G) = k_1/(2*pi*epsilon**2)
  [PASS] delta(1/G) = k_1/(2 pi eps^2)
  [PASS] minimal scalar: 1/G_ren = 1/G_B + 1/(12 pi eps^2)  [Solodukhin eq. 115]
  species                 k_1            delta(1/G)         delta(1/(4G))               S_ent/A              mismatch
  real scalar xi=0        1/6       1/(12 pi eps^2)       1/(48 pi eps^2)       1/(48 pi eps^2)                     0
  real scalar xi=1/6        0                     0                     0       1/(48 pi eps^2)      -1/(48 pi eps^2)
  Dirac fermion           1/3        1/(6 pi eps^2)       1/(24 pi eps^2)       1/(24 pi eps^2)                     0
  Weyl fermion            1/6       1/(12 pi eps^2)       1/(48 pi eps^2)       1/(48 pi eps^2)                     0
  Maxwell + ghosts       -2/3       -1/(3 pi eps^2)      -1/(12 pi eps^2)       1/(24 pi eps^2)       -1/(8 pi eps^2)
  [PASS] k_1(scalar, xi) = 1/6 - xi  [Solodukhin eq. 141; Kabat c_1 = 1/6 at xi = 0]
  [PASS] k_1(Dirac) = +1/3  [Solodukhin eq. 136 with D_1/2 = 2 and fermionic sign; Kabat c_1 = 2^[d/2]/12]
  [PASS] k_1(Maxwell+ghosts) = D_1/6 - 1 = -2/3  [Solodukhin eq. 147; Kabat c_1 = (d-2)/6 - 1]
  [PASS] Dirac delta(1/(4G)) = D_1/2/(48 pi eps^2) = 1/(24 pi eps^2)  [Solodukhin eq. 137]
  [PASS] Maxwell delta(1/(4G)) = (D_1/6 - 1)/(8 pi eps^2) = -1/(12 pi eps^2)  [Solodukhin eq. 148]
  [PASS] Kabat eq. (1.1): 2 pi c_1 A int ds/(4 pi s)^2 = c_1 A/(8 pi eps^2) = delta(1/(4G)) A
  [PASS] 2D: Larsen-Wilczek c(fermion)/c(boson) = 1/2 per component -> Dirac(2 comp.) = scalar = D_1/2(d=2) = 1

== B. Matching and mismatch, species by species; contact (Wald) terms ==
  [PASS] minimal scalar: delta(1/(4G)) = S_ent/A = 1/(48 pi eps^2)
  [PASS] Dirac fermion:  delta(1/(4G)) = S_ent/A = 1/(24 pi eps^2)
  [PASS] Weyl fermion:   delta(1/(4G)) = S_ent/A = 1/(48 pi eps^2)
  [PASS] conformal scalar xi=1/6: delta(1/(4G)) = 0 while S_ent/A = 1/(48 pi eps^2): MISMATCH = -1/(48 pi eps^2)
  [PASS] generic xi: mismatch = -xi/(8 pi eps^2)  [Solodukhin eq. 142]
  Wald contact term  -2 pi xi <phi^2>_eps = -xi/(8*pi*epsilon**2)   (xi = 1/6: -1/(48*pi*epsilon**2))
  [PASS] -2 pi xi <phi^2>_eps = -xi/(8 pi eps^2) = the mismatch (eq. 142)
  [PASS] S_ent + S_contact = (1/6 - xi)/(8 pi eps^2) = delta(1/(4G))  [Solodukhin eq. 144 = eq. 141]
  [PASS] xi = 1/6: S_ent + S_contact = 1/(48 pi) - 1/(48 pi) = 0 = delta(1/(4G))
  [PASS] replica trick on eq. (143) surface term gives the same -xi A/(8 pi eps^2)
  Kabat contact entropy (beta d/dbeta - 1)[A (2 pi - beta) <phi^2>_eps]|_(2 pi) = -A/(8*pi*epsilon**2)
  [PASS] Maxwell contact term = -A/(8 pi eps^2) = -2 pi <phi^2>_eps A  (Tr X^(1)/R = 1: 'xi_eff = 1')
  [PASS] Maxwell: S_ent/A + S_contact/A = 1/(24 pi) - 1/(8 pi) = -1/(12 pi eps^2) = delta(1/(4G))  [Kabat sec. 1; Solodukhin eq. 148]
  [PASS] bosonic rule  delta(1/(4G)) - S_ent/A = -(Tr X/R)/(8 pi eps^2)  for real scalar xi=0
  [PASS] bosonic rule  delta(1/(4G)) - S_ent/A = -(Tr X/R)/(8 pi eps^2)  for real scalar xi=1/6
  [PASS] bosonic rule  delta(1/(4G)) - S_ent/A = -(Tr X/R)/(8 pi eps^2)  for Maxwell + ghosts
  (Fermions: Kabat shows no contact term for spin 1/2; the X = R/4 term is absorbed by the statistics sign, eq. 136.)

== C. Species sums, Standard Model content; Sakharov 1/G_ind = (sum k_1) Lambda_UV^2/(2 pi), Lambda_UV = 1/eps ==
  Minimal-scalar units: N_eff := 6 sum k_1, so 1/G_ind = N_eff Lambda^2/(12 pi) and Lambda = sqrt(12 pi/N_eff) M_P.
  SM (N_s=4, N_W=45, N_V=12)   xi = 0   :  sum k_1 =   1/6  -> 1/G_ind = (1/12) Lambda^2/pi ;  N_eff = 1 ;  Lambda_UV = sqrt(12 pi/1) M_P = 6.140 M_P = 7.50e+19 GeV
  SM (N_s=4, N_W=45, N_V=12)   xi = 1/6 :  sum k_1 =  -1/2  -> 1/G_ind = (-1/4) Lambda^2/pi ;  N_eff = -3 ;  1/G_ind < 0: SIGN PROBLEM, no real Lambda_UV (induced G would be negative)
  SM (N_s=4, N_W=45, N_V=12)   entropy side: sum S_ent/A = 73/(48 pi eps^2)  (N_eff^(S) = 73, all species positive)  -> if matched alone: Lambda = sqrt(12 pi/73) M_P = 0.719 M_P = 8.77e+18 GeV
  SM + 3 nu_R (N_W=48)         xi = 0   :  sum k_1 =   2/3  -> 1/G_ind =  (1/3) Lambda^2/pi ;  N_eff = 4 ;  Lambda_UV = sqrt(12 pi/4) M_P = 3.070 M_P = 3.75e+19 GeV
  SM + 3 nu_R (N_W=48)         xi = 1/6 :  sum k_1 =     0  -> 1/G_ind =                  0 ;  N_eff = 0 ;  1/G_ind = 0 at one loop: NO induced G (Lambda undefined)
  SM + 3 nu_R (N_W=48)         entropy side: sum S_ent/A = 19/(12 pi eps^2)  (N_eff^(S) = 76, all species positive)  -> if matched alone: Lambda = sqrt(12 pi/76) M_P = 0.704 M_P = 8.60e+18 GeV
  [PASS] SM, xi = 0:   sum k_1 = (4 + 45 - 48)/6 = +1/6  -> 1/G_ind = Lambda^2/(12 pi), N_eff = 1
  [PASS] SM, xi = 0:   Lambda_UV = sqrt(12 pi) M_P = 6.14 M_P  (above M_P: one-loop dominance fails)
  [PASS] SM, xi = 1/6: sum k_1 = (0 + 45 - 48)/6 = -1/2 < 0  -> SIGN PROBLEM
  [PASS] SM entropy side: sum S_ent/A = 73/(48 pi eps^2) > 0
  [PASS] SM + 3 nu_R: xi = 0 gives N_eff = 4; xi = 1/6 gives sum k_1 = 0 exactly (vectors cancel 48 Weyl)
  Reading: with the contact terms included (which is what renormalises G), the SM's 12 vectors (k_1 = -2/3 each, total -8)
  nearly cancel its 45 Weyl fermions (+7.5) and 4 scalars (+2/3 at xi=0, 0 at xi=1/6). The induced 1/G is tiny or negative.
  Sakharov's programme at one loop with SM content does NOT give G: either a bare 1/G_B, more species, or a UV completion.
  (Visser gr-qc/0204062 secs. 3-5 reaches the same conclusion: one-loop induced gravity needs BSM content.)

== D. Regulator dependence: the same scalar area law in a different cutoff ==
  heat-kernel: S/A = 1/(48 pi eps^2);  Susskind-Uglum eq. (2.24): S/A = 1/(360 pi eps_SU^2);  ratio = 15/2 = (eps_SU/eps)^2 if identified
  [PASS] both are area laws c A/eps^2 with c > 0; coefficients differ (48 vs 360) -> the COEFFICIENT is regulator-dependent
  [PASS] Cooperman-Luty eq. (1.5): S_ent = -4 pi c_2 Lambda^2 A with EH coefficient -M_P^2/(16 pi) -> S_ent = (1/4) delta(M_P^2) A  [scheme-free statement]

V9 OVERALL: PASS
Verdict for O1: species-by-species matching delta(1/(4G)) = (S_ent + S_contact)/A holds in the heat-kernel scheme
(theorem for free matter loops); the coefficients are scheme-dependent; the finite G is NOT determined. P2 is consistent, not derived.
```

## `v10_graviton_sector_memory.py`

```text
== 10.1 worldline invariant ==
  (1-Z) - 2 sinh^2(H dt/2) = 0 -> PASS
  flat limit coefficient of 1/dt^2: 1/(4*pi**2) -> PASS
== 10.2 stationarity of the velocity correlator ==
  d_t d_t' G depends on tbar? False -> PASS
== 10.3 Brownian covariance at large separation ==
  G - (H^2/8pi^2)(2H tbar - H dt) ->  H**2*log(2)/(8*pi**2)  as dt->oo (a constant: PASS)
  => G ~ D * min(t,t') + const, D = H**3/(4*pi**2) -> PASS
== 10.4 secular growth of the renormalised variance ==
  d<phi^2>/dt = H**3/(4*pi**2) -> PASS
== 10.5 tensor translation ==
  d<h_ij h^ij>/dt = 2*H**3/(pi**2*M_P**2) ;  H*P_T = 2*H**3/(pi**2*M_P**2) -> PASS
== 10.6 KMS periodicity of the stationary part ==
  stationary part shifted by 2 pi i/H minus itself: 0 -> PASS

V10 OVERALL: PASS
```

## `v11_quasi_ds_kms_defect.py`

```text
== A. power law: exact first-order thermality and the second-order coefficient (sympy, explicit expansion) ==
  O(1): 0   O(1/q): 0 -> PASS (no first-order term)
  O(1/q^2): x**2*(x/tanh(x/2) - 3)/12  vs -x^2/4 + (x^3/12)coth(x/2): PASS
  [PASS]   exact ln R * q^2 at q=1e5, x=0.7: got -0.0375255308 want -0.0375255308 rel 7.61e-12
  [PASS]   exact ln R * q^2 at q=1e5, x=2.0: got -0.124643143 want -0.124643143 rel 4.43e-11
  => W_eps/W_dS(H_eff) = 1 + eps^2 [x^2/4 - (x^3/12) coth(x/2)] + O(eps^3), H_eff = Hbar(1-eps).  Reviewers confirmed.
  f2(x - 2 pi i) - f2(-x) = -(-x**3*(exp(x) + 1) + (x - 2*I*pi)**3*(exp(x) + 1) + 3*(x**2 - (x - 2*I*pi)**2)*(exp(x) - 1))/(12*exp(x) - 12)  (nonzero => KMS broken at O(eps^2))
== B. general a(t) to second adiabatic order (sympy construction + numerical verification) ==
  max |first-order difference from thermal form at H_eff = H(1-eps)| over random points: 2.17e-15 (float inputs; zero to machine precision) -> PASS
  max |second-order (eta=0) - power-law coefficient| over random points: 9.02e-16 -> PASS
  eps*eta term of the second-order relative remainder (general a(t)): epsilon*(x**2*(x + 6) - x**2*(x*exp(x) - 6*exp(x) + 12)*exp(x) - 24*x*exp(2*x) + 24*x + 48*(exp(x) - 2)*exp(x) + 48)*exp(-x)/(48*(cosh(x) - 1))
  THEOREM (any a(t)): to first adiabatic order the comoving conformal vacuum is thermal at T_eff = (H + Hdot/H)/2pi = (d/dt ln a_dot)/2pi.
== C. numerics: Richardson-extrapolated first- and second-order KMS defects (power law) ==
  D/eps relative to Hbar, Richardson over p = [200, 400, 800, 1600]  (expected 2 pi w/Hbar, i.e. c = 1 exactly):
  [PASS]    lim D/eps at w/H=0.5 = 2 pi w/H: got 3.141545863 want 3.141592654 rel 1.49e-5
  [PASS]    lim D/eps at w/H=1.0 = 2 pi w/H: got 6.283026938 want 6.283185307 rel 2.52e-5
  [PASS]    lim D/eps at w/H=2.0 = 2 pi w/H: got 12.56598506 want 12.56637061 rel 3.07e-5
  [PASS]    lim D/eps at w/H=3.0 = 2 pi w/H: got 18.84894548 want 18.84955592 rel 3.24e-5
  second-order defect D2(w) = lim [ln(G+/G-) - 2 pi w/H_eff]/eps^2  (the genuine non-thermal content):
    w/H = 0.5:  D2 = 2.501512   (p=1600 raw 2.508621)   D2/(2 pi w/H) = 0.796256
    w/H = 1.0:  D2 = 10.99629   (p=1600 raw 11.02289)   D2/(2 pi w/H) = 1.75011
    w/H = 2.0:  D2 = 31.47169   (p=1600 raw 31.53821)   D2/(2 pi w/H) = 2.50444
    w/H = 3.0:  D2 = 52.14107   (p=1600 raw 52.24711)   D2/(2 pi w/H) = 2.76617
  non-linearity of D2 in w (beta_eff variation per eps^2 between w=0.5H and 3H): 1.96991  -> frequency-dependent KMS violation is O(eps^2)
== D. the present universe: exact LambdaCDM worldline correlator vs the first-order thermal form (time domain) ==
  t0 = 0.950985/H0 (a worldline exists only for |Delta| < 2 t0: no full Fourier spectrum at t0 without a window)
  eps_0 = 0.4725, H_eff(t0) = H0 + Hdot/H0 = a''/a' = 0.5275 H0 -> T_eff = 0.08395 H0/2pi (first-order thermal form)
  [PASS]   eps_0 from a(t) directly: got 0.4725 want 0.4725 rel 2.09e-31
  Delta*H0   W_LCDM/W_dS(H_eff) - 1    eps0^2 f2(x_eff) [2nd-order prediction]   actual/2nd-order
    0.25         0.0012174             0.00032262                  3.774
     0.5          0.005563              0.0012792                  4.349
     1.0          0.036387              0.0049379                  7.369
     1.5            0.1916               0.010445                  18.34
  => at eps_0 = 0.47 the deviation from the first-order thermal form is 0.1% (Delta=0.25/H0) to 19% (1.5/H0) and exceeds the
     second-order prediction by factors 4-18: the adiabatic expansion has broken down; the present-day comoving vacuum is
     not thermal at any single temperature and must be treated with the exact kernel (no eps-expansion).
  (formal baseline for the vacuum medium; the physical photon field is dominated by the 2.7 K CMB bath at these scales)

V11 OVERALL: PASS
```

## `v12_two_temperatures.py`

```text
== 12.1 the two temperatures ==
  T_KH  = H*(2 - epsilon)/(4*pi)    T_det = H*(1 - epsilon)/(2*pi)
  T_KH = (H/2pi)(1 - eps/2), T_det = (H/2pi)(1 - eps): PASS
  T_det/T_KH = -epsilon**2/4 - epsilon/2 + 1  -> equal iff eps = 0: PASS
  a''/a' - (H + Hdot/H) = 0 -> PASS
== 12.2 Killing character of the apparent horizon ==
  R_A_dot = eps: the apparent horizon is stationary (Killing) iff eps = 0. For eps != 0 the modular flow of the
  restricted state is not geometric, P1 assigns no temperature to the region inside R_A, and T_KH is a surface-gravity label.
== 12.3 the apparent-horizon first law is an identity of the Einstein equations (flat FRW, perfect fluid) ==
  Cai-Kim:  A(rho+p)H R_A - T0 dS/dt = 0 -> PASS (identity with T0 = 1/2 pi R_A)
  Akbar-Cai: dE/dt - [(kappa/2pi) dS/dt + W dV/dt] = 0 -> PASS (identity with the signed Hayward kappa, |kappa|/2pi = T_KH)
  => both horizon 'first laws' are identities of the Einstein equations GRUT III already has from D3; T_KH (with its sign)
     is the bookkeeping temperature that makes Hayward's unified first law close, not a KMS parameter supplied by P1.
     (The sign/'heat-out' convention is the 'temperature confusion' of Tian-Booth 2014.)
== 12.4 lemma: massless field in power-law FRW == massive field in de Sitter (same nu) ==
  a''/a = (nu**2 - 1/4)/etan**2  vs (nu^2-1/4)/eta^2: PASS
  nu(eps) = (epsilon - 3)/(2*(epsilon - 1))  = 3/2 + eps/(1-eps): PASS
  dS analogue mass: m^2/H^2 = 9/4 - nu^2 = epsilon*(2*epsilon - 3)/(epsilon**2 - 2*epsilon + 1)  (< 0: the dS analogue is tachyonic; the IR needs care, the velocity correlator does not)

CONCLUSION (GRUT III's answer to the two-temperature question):
  (b') They differ because they are different objects. P1's bath temperature is the KMS parameter of the restricted
  state where the modular flow is geometric; for the comoving observer in the conformal vacuum that is T_det = (a''/a')/2pi
  at first order (derived, v11). The apparent horizon is not a Killing horizon for eps != 0, so P1 assigns it no
  temperature; T_KH is Hayward's dynamical surface gravity, used by Cai-Kim as the bookkeeping that turns the horizon
  first law into an identity of the Einstein equations (12.3). No heat-flow law is implied; the O(eps) mismatch equals
  the apparent horizon's expansion rate R_A_dot = eps. Testable statement: a comoving detector coupled to a conformal field
  in the conformal vacuum thermalises (to first order) at (a''/a')/2pi, not at T_KH -- conventional QFT, checkable by anyone.

V12 OVERALL: PASS
```
