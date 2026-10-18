"""Frozen exact arithmetic controls for the actual-source analytic proof.

No trajectories, spectral quadrature, parameter searches or floating-point
decisions. Floats in JSON are display only. All comparisons use Fraction.
"""
from fractions import Fraction as Q
from math import isqrt
import json
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

checks = []


def check(name, condition):
    checks.append({"name": name, "pass": bool(condition)})


def record(x):
    x = Q(x)
    return {"numerator": str(x.numerator), "denominator": str(x.denominator),
            "display": float(x)}


class Dual:
    def __init__(self, v, d=0):
        self.v, self.d = Q(v), Q(d)

    def __add__(self, y):
        y = y if isinstance(y, Dual) else Dual(y)
        return Dual(self.v+y.v, self.d+y.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, y):
        return self+-y if isinstance(y, Dual) else self+(-Q(y))

    def __rsub__(self, y):
        return -self+y

    def __mul__(self, y):
        y = y if isinstance(y, Dual) else Dual(y)
        return Dual(self.v*y.v, self.d*y.v+self.v*y.d)

    __rmul__ = __mul__

    def __truediv__(self, y):
        y = y if isinstance(y, Dual) else Dual(y)
        return Dual(self.v/y.v, (self.d*y.v-self.v*y.d)/y.v**2)

    def __rtruediv__(self, y):
        return Dual(y)/self

    def __pow__(self, n):
        if n == 0:
            return Dual(1)
        if n < 0:
            return Dual(1)/(self**(-n))
        return Dual(self.v**n, n*self.v**(n-1)*self.d)


def shell_functions(t):
    s0 = 9*t**2+t**3*(7/(1-t)+2*t/(1-t)**2)
    s1 = 18*t**2+t**3/2*(3*t*(1+t)/(1-t)**3
                         +23*t/(1-t)**2+44/(1-t))
    s11 = 36*t**2+t**3*(t*(1+4*t+t**2)/(1-t)**4
                         +12*t*(1+t)/(1-t)**3
                         +48*t/(1-t)**2+64/(1-t))
    return s0, s1, s11


def sqrt_up(x, digits=12):
    """Exact outward rational square-root enclosure, no float root."""
    x = Q(x)
    d = 10**digits
    n = isqrt(x.numerator*d*d//x.denominator)
    if n*n*x.denominator < x.numerator*d*d:
        n += 1
    y = Q(n, d)
    check("outward_sqrt_"+str(len(checks)),
          y*y >= x and (n == 0 or (y-Q(1, d))**2 < x))
    return y


jp, jm, a0 = Q(1, 64), Q(1, 96), Q(5, 4)
s, p = jp*jp+jm*jm, jp*jm
J, A = 2*p/a0, s/a0
vstar = (jp+jm)**2
u = 4*vstar/a0**2
amax = a0+5*J
check("nominal_J", J == Q(1, 3840))
check("closed_series_u", u == Q(1, 576))
check("positive_dispersion_correction", (2*vstar)**2/(a0**3*(1-u))/J == Q(1, 276))
check("refined_closed_arm_rest", 2*vstar**2/(a0**2*2*(1-u))/J == Q(5, 4416))
check("dispersion_improves_J_over_250", max(Q(1, 276), Q(5, 4416)) < Q(1, 250))

eps = (1-u)**(-2)-1
l = Q(501, 500)
C = l*l
cmin = 1/(1+eps)
j0 = Q(960, 961)*J
check("Wiener_log_converges", eps < 1)
check("one_sided_Wiener_factor", l*l*(1-eps) > 1)
check("J_E_lower", a0/(a0+5*J) == Q(960, 961))

S0, S1, S11 = shell_functions(u)
_, S1dual, _ = shell_functions(Dual(u, 1))
dS1 = S1dual.d
d2 = a0*a0*2
W0num = 2*s*s+4*p*p+4*p*s
W1num = 3*s*s+2*p*p+6*p*s
W11num = 5*s*s-2*p*p+8*p*s
Wrnum = s*s-2*p*p+2*p*s
Cb = Q(1, 460)
M0 = 4*A*J+W0num/d2+Cb*S0
M1 = 4*A*J+W1num/d2+Cb*S1
M11 = 4*A*J+W11num/d2+Cb*S11
Mr = Wrnum/d2+Cb*(S1-S0)
check("direct_row_moment_majorant", W1num-W0num == Wrnum)
for r in range(9):
    pairs = [(n, m) for n in range(r+1) for m in range(r+1) if max(n, m) == r]
    check("shell_count_"+str(r), len(pairs) == 2*r+1)
    check("shell_first_moment_"+str(r),
          sum(m+1 for n, m in pairs) == Q((r+1)*(3*r+2), 2))
    check("shell_product_moment_"+str(r),
          sum((m+1)*(n+1) for n, m in pairs) == (r+1)**3)
check("first_nine_remainder_entries", sum(m+1 for n in range(3) for m in range(3)) == 18)
check("product_first_nine_entries", sum((n+1)*(m+1) for n in range(3) for m in range(3)) == 36)
check("positive_Mr", Mr > 0)
check("M0_small", M0 < J/400)
check("M1_small", M1 < J/390)
check("M11_small", M11 < J/300)
eta = C*M1/j0
check("localized_boundary_Neumann", eta < Q(1, 300))

# Positive-series majorant derivatives F(E,v)=2v/a/(1-4v/a^2).
F10 = 2*vstar/a0**2*(1+u)/(1-u)**2
F01 = 2/a0/(1-u)**2
F20 = 4*vstar/a0**3*(1+3*u)/(1-u)**3
F11 = 2/a0**2*(1+3*u)/(1-u)**3
F02 = 16/a0**3/(1-u)**3
F30 = 12*vstar/a0**4*(1+6*u+u*u)/(1-u)**4
F21 = 4/a0**3*(1+8*u+3*u*u)/(1-u)**4
F12 = 48/a0**4*(1+u)/(1-u)**4
F03 = 192/a0**5/(1-u)**4
P, bmin = 1+F10, 2/amax
RE = (2/a0)*(F20/bmin+P*F11/bmin**2)
Rv = (2/a0)*(F11/bmin+P*F02/bmin**2)
REE = (2/a0)*(F30/bmin+(2*F20*F11+P*F21)/bmin**2
                       +2*P*F11**2/bmin**3)
REv = (2/a0)*(F21/bmin+(F20*F02+F11*F11+P*F12)/bmin**2
                       +2*P*F11*F02/bmin**3)
Rvv = (2/a0)*(F12/bmin+(2*F11*F02+P*F03)/bmin**2
                       +2*P*F02**2/bmin**3)
Rmin = 1/(1+eps)
Rmax = amax/a0*(1+F10)
V = Q(201, 100)*J
check("band_velocity_upper", 2/Rmin < Q(201, 100))
check("band_velocity_lower", 2/Rmax > Q(199, 100))
R1 = RE*V+Rv*2*p
E2 = 2*J*(1/Rmin+R1/Rmin**2)
R2 = REE*V**2+2*REv*V*(2*p)+Rvv*(2*p)**2+RE*E2+Rv*(2*p)
check("R_first_derivative", R1 < Q(1, 100))
check("R_second_derivative", R2 < Q(1, 100))
check("velocity_increases_left", Q(1, 2)/Rmax-R1/Rmin**2 > 0)
check("velocity_decreases_right", Q(1, 2)/Rmax-R1/Rmin**2 > 0)
check("one_velocity_maximum_middle", Q(4, 5)/Rmax
      -2*R1/Rmin**2-R2/Rmin**2-2*R1**2/Rmin**3 > 0)
Iconst = (2*Rmax+5*J*RE+4*p*Rv)/Q(31, 10)
check("actual_dispersion_return_integral", Iconst < 1)

Mdot1 = V*(s/a0**2+W1num/d2*(2/a0+Q(1, 2))
           +Cb/a0*((1+2*u/(1-u))*S1+2*u*dS1))
check("boundary_parameter_derivative", Mdot1 < J/2000)
AE = 2/a0*((1+u)/(1-u)**3-1)
Av = u/(vstar*(1-u)**3)
Am = 2*p*Av
Adot = AE*V+Av*2*p
Ldot = l*Adot/(2*(1-eps))
Lm = l*Am/(2*(1-eps))
cdot = Adot/(2*(1-eps))
jrel = V/a0
Cder = 2*l*Ldot+C*jrel
B = (Cder*M1+C*M11+C*Mdot1)/j0
s0 = Q(13, 20)
tdot = Ldot+l*(1/s0+1/s0**2)
x = l/(1-eta)
xdot = (tdot+B*x)/(1-eta)
z0 = eta*x/cmin
z1 = (B*x+eta*xdot)/cmin+z0*cdot/cmin
check("boundary_jet_finite", B < Q(1, 200))
check("relative_source_jump", z0 < Q(1, 300))
check("relative_source_jump_derivative_central", z1 < Q(1, 50))

R0 = C*C*M0/(j0*j0*(1-eta))
alphaerr = max(eps/(1+eps), 4*J)
derr = Q(1, 276)+alphaerr+J*R0
check("fuel_denominator_nonzero_whole_band", derr < Q(1, 100))
fmax = 1/(1-derr)
Rdiff_dot = (B*C/(1-eta)**2+eta*(C+Cder)/(1-eta))/j0
alphadot = (1+4*J)*(jrel+Adot/(1-eps))
banddot = 2*max(1/Rmin-1, 1-1/Rmax)
dder = banddot+alphaerr+alphadot+J*Rdiff_dot
fprime = (dder+derr)/(1-derr)**2
Tprime = Adot+Am+(1+eps)*(F20*V+F11*(2*p))
Zmax = (1+eps)*fmax*l*(1+z0)
Zprime = (Tprime*fmax*l*(1+z0)
          +(1+eps)*fprime*l*(1+z0)
          +(1+eps)*fmax*(Ldot+Lm+cdot*l)*(1+z0)
          +(1+eps)*fmax*l*z1)
rmin = (1-z0)**2/((1+F10)*(1+eps)*(1+derr)**2)
rmax = (1+eps)*fmax*fmax*(1+z0)**2
rprime = (Tprime*fmax*fmax*(1+z0)**2
          +(1+eps)*2*cdot*fmax*fmax*(1+z0)**2
          +(1+eps)*2*fmax*fprime*(1+z0)**2
          +(1+eps)*fmax*fmax*2*(1+z0)*z1)
check("actual_density_lower", rmin > Q(24, 25))
check("actual_density_upper", rmax < Q(26, 25))
check("actual_density_parameter_derivative", rprime < Q(1, 8))
check("return_for_full_retention_window", Q(2, 25)+Q(1, 2) == Q(29, 50) < Q(2, 3))

# Exact central interval and mass coverage; alternating Taylor bounds.
kstart = Q(3, 4)
cos_upper = 1-kstart*kstart/2+kstart**4/24
sin_arg = Q(3, 2)
sin_lower = sin_arg-sin_arg**3/6+sin_arg**5/120-sin_arg**7/5040
pi_lo, pi_hi = Q(31, 10), Q(22, 7)
check("central_interval_contains_predeclared_k_segment", 2*cos_upper+Q(1, 276) < Q(3, 2))
free_central_mass_lower = 1-(Q(3, 2)-sin_lower)/pi_lo
check("central_free_mass", free_central_mass_lower > Q(5, 6))
check("actual_central_mass", Q(24, 25)*free_central_mass_lower > Q(4, 5))
central_cos = Q(94, 125)
check("central_sine_bound", 1-central_cos**2 > s0*s0)
check("central_velocity_1_3", Q(199, 100)**2*(1-central_cos**2) > Q(13, 10)**2)

Mb = Zmax/pi_lo
TVb = (Zmax*(2-2*s0)+2*Zprime)/pi_lo
phase_constant = 4*Mb+TVb
check("phase_sensitive_amplitude_variation", phase_constant < Q(8, 5))
Mw = 2*rmax/pi_lo
TVw = 2/pi_lo*(2*rmax+rprime*pi_hi)
fuel_central_constant = (4*Mw+TVw)/Q(13, 10)
check("central_fuel_return", fuel_central_constant < 10)
Rloc = Lm*(1+z0)/j0+C/(j0*s0)*C*Mr/(j0*(1-eta))
local_packet = 3*J/pi_lo*fmax*Rloc
check("global_l2_local_packet_remainder", local_packet < Q(1, 200))

# Check the exact free remainder identity and derivative phase moment.
for n in range(4):
    for m in range(4):
        actual = {n+m+1-2*j: Q(-1) for j in range(min(n, m)+1)}
        outgoing = {n+m+1-2*j: Q(-1) for j in range(m+1)}
        if m > n:
            for j in range(m-n):
                exponent = m-n-1-2*j
                outgoing[exponent] = outgoing.get(exponent, Q(0))+1
        outgoing = {e: c for e, c in outgoing.items() if c}
        check("exact_free_global_remainder_"+str(n)+"_"+str(m), actual == outgoing)
        check("free_derivative_product_moment_"+str(n)+"_"+str(m),
              sum(n+m+1-2*j for j in range(min(n, m)+1)) == (n+1)*(m+1))

theta = Q(252)
out_tail_square = 1/(Q(3, 10)*theta-2)+1/(Q(6, 25)*theta+1)
in_tail_square = 1/(Q(13, 10)*theta+1)
tail_upper = Q(8, 5)*(sqrt_up(out_tail_square)+sqrt_up(in_tail_square))+Q(1, 200)
central_loss_upper = tail_upper**2+(10/theta)**2+Q(1, 362)
check("continuous_scaling_window_capture", central_loss_upper < Q(3, 20))
check("coherent_projection_first_root", Q(4, 5)-Q(3, 20) > Q(4, 5)**2)
check("coherent_projection_second_root", (1-Q(4, 5))*Q(3, 20)/Q(4, 5) < Q(1, 5)**2)
source_A_lower = (Q(4, 5)-Q(1, 5))**2
check("actual_original_source_capture", source_A_lower == Q(9, 25) > Q(17, 50))

sigma = Q(6, 27613)
omega = sigma*sigma
delta = Q(1, 1000)
check("one_source_to_joint_amplitude", source_A_lower-2*sigma-delta > Q(717, 2000))
check("retention_original_source_account", 1-omega-(Q(29, 50)+omega)**2-Q(1, 1000)**2 > Q(33, 50))
joint_lower = Q(717, 2000)**2
retention_lower = Q(659, 1000)**2
check("joint_probability_margin", joint_lower == Q(514089, 4000000) > Q(1, 10))
check("retention_probability_margin", retention_lower == Q(434281, 1000000) > Q(1, 4))

bounds = {name: record(value) for name, value in {
    "epsilon_Wiener": eps, "M0_over_J": M0/J, "M1_over_J": M1/J,
    "M11_over_J": M11/J, "Mr_over_J": Mr/J, "Mdot1_over_J": Mdot1/J,
    "eta": eta, "Rmin": Rmin, "Rmax": Rmax, "R1": R1, "R2": R2,
    "return_integral_constant": Iconst, "L_dot": Ldot, "L_first_moment": Lm,
    "boundary_jet_B": B, "zeta_max": z0, "zeta_prime": z1,
    "fuel_denominator_error": derr, "fuel_denominator_error_derivative": dder,
    "Z_max": Zmax, "Z_prime": Zprime, "density_ratio_lower": rmin,
    "density_ratio_upper": rmax, "density_ratio_prime": rprime,
    "central_mass_lower": Q(24, 25)*free_central_mass_lower,
    "phase_constant": phase_constant, "fuel_central_constant": fuel_central_constant,
    "global_local_packet_remainder": local_packet, "central_loss_upper": central_loss_upper,
    "one_source_A_lower": source_A_lower, "joint_probability_lower": joint_lower,
    "retention_probability_lower": retention_lower,
}.items()}
ok = all(c["pass"] for c in checks)
print(json.dumps({"status": "all_frozen_controls_pass" if ok else "frozen_control_failure",
                  "arithmetic": "exact Fractions; rational outward square roots; display floats noncertifying",
                  "role": "analytic constant and identity controls, no physical trajectory",
                  "checks": checks, "bounds": bounds,
                  "operating_gate": "author_derived_pass_pending_independent_proof_review" if ok else "unresolved"},
                 indent=2, sort_keys=True))
raise SystemExit(0 if ok else 1)
