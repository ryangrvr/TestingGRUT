"""Fixed ALP cogenesis EFT and the homogeneous pendulum invariants."""
from dataclasses import dataclass, replace
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import ellipe, ellipeinc, ellipk, ellipkinc

HERE = Path(__file__).resolve().parent
Q = 4/np.pi
MASS_EV = 110e-6
T_EW = 130.
Y_B = 8.7e-11
DARK_PER_ENTROPY_EV = .44
C_B = 12/79
A_S = 2.1e-9
G_STAR = 61.75
REDUCED_PLANCK_GEV = 2.435e18
ALPHA = 1/137.035999


@dataclass(frozen=True)
class Config:
    h: float = .01
    n: int = 48
    length: float = 32.
    dt: float = .05
    seed: int = 314
    seed_scale: float = 1.
    a_start: float = .6
    a_end: float = 4.


BASE = Config()
CASES = dict(
    base=BASE,
    spacing=replace(BASE, n=64),
    volume=replace(BASE, n=72, length=48.),
    time=replace(BASE, dt=.025),
    seed=replace(BASE, seed=2718),
    lower_seed=replace(BASE, seed_scale=.1),
    faster=replace(BASE, h=.02),
    slower=replace(BASE, h=.005),
)


def rotation_energy(a):
    """Dimensionless rho/(2 f^2 m^2) at fixed rotational action Q."""
    if not 0 < a < 1:
        raise ValueError("The rotational invariant is used below the reference separatrix.")
    return brentq(lambda r: np.sqrt(r)*ellipe(1/r)-a**-3, 1., 1+a**-6,
                  xtol=3e-14)


def rotation_speed(a, phase=0.):
    return np.sqrt(4*rotation_energy(a)-2*(1-np.cos(phase)))


def libration_action(a, r):
    if not 0 <= r < 1:
        raise ValueError("Libration action requires energy below the barrier.")
    return 8*a**3/np.pi*(ellipe(r)-(1-r)*ellipk(r))


def floquet(momentum, r, rotating=True):
    """Exact static n=1 Lamé growth exponent in units of the ALP mass."""
    p, r = np.broadcast_arrays(momentum, r)
    z = r-p*p if rotating else 1-p*p/r
    keep = (z > 0) & (z < 1)
    out = np.zeros(p.shape)
    parameter = 1/r[keep] if rotating else r[keep]
    amplitude = np.arcsin(np.sqrt(z[keep]))
    zeta = ellipeinc(amplitude, parameter) - (
        ellipe(parameter)/ellipk(parameter))*ellipkinc(amplitude, parameter)
    out[keep] = np.maximum(0., np.sqrt(r[keep])*zeta if rotating else zeta)
    return out


def matching(conversion):
    y_theta = DARK_PER_ENTROPY_EV/(MASS_EV*conversion)
    fa = T_EW*np.sqrt(C_B*y_theta/Y_B)
    return dict(conversion=conversion, Y_theta=y_theta, f_a_gev=fa,
                g_agamma_per_gev=ALPHA/(2*np.pi*fa))


def physical_scales():
    m = MASS_EV*1e-9
    entropy = Q*C_B*T_EW**2*m/Y_B
    temperature = (45*entropy/(2*np.pi**2*G_STAR))**(1/3)
    hubble = np.pi*np.sqrt(G_STAR/90)*temperature**2/REDUCED_PLANCK_GEV
    s_ew = (2*np.pi**2/45)*106.75*T_EW**3
    a_ew = (entropy/s_ew)**(1/3)
    h_ew = np.pi*np.sqrt(106.75/90)*T_EW**2/REDUCED_PLANCK_GEV
    return dict(T_star_gev=temperature, s_star_gev3=entropy, H_star_gev=hubble,
                H_star_over_m=hubble/m, m_over_H_star=m/hubble, a_EW=a_ew,
                aH_EW_over_m=a_ew*h_ew/m, theta_dot_EW_gev=Y_B*s_ew/(C_B*T_EW**2))


def window(k):
    out = np.zeros_like(k)
    nonzero = k > 0
    out[nonzero] = np.exp(-(k[nonzero]/4)**4-(.30/k[nonzero])**4)
    return out


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)
