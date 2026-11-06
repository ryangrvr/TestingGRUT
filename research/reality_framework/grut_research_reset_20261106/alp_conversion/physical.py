"""Physical hierarchy and linear onset diagnostics before backreaction."""
import json

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import logsumexp

from model import (A_S, C_B, HERE, MASS_EV, floquet, matching, physical_scales,
                   rotation_energy, rotation_speed, window, write_json)


def onset(na, nk):
    scales = physical_scales()
    h = scales["H_star_over_m"]
    a = np.linspace(.25, .76, na)
    k = np.geomspace(.2, 12, nk)
    logk = np.log(k)
    weights = .5*(np.r_[np.diff(logk), 0.]+np.r_[0., np.diff(logk)])
    rho = np.array([rotation_energy(value) for value in a])
    omega = np.sqrt(1+(k[None, :]/a[:, None])**2)
    growth = floquet(k[None, :]/a[:, None], rho[:, None])
    exponent = cumulative_trapezoid(2*a[:, None]*growth/h, x=a, axis=0, initial=0.)
    with np.errstate(divide="ignore"):
        # A phase-mixed covariance gains cosh(2G)-1, with G=int mu dt.
        log_gain = exponent-np.log(2.)+2*np.log(-np.expm1(-exponent))
    initial_v = rotation_speed(a[0])
    log_common = np.log(omega*(a[0]/a[:, None])**3)+log_gain
    seeds = {"classical": A_S*initial_v**2*window(k)/omega[0]}
    for c in (.5, 2., 3.):
        fa = matching(c)["f_a_gev"]
        mass_gev = MASS_EV*1e-9
        seeds[f"vacuum_C{c:g}"] = (mass_gev/fa)**2*k**3/(4*np.pi**2*a[0]**3)
    rows = {}
    for name, seed in seeds.items():
        density_per_log = log_common+np.log(seed)
        fraction = logsumexp(density_per_log+np.log(weights), axis=1)-np.log(2*rho)
        crossing = int(np.flatnonzero(fraction >= np.log(.01))[0])
        peak = int(np.argmax(density_per_log[crossing]))
        rows[name] = dict(a_onset=float(np.interp(np.log(.01), fraction[crossing-1:crossing+1],
                                                 a[crossing-1:crossing+1])),
                          peak_comoving_k_over_m=float(k[peak]),
                          peak_physical_p_over_m=float(k[peak]/a[crossing]),
                          rotation_energy_over_barrier=float(rho[crossing]),
                          log_power_gain=float(exponent[crossing, peak]),
                          peak_k_over_aH_at_EW=float(k[peak]/scales["aH_EW_over_m"]))
    return rows


def check_sweep(k, h):
    def edge(a):
        return rotation_energy(a)-(k/a)**2

    enter = brentq(lambda a: edge(a)-1, .25, .8)
    leave = brentq(edge, enter, .8)
    width = leave-enter
    start, stop = enter-.6*width, leave+.6*width
    tmax = (stop-start)/h
    wi = np.sqrt(k*k+start*start)
    initial = [0., start*rotation_speed(start), wi**-.5, 0., 0., wi**.5]

    def rhs(t, y):
        a = start+h*t
        theta, speed, x1, v1, x2, v2 = y
        frequency2 = k*k+a*a*np.cos(theta)
        return [speed, -2*h/a*speed-a*a*np.sin(theta),
                v1, -frequency2*x1, v2, -frequency2*x2]

    sol = solve_ivp(rhs, (0., tmax), initial, t_eval=[tmax], method="DOP853",
                    rtol=3e-11, atol=1e-13, max_step=.2)
    if not sol.success:
        raise RuntimeError(sol.message)
    _, _, x1, v1, x2, v2 = sol.y[:, -1]
    wf = np.sqrt(k*k+stop*stop)
    direct = np.log(.5*((x1*x1+x2*x2)*wf+(v1*v1+v2*v2)/wf))
    exponent = 2*quad(lambda a: a/h*float(floquet(k/a, rotation_energy(a))),
                      enter, leave, epsabs=1e-8)[0]
    return dict(k_over_m=k, a_enter=enter, a_leave=leave, direct_log_energy_gain=direct,
                adiabatic_log_energy_gain=float(np.logaddexp(exponent, -exponent)-np.log(2.)),
                period_count_approx=float((leave-enter)/h*k/np.pi))


def radiative_size():
    fa = matching(2.)["f_a_gev"]
    alpha_w = .0338
    # Diagnostic only: if the zero fermion couplings were imposed at fa.
    c_left = -(9/2)*(alpha_w/(4*np.pi))**2*np.log(fa/246.)
    delta_cb = 3*(18/79+25/237)*c_left
    m = MASS_EV*1e-9
    return dict(C2_scale_gev=fa, alternate_UV_zero_cq_cl=c_left,
                alternate_UV_relative_cB_shift=delta_cb/C_B,
                electron_derivative_loop_kinematic_ratio=(m/.00051099895)**2,
                W_loop_kinematic_ratio=(m/80.4)**2,
                boundary_condition="c_f=0 at mu=130 GeV; the UV log is a size diagnostic")


def main():
    first, second = onset(4097, 401), onset(8193, 801)
    scales = physical_scales()
    k = second["vacuum_C2"]["peak_comoving_k_over_m"]
    result = dict(scales=scales, coarse=first, fine=second, sweep=check_sweep(k, scales["H_star_over_m"]),
                  radiative_size=radiative_size())
    write_json(HERE / "results" / "physical_linear.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
