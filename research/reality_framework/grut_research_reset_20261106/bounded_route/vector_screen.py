"""Radiation transfer for the frozen inflation-produced Proca comparator."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.special import gamma

DIRECTORY = Path(__file__).resolve().parent / "results"
MASS_EV = .001
MASS_GEV = MASS_EV*1e-9
MPL = 2.435e18
G = 106.75
SCALAR_POWER = 2.1e-9
CDM_ENERGY_PER_ENTROPY = .44e-9
PIVOT_MPC = .05
GEV_PER_MPC_INVERSE = 1/1.563738305e38
T_CMB_GEV = 2.34865e-13
G_ENTROPY_TODAY = 3.91


def radiation_transfer(q, start, stops, rtol):
    """u=A_L/A_L,in, in x=a/a_* with H(a_*)=m."""
    k2 = q*q

    def rhs(a, flat):
        value, speed = flat.reshape(2, -1)
        total = a*a+k2
        return np.r_[speed, -2*k2/(a*total)*speed-total*value]

    initial = np.r_[1-k2*start**2/6, -k2*start/3]
    sol = solve_ivp(rhs, (start, max(stops)), initial, t_eval=stops, method="DOP853",
                    rtol=rtol, atol=rtol/50)
    if not sol.success:
        raise RuntimeError(sol.message)
    action = []
    for a, y in zip(stops, sol.y.T):
        value, speed = y.reshape(2, -1)
        total = a*a+k2
        canonical = a/np.sqrt(total)*value
        canonical_speed = a/np.sqrt(total)*speed+k2/total**1.5*value
        frequency = np.sqrt(total+3*k2/total**2)
        action.append((canonical_speed**2+frequency**2*canonical**2)/frequency)
    return np.array(action), sol.nfev


def screen(q_min, q_max, points, start, stops, rtol):
    q = np.geomspace(q_min, q_max, points)
    twice_action, evaluations = radiation_transfer(q, start, stops, rtol)
    ir_limit = 2*gamma(.75)**2/np.pi
    ir_tail = .5*ir_limit*q_min**2
    uv_tail = 1/q_max
    densities = q*q*twice_action
    integral = simpson(densities, x=np.log(q), axis=1)+ir_tail+uv_tail
    final = densities[-1]
    peak = int(np.argmax(final))
    row = dict(q_min=q_min, q_max=q_max, points=points, starting_a=start,
               readout_a=stops, evaluations=evaluations, transfer_integrals=integral.tolist(),
               ir_coefficient=float(twice_action[-1, 0]), ir_reference=float(ir_limit),
               uv_coefficient=float(q_max**3*twice_action[-1, -1]),
               ultraviolet_tail=uv_tail, peak_q=float(q[peak]))
    return row, q, final


def main():
    coarse, _, _ = screen(.001, 64., 321, 1e-5, (12., 18.), 4e-9)
    fine, q, spectrum = screen(.0005, 128., 641, 5e-6, (16., 24.), 3e-10)
    transfer = fine["transfer_integrals"][-1]
    error = abs(coarse["transfer_integrals"][-1]/transfer-1)
    late_drift = abs(fine["transfer_integrals"][0]/transfer-1)
    ir_error = abs(fine["ir_coefficient"]/fine["ir_reference"]-1)
    uv_error = abs(fine["uv_coefficient"]-1)
    if max(error, late_drift) >= .005 or max(ir_error, uv_error) >= .015:
        raise RuntimeError("The declared transfer comparison missed.")
    temperature = np.sqrt(MASS_GEV*MPL/(np.pi*np.sqrt(G/90)))
    entropy = 2*np.pi**2*G*temperature**3/45
    h_req = np.sqrt(8*np.pi**2*CDM_ENERGY_PER_ENTROPY*entropy/(MASS_GEV**2*transfer))
    r_req = 2*h_req*h_req/(np.pi**2*MPL*MPL*SCALAR_POWER)
    a_today_norm = T_CMB_GEV/temperature*(G_ENTROPY_TODAY/G)**(1/3)
    k_star = a_today_norm*MASS_GEV/GEV_PER_MPC_INVERSE
    h_paper = 1e14*(6e-6/MASS_EV)**.25
    result = dict(mass_eV=MASS_EV, assumed_scalar_power=SCALAR_POWER,
                  coarse=coarse, fine=fine, grid_shift=error, late_shift=late_drift,
                  infrared_error=ir_error, ultraviolet_error=uv_error,
                  production_temperature_gev=temperature, k_star_per_mpc=k_star,
                  log_kstar_over_pivot=float(np.log(k_star/PIVOT_MPC)),
                  required_flat_H_gev=h_req, tensor_ratio_for_flat_H=r_req,
                  published_fit_H_gev=h_paper, independent_H_ratio=h_req/h_paper)
    np.savez_compressed(DIRECTORY / "vector_spectrum.npz", q=q, spectrum=spectrum)
    (DIRECTORY / "vector_transfer.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
