"""Scale and stress estimates for the frozen production histories."""
import json

import numpy as np
from scipy.integrate import simpson

from transfer import (A_S, DM_PER_ENTROPY, GEV_TO_INVERSE_MPC, G_S0, HISTORIES,
                      IR_RADIATION, K_PIVOT, M, MPL, RESULTS, T_CMB, T_STAR,
                      entropy, write)


def estimates(history, transfer, fit):
    name = history.name
    q_end = fit["a_end"]*fit["H_end_gev"]*GEV_TO_INVERSE_MPC/history.k_star
    h_reheat = history.h_reheat or fit["H_end_gev"]/M
    q_reheat = (history.x_reheat or h_reheat**(-.5))*h_reheat
    integral = transfer["fine"]["integral"][-1]
    r_flat, r = transfer["flat_r"], fit["tensor_ratio"]
    entropy_readout_temperature = history.reheat_temperature if h_reheat < 1 else T_STAR
    x = history.x_reheat if h_reheat < 1 else 1.
    plateau_h = MPL*np.sqrt(fit["V0_over_Mpl4"]/3)
    with np.load(RESULTS / (name+"_transfer.npz")) as packet:
        q, weight = packet["q"], packet["weight"]
    exponent, coefficient = history.uv_power
    boosted = simpson(q*weight, x=np.log(q))
    if exponent == 1:
        boosted_tail = coefficient*np.log(q_end/q[-1])
    else:
        beta = 1-exponent
        boosted_tail = coefficient*(q_end**beta-q[-1]**beta)/beta
    # Replace the true red spectrum with its pivot value for the boost estimate.
    # Double the transfer-tail coefficient to allow a margin near matching.
    stress_ratio = ((4*DM_PER_ENTROPY/(3*entropy_readout_temperature))
                    * (1+(r/r_flat)*(boosted+2*boosted_tail)/(integral*x)))
    tensor_stress = None
    if history.w > 1/3:
        beta = (3*history.w-1)/(history.p-1)
        tensor_stress = r*A_S/12*((q_end/q_reheat)**beta-1)/beta
    return dict(average_tensor_power_over_pivot=r_flat*fit["density_ratio"]/r,
                tensor_weight_above_reheat_RD_estimate=float(1/(integral*q_reheat))
                if h_reheat > 1 else None,
                k_star_over_pivot=history.k_star/K_PIVOT, q_end=float(q_end),
                q_reheat=float(q_reheat), H_reheat_over_m=h_reheat,
                min_omega_over_H_at_reheat=(1/h_reheat if h_reheat < 1 else None),
                N_from_end_to_reheat=(np.log(fit["H_end_gev"]/(M*h_reheat))/history.p),
                plateau_H_gev=float(plateau_h), plateau_H4_over_total=(
                    plateau_h**2/(3*MPL**2)),
                vector_stress_ratio_estimate=float(stress_ratio),
                stress_temperature_gev=entropy_readout_temperature,
                stiff_tensor_stress_ratio_estimate=tensor_stress)


def main():
    transfer = json.loads((RESULTS / "transfer.json").read_text())["histories"]
    inflation = json.loads((RESULTS / "inflation.json").read_text())
    late_extension = json.loads((RESULTS / "late_stiff_extension.json").read_text())
    fits = {name: row["fine"] for name, row in inflation["fitted_histories"].items()
            if row["fine"]["root_bracketed"]}
    if late_extension["passed"]:
        fits["late_stiff"] = late_extension["result"]["fine"]
    records = {name: estimates(HISTORIES[name], transfer[name], fit)
               for name, fit in fits.items()}
    # Unit conversion from Omega_c h^2=0.120 to energy per current entropy.
    hubble_100 = 100/3.0856775814913673e19*6.582119569e-25  # km/s/Mpc to GeV
    critical_density = 3*MPL*MPL*hubble_100*hubble_100
    present_entropy = entropy(T_CMB, G_S0)
    prompt_integral = transfer["prompt"]["fine"]["integral"][-1]
    earlier_fraction = (IR_RADIATION/2*(K_PIVOT/HISTORIES["prompt"].k_star)**2
                        / prompt_integral)
    result = dict(histories=records, H100_gev=hubble_100,
                  rho_critical_over_h2_gev4=critical_density,
                  present_entropy_gev3=present_entropy,
                  energy_per_entropy_at_Omega_c_h2_0120_eV=.12*critical_density
                  /present_entropy*1e9,
                  pre_pivot_RD_weight_fraction=float(earlier_fraction),
                  pre_pivot_absolute_r_correction_if_Pt_lt_1=float(earlier_fraction/A_S))
    write(RESULTS / "diagnostics.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
