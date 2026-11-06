"""Retrospective, leading-order ALP cogenesis screen, with no detector fit."""
import json
from math import isclose, pi, sqrt
from pathlib import Path

MASS_UEV = 110.
ELECTROWEAK_TEMPERATURE = 130.
BARYON_YIELD = 8.7e-11
DARK_MASS_PER_ENTROPY_EV = .44
TRANSFER = 2.
BARYON_COEFFICIENT = 12/79
PHOTON_COEFFICIENT = 1.
INVERSE_ALPHA = 137.035999
LOCAL_DARK_DENSITY = .45
REFERENCE_COUPLING_LIMIT = 4e-12
RELATIVISTIC_DEGREES = 106.75


def point(cosmic_fraction, local_fraction):
    mass_ev = MASS_UEV*1e-6
    charge_yield = cosmic_fraction*DARK_MASS_PER_ENTROPY_EV/(TRANSFER*mass_ev)
    f_a = ELECTROWEAK_TEMPERATURE*sqrt(BARYON_COEFFICIENT*charge_yield/BARYON_YIELD)
    coupling = PHOTON_COEFFICIENT/(2*pi*INVERSE_ALPHA*f_a)
    entropy = (2*pi*pi/45)*RELATIVISTIC_DEGREES*ELECTROWEAK_TEMPERATURE**3
    radiation = (pi*pi/30)*RELATIVISTIC_DEGREES*ELECTROWEAK_TEMPERATURE**4
    speed = charge_yield*entropy/f_a**2
    return dict(cosmic_fraction=cosmic_fraction, local_fraction=local_fraction,
                charge_yield=charge_yield, decay_constant_gev=f_a,
                photon_coupling_per_gev=coupling,
                rotation_rate_at_ew_gev=speed, rotation_rate_over_mass=speed/(mass_ev*1e-9),
                rotation_energy_fraction_at_ew=.5*f_a*f_a*speed*speed/radiation,
                halo_power_factor=coupling**2*local_fraction*LOCAL_DARK_DENSITY/(mass_ev*1e-9),
                power_over_quoted_sensitivity=local_fraction*(coupling/REFERENCE_COUPLING_LIMIT)**2)


def main():
    sample = point(1., 1.)
    fractions = [point(f, f) for f in (1., .3, .1)]
    assert all(isclose(p["halo_power_factor"], sample["halo_power_factor"], rel_tol=1e-14)
               for p in fractions)
    coupling = sample["photon_coupling_per_gev"]
    result = dict(mass_uev=MASS_UEV, frequency_ghz=MASS_UEV*1e-6*241798.9242,
                  c_W=-1., c_Y=0., c_B=BARYON_COEFFICIENT, c_gamma=-PHOTON_COEFFICIENT,
                  T_EW_gev=ELECTROWEAK_TEMPERATURE, baryon_yield=BARYON_YIELD,
                  dark_mass_per_entropy_ev=DARK_MASS_PER_ENTROPY_EV, transfer=TRANSFER,
                  g_limit_screen_per_gev=REFERENCE_COUPLING_LIMIT,
                  rho_dm_local_gev_per_cm3=LOCAL_DARK_DENSITY,
                  all_dark_matter=sample, comoving_fraction_controls=fractions,
                  transfer_at_quoted_sensitivity=TRANSFER*(REFERENCE_COUPLING_LIMIT/coupling)**2)
    path = Path(__file__).resolve().parent / "results" / "cogenesis_screen.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
