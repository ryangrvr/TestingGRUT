"""Fixed instant-preheating analytic screen, in reduced Planck units."""
from math import exp, expm1, log, pi, sqrt
from pathlib import Path
import json

from scipy.integrate import quad
from scipy.optimize import brentq

M_PLANCK_GEV = 2.435e18
H_CROSS = 1.44e-6
H_PIVOT = 5.48e-6
TILT = 0.968
G = 1e-5
WIDTH_AT_GM = 1e-10
G_THERMAL = 106.75
PHI_SPEED_PER_EFOLD = sqrt(6.)
PRODUCTION_LIMIT = .1

RHO_CROSS = 3*H_CROSS**2
K_STAR = sqrt(G*PHI_SPEED_PER_EFOLD*H_CROSS)
CHI_NUMBER = K_STAR**3/(8*pi**3)
YUKAWA = sqrt(8*pi*WIDTH_AT_GM/G)


def mass(efolds):
    return G*PHI_SPEED_PER_EFOLD*efolds


def optical_depth(efolds):
    exponent = 3*efolds
    return (WIDTH_AT_GM*PHI_SPEED_PER_EFOLD/(9*H_CROSS)
            * (exponent*exp(exponent)-expm1(exponent)))


def optical_depth_slope(efolds):
    return WIDTH_AT_GM*PHI_SPEED_PER_EFOLD/H_CROSS*efolds*exp(3*efolds)


def efolds_at_depth(target):
    return brentq(lambda efolds: optical_depth(efolds)-target, 0., 10., xtol=1e-13)


def radiation_comoving(depth):
    """Deposit the rest energy at decay and redshift daughters as radiation."""
    def integrand(tau):
        efolds = efolds_at_depth(tau)
        return mass(efolds)*exp(efolds-tau)
    return CHI_NUMBER*quad(integrand, 0., depth, epsabs=1e-14, epsrel=1e-10)[0]


def equality(comoving_radiation):
    efolds = .5*log(RHO_CROSS/comoving_radiation)
    radiation_density = comoving_radiation*exp(-4*efolds)
    formal_temperature = M_PLANCK_GEV*(30*radiation_density/(pi*pi*G_THERMAL))**.25
    return efolds, formal_temperature


def loop_to_kinetic(efolds):
    """Real-scalar MS-bar term at the declared fixed mu=g Mbar."""
    mass_chi = mass(efolds)
    vacuum_term = mass_chi**4/(64*pi*pi)*(2*log(mass_chi/G)-1.5)
    return vacuum_term/(RHO_CROSS*exp(-6*efolds))


def main():
    max_depth = 60.
    comoving = radiation_comoving(max_depth)
    end = efolds_at_depth(max_depth)
    direct = CHI_NUMBER*quad(
        lambda n: mass(n)*exp(n-optical_depth(n))*optical_depth_slope(n),
        0., end, epsabs=1e-14, epsrel=1e-10)[0]
    integral_difference = abs(direct/comoving-1)
    assert integral_difference < 1e-8

    eq, temperature = equality(comoving)
    median = efolds_at_depth(log(2))
    nonrel = brentq(lambda n: mass(n)-K_STAR*exp(-n), 0., 3.)
    loop = brentq(lambda n: loop_to_kinetic(n)-1, 1., eq)
    tau_at_loop = optical_depth(loop)
    radiation_at_loop = radiation_comoving(tau_at_loop)*exp(2*loop)/RHO_CROSS
    chi_at_loop = CHI_NUMBER*mass(loop)*exp(3*loop-tau_at_loop)/RHO_CROSS

    # Recover the source's constant-mass sudden-decay approximation separately.
    sudden_decay = log(H_CROSS/WIDTH_AT_GM)/3
    constant_width_comoving = CHI_NUMBER*G*exp(sudden_decay)
    _, source_temperature = equality(constant_width_comoving)

    eps_pivot = 3*(1-TILT)**2/16
    ratio = sqrt(3*eps_pivot/4)
    pivot_field = log(ratio/(1+ratio))/sqrt(2/3)
    chi_count_hubble = 4*pi*CHI_NUMBER/(3*H_CROSS**3)
    count_decays_per_hubble_at_loop = (chi_count_hubble*exp(6*loop)
                                       * -expm1(-tau_at_loop))
    result = dict(
        source=dict(reference="arXiv:2004.11843, sections 3 and 4.3",
                    H_cross_over_Mbar=H_CROSS, H_pivot_over_Mbar=H_PIVOT,
                    tree_pivot_tilt=TILT, coupling_g=G, width_at_m_gM_over_Mbar=WIDTH_AT_GM,
                    yukawa_y=YUKAWA),
        production=dict(k_star_over_H_cross=K_STAR/H_CROSS,
                        crossing_time_times_H=H_CROSS/K_STAR,
                        flat_burst_control_pass=H_CROSS/K_STAR < PRODUCTION_LIMIT,
                        initial_chi_count_per_hubble_sphere=chi_count_hubble,
                        mediator_mass_over_H_at_pivot=G*abs(pivot_field)/H_PIVOT),
        tree_level_formal_history=dict(
            median_decay_efolds=median,
            optical_depth_before_typical_momentum_equals_mass=optical_depth(nonrel),
            typical_momentum_over_mass_at_median=K_STAR*exp(-median)/mass(median),
            comoving_daughter_density_over_initial_kinetic_density=comoving/RHO_CROSS,
            equality_efolds=eq, thermal_equivalent_at_equality_GeV=temperature,
            constant_width_source_approximation_GeV=source_temperature,
            radiation_integral_change_of_variable_error=integral_difference),
        loop_control=dict(
            finite_term="Real scalar MS-bar at mu=g Mbar with no compensating finite tail terms",
            unity_efolds=loop, ratio_on_formal_tree_equality=loop_to_kinetic(eq),
            m_chi_over_H_at_unity=mass(loop)/(H_CROSS*exp(-3*loop)),
            decayed_fraction_at_unity=-expm1(-tau_at_loop),
            decays_in_comoving_volume_of_hubble_sphere_at_unity=count_decays_per_hubble_at_loop,
            tree_radiation_to_kinetic_at_unity=radiation_at_loop,
            tree_mediator_to_kinetic_at_unity=chi_at_loop,
            tree_trajectory_control_pass=loop >= eq),
    )
    path = Path(__file__).with_name("results.json")
    path.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(path.read_text(), end="")


if __name__ == "__main__":
    main()
