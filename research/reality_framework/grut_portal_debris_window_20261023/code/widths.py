#!/usr/bin/env python3
"""Declared SM partial-width models and the action's hidden decays."""
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HBAR, HBARC = 6.582119569e-25, 1.973269804e-16  # GeV s, GeV m
VH, MH = 246., 125.
ME, MMU = .00051099895, .1056583755
ALPHA = 1/137.035999084
G, HD = .1, .01


class Widths:
    def __init__(self):
        d = ROOT/'data'
        self.hips = {(p, m): np.loadtxt(d/f'hipso_{p}_{m}.csv', delimiter=',', skiprows=1)
                     for p in ('pi', 'K') for m in ('DGL', 'BTPZ')}
        self.foresee = {name: np.loadtxt(d/f'foresee_{name}.txt')
                        for name in ('ctau', 'br_pi_pi', 'br_K_K', 'br_e_e', 'br_mu_mu')}

    def lepton(self, mass):
        return sum(ml**2*mass/(8*np.pi*VH**2)*max(1-4*ml**2/mass**2, 0)**1.5 for ml in (ME, MMU))

    def unit(self, mass, method='DGL', band='central', photon_coefficient=10.):
        if not .041 <= mass <= 2:
            raise ValueError('Outside the declared SM-width mass range')
        column = {'central': 1, 'lower': 2, 'upper': 3}[band]
        hadronic = 0.
        for p, threshold in (('pi', .268), ('K', .994)):
            a = self.hips[p, method]
            # Preserve the package's isospin-limit masses; do not smear a
            # nonzero tabulated width below its actual threshold.
            hadronic += np.interp(mass, a[:, 0], a[:, column]) if mass > threshold else 0.
        tab = self.foresee['ctau']
        full_old = HBARC/np.exp(np.interp(np.log(mass), np.log(tab[:, 0]), np.log(tab[:, 1])))
        branch_sum = sum(np.interp(mass, self.foresee[n][:, 0], self.foresee[n][:, 1])
                         for n in ('br_pi_pi', 'br_K_K', 'br_e_e', 'br_mu_mu'))
        # Phenomenological multihadron remainder, not a new dispersive result.
        remainder = full_old*max(1-branch_sum, 0.) if mass > .5 else 0.
        photons = photon_coefficient**2*ALPHA**2*mass**3/(256*np.pi**3*VH**2)
        return self.lepton(mass)+hadronic+remainder+photons

    def ensemble(self, mass):
        vals = [self.unit(mass, m, b, 0 if b == 'lower' else 10)
                for m in ('DGL', 'BTPZ') for b in ('lower', 'central', 'upper')]
        return min(vals), self.unit(mass), max(vals)


def hidden(v, lam):
    ms = np.sqrt(2*lam)*v
    mx, md = G*v, np.sqrt(G**2+HD**2)*v
    # s X X*: distinct final particles; s D D: identical real daughters.
    # Evaluate thresholds from dimensionless couplings so roundoff in v does
    # not spuriously open a decay exactly on its boundary.
    gx = G**4*v*v/(4*np.pi*ms)*np.sqrt(max(1-2*G**2/lam, 0.))
    gd = HD**4*v*v/(8*np.pi*ms)*np.sqrt(max(1-2*(G**2+HD**2)/lam, 0.))
    return float(gx), float(gd)


def visible_coefficient(v, lam, unit_width):
    ms2 = 2*lam*v*v
    return (v*VH/(MH*MH-ms2))**2*unit_width


def population(q, energy, v, lam, kappa, unit_width, source_gamma=1.):
    """Conditional: all net charge in these objects, release as free s, no
    subsequent inverse reactions/accretion, common halo, cold asymptotic DM.
    Returns late number/rest-mass accounting; not a formation calculation.
    """
    ms = np.sqrt(2*lam)*v
    mx, md = G*v, np.sqrt(G**2+HD**2)*v
    gx, gd = hidden(v, lam)
    gv = visible_coefficient(v, lam, unit_width)*kappa*kappa
    total = gx+gd+gv
    bx, bd, bv = gx/total, gd/total, gv/total
    release = G*v*q-energy
    ns = release/(ms*source_gamma)
    nx, nd = 2*bx*ns, 2*bd*ns
    late_mass = energy+nx*mx+nd*md
    yq = 4.4e-10*q/late_mass
    return dict(source_release_GeV=release, source_scalar_number_per_object=ns,
                hidden_X_branch=bx, hidden_D_branch=bd, visible_branch=bv,
                total_width_GeV=total, scalar_lifetime_seconds=HBAR/total,
                late_X_particles_per_object=nx, late_D_particles_per_object=nd,
                late_mass_per_object_GeV=late_mass,
                object_halo_fraction=energy/late_mass,
                required_charge_yield=yq,
                old_fixed_charge_total_DM_ratio=4.8294009076988175e-9/q*late_mass/4.4e-10,
                scalar_release_GeV_per_entropy=yq/q*release,
                scalar_visible_energy_GeV_per_entropy=yq/q*release*bv,
                charge_yield_residual=(yq/q*q-yq)/yq,
                rest_mass_normalization_residual=yq/q*late_mass/4.4e-10-1.)
