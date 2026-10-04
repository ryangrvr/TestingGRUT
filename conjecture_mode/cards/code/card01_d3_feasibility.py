#!/usr/bin/env python3
"""CARD-01 D3 interface feasibility + C1 early-time check.  MODEL ONLY: no likelihood, no data file,
no DESI product.  Runs CAMB (pip camb 2.0.4) with the exact Card #1 v1 w(a) supplied as a table
(DarkEnergyPPF.set_w_a_table), NOT replaced by its CPL projection.

Checks:
  F1  CAMB background H(z) with the tabulated w(a) agrees with an independent quadrature of
      rho_DE(z)/rho_DE0 = exp(3 eps int_0^z R(z')/(1+z') dz') at the same physical parameters;
  F2  eps = 0 reproduces CAMB LambdaCDM;
  F3  CMB TT/EE spectra compute without error on the primary branch (perturbations via PPF; the
      model never touches w = -1 for eps < 0, so no divide-crossing pathology);
  E1  early-time DE fraction Omega_DE(z) at z = 1100 and z = 1e5 for primary and control branches.
Physical parameters used here are a FIXED ILLUSTRATIVE POINT (Planck-like round numbers), not fitted
to anything; they only exercise the code path.
"""
import numpy as np
import camb
from scipy.integrate import quad

OMR, OLR = 0.31, 0.69          # E_ref (frozen historical reference background, NOT the varied cosmology)


def R_of_z(z):
    e2 = OMR * (1 + z) ** 3 + OLR
    return e2 / (1 + e2)


def w_of_a(a, eps):
    return -1.0 + eps * R_of_z(1.0 / a - 1.0)


def table(eps, n=3000):
    a = np.logspace(-9, 0, n)
    return a, w_of_a(a, eps)


def params(eps, lmax=2500):
    p = camb.CAMBparams()
    p.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12, mnu=0.06, omk=0, tau=0.055)
    p.InitPower.set_params(As=2.1e-9, ns=0.965)
    p.set_for_lmax(lmax, lens_potential_accuracy=0)
    if eps != 0.0:
        de = camb.dark_energy.DarkEnergyPPF()
        a, w = table(eps)
        de.set_w_a_table(a, w)
        p.DarkEnergy = de
    return p


def rho_ratio(z, eps):
    I = quad(lambda zz: R_of_z(zz) / (1 + zz), 0, z, limit=400)[0]
    return np.exp(3 * eps * I)


def main():
    zs = np.array([0.0, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.33, 3.0, 10.0, 1100.0])
    print("CARD-01 D3 feasibility (model only).  camb", camb.__version__)
    for eps in [0.0, -0.05, -0.2, -0.5, +0.2]:
        tag = "PRIMARY" if eps < 0 else ("NULL" if eps == 0 else "CONTROL - NOT A GRUT CLAIM")
        p = params(eps)
        bg = camb.get_background(p)
        H = np.array([bg.hubble_parameter(z) for z in zs])
        # independent background: same physical densities, rho_DE from quadrature
        H0 = p.H0
        rho = bg.get_background_redshift_evolution if False else None
        dens0 = bg.get_background_densities(1.0, vars=['tot', 'de'])
        frac_de0 = dens0['de'][0] / dens0['tot'][0]
        Hind = []
        for z in zs:
            a = 1 / (1 + z)
            d = bg.get_background_densities(a, vars=['tot', 'de'])
            nonde = d['tot'][0] - d['de'][0]                     # a^4 * 8piG rho, non-DE part from CAMB
            de = dens0['de'][0] * rho_ratio(z, eps) * a ** 4     # a^4 * 8piG rho_DE from independent quadrature
            Hind.append(H0 * np.sqrt((nonde + de) / dens0['tot'][0]) / a ** 2)
        Hind = np.array(Hind)
        dev = np.max(np.abs(H / Hind - 1))
        d1100 = bg.get_background_densities(1 / 1101, vars=['tot', 'de'])
        d1e5 = bg.get_background_densities(1e-5, vars=['tot', 'de'])
        print(f"\n eps={eps:+.2f} [{tag}]  Omega_DE0={frac_de0:.5f}")
        print(f"   F1 max|H_CAMB/H_indep - 1| over z in {list(zs)} = {dev:.2e}")
        print(f"   E1 Omega_DE(z=1100) = {d1100['de'][0]/d1100['tot'][0]:.3e}   Omega_DE(z=1e5) = {d1e5['de'][0]/d1e5['tot'][0]:.3e}")
        print(f"   w(z=0)={w_of_a(1.0, eps):+.6f}  w(z=2)={w_of_a(1/3, eps):+.6f}  w(z=1100)={w_of_a(1/1101, eps):+.6f}")
        if eps in (0.0, -0.2, +0.2):
            res = camb.get_results(p)
            cl = res.get_cmb_power_spectra(p, CMB_unit='muK')['total']
            print(f"   F3 CMB spectra OK: D_l^TT(l=220)={cl[220,0]:.2f}  D_l^EE(l=1000)={cl[1000,1]:.4f}  "
                  f"theta*={res.get_derived_params()['thetastar']:.6f}  rdrag={res.get_derived_params()['rdrag']:.4f}")
    # F2: eps=0 vs pure LCDM path, and tiny eps continuity
    p0 = params(0.0); pm = params(-1e-6)
    r0 = camb.get_background(p0); rm = camb.get_background(pm)
    print(f"\n F2 continuity eps=-1e-6 vs LCDM: max|dH/H| = "
          f"{max(abs(rm.hubble_parameter(z)/r0.hubble_parameter(z)-1) for z in zs):.2e}")


if __name__ == "__main__":
    main()
