"""production_map / b2: the spectator-condensate bath as a hidden bath (does it restore Z2?), per H_inf.

Output: results/b2_spectator_bath.json.  Needs results/b1_spectator.json and a1 (imports a1 helpers).  < 10 s.

For instant reheating and minimal couplings (b1): eps = c_eps (H/M_P)^2, Y_N in [Y_floor, Y_J].
The gas part of an overpopulated bath gives Pi >= Pi_th(T_E) (a1 section 7), so the restored phase lasts at least
until T_SM = T_c(1)/xi_E [D].  With Pi = A H in RD, the zero mode responds iff A > 2 sqrt(lphi) v [D].
"""
import math

from common import v, lphi, MP, Y_MAX_STRICT, COUP_PHI, gamma_el, H_of_T, load, dump

b1 = load("b1_spectator.json")
a1 = load("a1_restoration.json")
inst = b1["instant reheating, minimal couplings [C]"]
c_eps = inst["eps coef from sampled initial data incl. slow roll [C]"]
cxi = inst["xi_E = coef (H/M_P)^(1/2) [C]"]
cF = inst["Y_floor = coef (H/M_P)^(3/2) (mu=0 number at the same energy) [C]"]
cJ = inst["Y_J = coef (H/M_P)^(3/2) (condensate quanta conserved) [C]"]
row1 = a1["Tc_of_z [C] (Hartree, m_phi=0 at boundary, Bose, thermal masses)"][0]
Tc1 = row1["Tc_over_v"] * v
TRH_COEF = (90 / (math.pi ** 2 * 106.75)) ** 0.25

rows = []
for H in (2e3, 1e4, 1e5, 1e6, 1e8, 1e9, 1e10, 3e10, 1e11, 2e11, 1e12, 1e13):
    xiE = cxi * math.sqrt(H / MP)
    Tend = Tc1 / xiE
    Hend = H_of_T(Tend)
    A = lphi * v ** 2 / Hend
    nphi = row1["n_phi_over_T3"] * Tc1 ** 3
    E = row1["E_mean_over_T"] * Tc1
    G = gamma_el(nphi, E, COUP_PHI) / Hend
    rows.append(dict(H_inf=H, eps=c_eps * (H / MP) ** 2, xi_E=xiE, Y_floor=cF * (H / MP) ** 1.5,
                     Y_J=cJ * (H / MP) ** 1.5, Y_floor_over_Ymax=cF * (H / MP) ** 1.5 / Y_MAX_STRICT,
                     Y_J_over_Ymax=cJ * (H / MP) ** 1.5 / Y_MAX_STRICT,
                     T_RH_instant=TRH_COEF * math.sqrt(H * MP),
                     T_SM_end_gas_bound=Tend, A_gas_GeV=A, zero_mode_responds=bool(A > 2 * math.sqrt(lphi) * v),
                     Gamma_phi_el_over_H_at_exit_if_thermal=G))
OUT = {"module": "production_map", "part": "B2 spectator bath as a hidden bath",
       "rows [D, C]": rows,
       "H_inf above which a mu=0-thermalized spectator bath is kinetically thermal at its exit (xi_E >= 6.6e-5)":
           MP * (6.6e-5 / cxi) ** 2,
       "H_inf above which the gas-bound zero mode responds (A > 200 GeV)":
           MP * (math.sqrt(2 * math.sqrt(lphi) * v / (lphi * v ** 2 / H_of_T(Tc1))) / cxi) ** 2,
       "note": "T_SM_end_gas_bound uses only the gas part (mu = 0 at the energy temperature); the condensate excess "
               "(J/n_eq > 1) prolongs restoration. Walls at the exit form with p = 1/2 by the P symmetry of P_eq; "
               "inflation_spectator already gives p_- = 1/2 for H_inf >~ 1.4-1.8 TeV [Q]."}
dump("b2_spectator_bath.json", OUT)
if __name__ == "__main__":
    for r in rows:
        print({k: ("%.3g" % val if isinstance(val, float) else val) for k, val in r.items()})
    print(OUT["H_inf above which a mu=0-thermalized spectator bath is kinetically thermal at its exit (xi_E >= 6.6e-5)"],
          OUT["H_inf above which the gas-bound zero mode responds (A > 200 GeV)"])
