# K1-HS final audit: Linder's separately stated sound-speed condition reduces to the Bellini-Sawicki/Peirone N.
# Units: rho/H^2 in m_p^2 units (Linder), R = m_p^2/M_*^2; ' = d/dln a.
import sympy as sp
aB, aM, aBp, R, h, rm, rde, w = sp.symbols('alpha_B alpha_M alpha_Bp R hprime rho_m_over_H2 rho_de_over_H2 w')
A = aB + 2*aM; D = (2 - aB)*A + 2*aBp
linder = (1 - aB/2)*A + (aBp + aB*h) + rm*(1 - R) + rde*(1 + w)          # (H aB)'/H = aB' + aB H'/H, h = H'/H
fried = sp.solve(sp.Eq(-2*h, rm + rde*(1 + w)), rde)[0]                     # -2H'/H = (rho_m + rho_de(1+w))/(m_p^2 H^2)
N_BS = D/2 - (2 - aB)*h - rm*R                                             # rho_m/(H^2 M_*^2) = R * rho_m/(H^2 m_p^2)
print("Linder alpha c_s^2 - (D/2 + aB H'/H + rho_m(1-R)/H^2 + rho_de(1+w)/H^2) =", sp.simplify(linder - (D/2 + aB*h + rm*(1 - R) + rde*(1 + w))))
print("Linder alpha c_s^2 (Friedmann-substituted) - N_BS =", sp.simplify(linder.subs(rde, fried) - N_BS))
print("=> D/2 is NOT alpha c_s^2 in general; the 'reading L' identification was mistaken; N is the unique stability combination")
