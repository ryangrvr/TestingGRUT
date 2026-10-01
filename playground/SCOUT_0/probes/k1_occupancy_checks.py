# K1 occupancy algebra (sympy). Literature inputs are IMPORTED-SECONDARY (primary arXiv text not reachable from
# this sandbox): Only Run Gravity (Linder 2020, alpha_B = 0): G_matter = R(1+a_M), G_light = R(1+a_M/2), R = m_p^2/M_*^2;
# small-alpha Horndeski (alpha_T = 0, M_* ~ m_p, quasi-static): (Sigma-1)/(mu-1) ~ (1+r)/(2+r), r = alpha_B/alpha_M.
import sympy as sp
aM, R, r, x, al, eps = sp.symbols('alpha_M R r x alpha epsilon')
mu, Sig = R*(1+aM), R*(1+aM/2)
off = sp.simplify((Sig-1) - (mu-1)/2)
print("Only Run: Sigma-1 - (mu-1)/2 =", off, " -> with R = 1+eps:", sp.expand(off.subs(R, 1+eps)))
print("Only Run: mu*eta = 2 Sigma - mu =", sp.simplify(2*Sig - mu), " (GRUT family: mu*eta = 1)")
ratio = (1+r)/(2+r)
print("Horndeski small-alpha ratio limits: f(R) r=-1 ->", ratio.subs(r, -1), "(Sigma=1);  alpha_M->0 (r->oo) ->",
      sp.limit(ratio, r, sp.oo), "(no slip, Sigma=mu);  alpha_B=0 (r=0) ->", ratio.subs(r, 0), "(GRUT line)")
print("GRUT family: Sigma =", sp.simplify((1+x*al)*(1+1/(1+x*al))/2), "; ratio (Sigma-1)/(mu-1) =",
      sp.simplify(((1+x*al)*(1+1/(1+x*al))/2 - 1)/(x*al)), "for all x, alpha")
print("Staying on the line at all a in Only Run needs eps(a) = 0 for all a => M_* constant => alpha_M = dln M_*^2/dln a = 0 => GR")
