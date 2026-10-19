#!/usr/bin/env python3
"""Exact one-loop scalar tensor/RG-identity check, using Fraction polynomials.

No third-party symbolic package, fitting, RG integration, or soliton determinant.
See the frozen rg_protocol.json for scope and normalization.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import sys
from fractions import Fraction as F
from pathlib import Path


NAMES = ("phi", "x", "y", "D", "C", "mp", "mx", "md", "lp", "lx", "ld", "c", "d", "alpha")
INDEX = {name: i for i, name in enumerate(NAMES)}
ZERO_POWER = (0,) * len(NAMES)
HERE = Path(__file__).resolve().parent


class Poly:
    def __init__(self, terms=None):
        self.terms = {p: F(c) for p, c in (terms or {}).items() if c}

    @classmethod
    def scalar(cls, value):
        value = F(value)
        return cls({ZERO_POWER: value}) if value else cls()

    @classmethod
    def symbol(cls, name):
        powers = list(ZERO_POWER)
        powers[INDEX[name]] = 1
        return cls({tuple(powers): F(1)})

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Poly) else Poly.scalar(other)

    def __add__(self, other):
        other = self.coerce(other)
        result = self.terms.copy()
        for powers, coefficient in other.terms.items():
            result[powers] = result.get(powers, F(0)) + coefficient
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({p: -c for p, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        result = {}
        for pa, ca in self.terms.items():
            for pb, cb in other.terms.items():
                powers = tuple(a + b for a, b in zip(pa, pb))
                result[powers] = result.get(powers, F(0)) + ca * cb
        return Poly(result)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (F(1) / F(scalar))

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("Polynomial powers must be nonnegative integers.")
        result = Poly.scalar(1)
        for _ in range(power):
            result = result * self
        return result

    def diff(self, name):
        i = INDEX[name]
        result = {}
        for powers, coefficient in self.terms.items():
            if powers[i]:
                new = list(powers)
                new[i] -= 1
                result[tuple(new)] = coefficient * powers[i]
        return Poly(result)

    def substitute(self, replacements):
        replacements = {INDEX[k]: F(v) for k, v in replacements.items()}
        result = {}
        for powers, coefficient in self.terms.items():
            new = list(powers)
            for i, value in replacements.items():
                coefficient *= value ** powers[i]
                new[i] = 0
            powers_new = tuple(new)
            result[powers_new] = result.get(powers_new, F(0)) + coefficient
        return Poly(result)

    def swap(self, a, b):
        ia, ib = INDEX[a], INDEX[b]
        result = {}
        for powers, coefficient in self.terms.items():
            new = list(powers)
            new[ia], new[ib] = new[ib], new[ia]
            result[tuple(new)] = coefficient
        return Poly(result)

    def constant(self):
        if any(p != ZERO_POWER for p in self.terms):
            raise ValueError(f"Not a constant: {self}")
        return self.terms.get(ZERO_POWER, F(0))

    def is_zero(self):
        return not self.terms

    def json_terms(self):
        return [{"coefficient": str(c), "powers": {n: p[i] for i, n in enumerate(NAMES) if p[i]}}
                for p, c in sorted(self.terms.items(), reverse=True)]

    def __str__(self):
        if not self.terms:
            return "0"
        pieces = []
        for powers, coefficient in sorted(self.terms.items(), reverse=True):
            factors = [name if powers[i] == 1 else f"{name}^{powers[i]}"
                       for i, name in enumerate(NAMES) if powers[i]]
            mag = abs(coefficient)
            value = "*".join(([str(mag)] if mag != 1 or not factors else []) + factors)
            pieces.append(("-" if coefficient < 0 else "+", value))
        first_sign, first = pieces[0]
        text = ("-" if first_sign == "-" else "") + first
        return text + "".join(f" {sign} {piece}" for sign, piece in pieces[1:])


Z = Poly()
P = {name: Poly.symbol(name) for name in NAMES}


def assert_equal(a, b, label, checks):
    residual = Poly.coerce(a) - Poly.coerce(b)
    if not residual.is_zero():
        raise AssertionError(f"{label}: nonzero residual {residual}")
    checks.append({"check": label, "passed": True, "residual_monomials": 0})


def tensor_from_couplings(couplings):
    """Explicit 1/4! tensor seeds, independent of Hessian differentiation."""
    lp, lx, ld, c, d, alpha = (couplings[k] for k in ("lp", "lx", "ld", "c", "d", "alpha"))
    seeds = [((0, 0, 0, 0), 6 * lp), ((1, 1, 1, 1), 6 * lx),
             ((2, 2, 2, 2), 6 * lx), ((3, 3, 3, 3), 6 * ld),
             ((0, 0, 1, 1), 2 * c), ((0, 0, 2, 2), 2 * c),
             ((0, 0, 3, 3), 2 * d), ((1, 1, 2, 2), 2 * lx),
             ((1, 1, 3, 3), 2 * alpha), ((2, 2, 3, 3), 2 * alpha)]
    tensor = {}
    for indices, value in seeds:
        for permuted in set(itertools.permutations(indices)):
            tensor[permuted] = value
    return tensor


def main():
    protocol_path = HERE / "rg_protocol.json"
    protocol_bytes = protocol_path.read_bytes()
    protocol = json.loads(protocol_bytes)
    output = HERE / "rg_results.json"
    if output.exists():
        raise SystemExit("RG results already exist; refusing to overwrite the recorded check.")
    checks = []
    phi, x, y, D = (P[k] for k in ("phi", "x", "y", "D"))
    C, mp, mx, md, lp, lx, ld, c, d, alpha = (P[k] for k in ("C", "mp", "mx", "md", "lp", "lx", "ld", "c", "d", "alpha"))
    r2 = x * x + y * y

    # Arithmetic/differentiation preflight, before the physics contractions.
    assert_equal((x + y) * (x - y), x*x - y*y, "polynomial multiplication preflight", checks)
    assert_equal(((x + y)**3).diff("x"), 3 * (x + y)**2, "polynomial differentiation preflight", checks)

    # Route A: explicit symmetric tensor contractions.
    tensor = tensor_from_couplings(P)
    beta_tensor = {}
    for ia, ib, ic, id_ in itertools.product(range(4), repeat=4):
        value = Z
        for ie, iff in itertools.product(range(4), repeat=2):
            value += tensor.get((ia, ib, ie, iff), Z) * tensor.get((ie, iff, ic, id_), Z)
            value += tensor.get((ia, ic, ie, iff), Z) * tensor.get((ie, iff, ib, id_), Z)
            value += tensor.get((ia, id_, ie, iff), Z) * tensor.get((ie, iff, ib, ic), Z)
        beta_tensor[(ia, ib, ic, id_)] = value
    B = {
        "lp": beta_tensor[(0, 0, 0, 0)] / 6,
        "lx": beta_tensor[(1, 1, 1, 1)] / 6,
        "ld": beta_tensor[(3, 3, 3, 3)] / 6,
        "c": beta_tensor[(0, 0, 1, 1)] / 2,
        "d": beta_tensor[(0, 0, 3, 3)] / 2,
        "alpha": beta_tensor[(1, 1, 3, 3)] / 2,
    }
    reconstructed_beta = tensor_from_couplings(B)
    for indices in itertools.product(range(4), repeat=4):
        assert_equal(beta_tensor[indices], reconstructed_beta.get(indices, Z), f"quartic closure {indices}", checks)

    mass = (mp, mx, mx, md)
    beta_mass = {}
    for ia, ib in itertools.product(range(4), repeat=2):
        beta_mass[(ia, ib)] = sum((tensor.get((ia, ib, j, j), Z) * mass[j] for j in range(4)), Z)
    B.update(mp=beta_mass[(0, 0)], mx=beta_mass[(1, 1)], md=beta_mass[(3, 3)],
             C=sum((m*m for m in mass), Z) / 2)
    for ia, ib in itertools.product(range(4), repeat=2):
        expected = (B["mp"], B["mx"], B["mx"], B["md"])[ia] if ia == ib else Z
        assert_equal(beta_mass[(ia, ib)], expected, f"mass closure {(ia, ib)}", checks)

    expected = {
        "lp": 18*lp**2 + 4*c**2 + 2*d**2,
        "lx": 20*lx**2 + 2*c**2 + 2*alpha**2,
        "ld": 18*ld**2 + 2*d**2 + 4*alpha**2,
        "c": c*(6*lp + 8*lx + 8*c) + 2*d*alpha,
        "d": d*(6*lp + 6*ld + 8*d) + 4*c*alpha,
        "alpha": 2*c*d + alpha*(8*lx + 6*ld + 8*alpha),
        "mp": 6*lp*mp + 4*c*mx + 2*d*md,
        "mx": 2*c*mp + 8*lx*mx + 2*alpha*md,
        "md": 2*d*mp + 4*alpha*mx + 6*ld*md,
        "C": (mp**2 + 2*mx**2 + md**2) / 2,
    }
    for name in expected:
        assert_equal(B[name], expected[name], f"projected beta formula {name}", checks)

    # Route B: a directly entered potential, differentiated to its Hessian.
    V0 = C + mp*phi**2/2 + mx*r2/2 + md*D**2/2
    V0 += lp*phi**4/4 + lx*r2**2/4 + ld*D**4/4
    V0 += c*phi**2*r2/2 + d*phi**2*D**2/2 + alpha*r2*D**2/2
    fields = ("phi", "x", "y", "D")
    H = [[V0.diff(ia).diff(ib) for ib in fields] for ia in fields]
    trace_H_squared = sum((H[i][j]*H[j][i] for i in range(4) for j in range(4)), Z)
    rg_left = sum((B[name] * V0.diff(name) for name in B), Z)
    assert_equal(rg_left, trace_H_squared/2, "independent full polynomial RG identity including C", checks)
    assert_equal(V0.swap("x", "y"), V0, "potential x-y interchange", checks)
    assert_equal(trace_H_squared.swap("x", "y"), trace_H_squared, "Hessian trace x-y interchange", checks)
    assert_equal(x*V0.diff("y") - y*V0.diff("x"), Z, "potential infinitesimal O(2) invariance", checks)
    assert_equal(x*trace_H_squared.diff("y") - y*trace_H_squared.diff("x"), Z,
                 "Hessian trace infinitesimal O(2) invariance", checks)

    # Tree-vacuum running identity, kept polynomial by multiplying by lp^2.
    tree_vac_running_times_lp2 = lp**2*B["C"] - mp*lp*B["mp"]/2 + mp**2*B["lp"]/4
    vacuum_trace_times_lp2 = ((-2*mp*lp)**2 + 2*(mx*lp-c*mp)**2 + (md*lp-d*mp)**2)/2
    assert_equal(tree_vac_running_times_lp2, vacuum_trace_times_lp2,
                 "tree minimum running equals half vacuum Hessian trace", checks)

    omission_residuals = {}
    for name, replacements in (("alpha", {"alpha": 0}), ("mx", {"alpha": 0, "mx": 0}), ("C", {})):
        residual = (rg_left - B[name]*V0.diff(name) - trace_H_squared/2).substitute(replacements)
        if residual.is_zero():
            raise AssertionError(f"Omission check unexpectedly vanishes for {name}")
        omission_residuals[name] = {"residual": str(residual), "terms": residual.json_terms()}

    match = {"lp": F(1,100), "lx": F(1,1000), "ld": F(1,1000), "c": F(1,100),
             "d": F(1,10000), "alpha": F(0), "mp": F(-1,100), "mx": F(0),
             "md": F(1,100), "C": F(1,400)}
    match_rows = {}
    for name, expression in B.items():
        numerator = expression.substitute(match).constant()
        match_rows[name] = {"B_exact": str(numerator), "B_float": float(numerator),
                            "beta": float(numerator)/(16*math.pi**2)}
    leading_logs = []
    for factor in (0.5, 1.0, 2.0):
        t = math.log(factor)
        leading_logs.append({"mu_over_mu0": factor, "log_ratio": t,
                             "parameters_v1": {name: float(value)+match_rows[name]["beta"]*t
                                                for name,value in match.items()}})
    Bvac_exact = tree_vac_running_times_lp2.substitute(match).constant() / match["lp"]**2
    vacuum_hessian = {"phi": -2*match["mp"], "x": match["mx"]-match["c"]*match["mp"]/match["lp"],
                      "y": match["mx"]-match["c"]*match["mp"]/match["lp"],
                      "D": match["md"]-match["d"]*match["mp"]/match["lp"]}
    vacuum = {"tree_minimum_value_v4": str(match["C"]-match["mp"]**2/(4*match["lp"])),
              "vacuum_Hessian_eigenvalues_over_v2": {k:str(v) for k,v in vacuum_hessian.items()},
              "B_tree_minimum_value_over_v4_exact": str(Bvac_exact),
              "beta_tree_minimum_value_over_v4": float(Bvac_exact)/(16*math.pi**2)}
    result = {
        "protocol_sha256": hashlib.sha256(protocol_bytes).hexdigest(),
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": sys.version, "platform": platform.platform(),
        "arithmetic": "exact fractions over sparse multivariate polynomials",
        "normalization": "B_z=16*pi^2*d z/d ln(mu), mp/mx/md are mass-squared coefficients",
        "checks_passed": len(checks), "checks": checks,
        "full_identity_residual_monomials": len((rg_left-trace_H_squared/2).terms),
        "identity_lhs_monomials": len(rg_left.terms), "identity_rhs_monomials": len((trace_H_squared/2).terms),
        "beta_B_expressions": {name:str(expression) for name,expression in B.items()},
        "matching_v1": {k:str(v) for k,v in match.items()}, "matching_beta": match_rows,
        "leading_logs_first_order_only": leading_logs, "vacuum_normalization": vacuum,
        "incorrect_omission_residuals": omission_residuals,
        "limits": protocol["limits"],
    }
    polynomials = {"symbols": NAMES, "V0": V0.json_terms(),
                   "Hessian": [[entry.json_terms() for entry in row] for row in H],
                   "trace_H_squared": trace_H_squared.json_terms(),
                   "RG_left_B": rg_left.json_terms(),
                   "RG_residual": (rg_left-trace_H_squared/2).json_terms(),
                   "beta_B": {name:expr.json_terms() for name,expr in B.items()}}
    output.write_text(json.dumps(result, indent=2)+"\n")
    (HERE/"rg_polynomials.json").write_text(json.dumps(polynomials, indent=2)+"\n")
    lines = ["B_z = 16*pi^2 * beta_z", *(f"B_{name} = {expr}" for name,expr in B.items()),
             f"Exact checks passed: {len(checks)}", "Full polynomial RG residual: 0",
             f"Identity monomials per side: {len(rg_left.terms)}"]
    (HERE/"rg_beta_functions.txt").write_text("\n".join(lines)+"\n")
    (HERE/"rg_execution.log").write_text("\n".join(lines)+"\n")
    print(json.dumps({"checks_passed":len(checks), "identity_residual_monomials":0,
                      "beta_B_expressions":result["beta_B_expressions"],
                      "matching_beta":match_rows,"vacuum_normalization":vacuum},indent=2))


if __name__ == "__main__":
    main()
