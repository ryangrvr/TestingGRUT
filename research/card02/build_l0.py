"""Disclosed arithmetic transcription builder; no candidate scoring or search."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def app(f, *args):
    for x in args:
        f = ["apply", f, x]
    return f


def lam(names, body):
    for name in reversed(names):
        body = ["lambda", name, body]
    return body


def call(name, *args):
    return app(["use", name], *args)


def summation(name, size, body):
    return ["sum_finite", ["range", size], ["lambda", name, body]]


def at(m, *indices):
    return app(m, *indices)


def build():
    defs = []
    def add(name, args, body):
        defs.append([name, lam(args, body)])
    def product_part(part):
        if part == "r":
            body = ["-", ["*", at("ar", "u", "k"), at("br", "k", "v")],
                         ["*", at("ai", "u", "k"), at("bi", "k", "v")]]
        else:
            body = ["+", ["*", at("ar", "u", "k"), at("bi", "k", "v")],
                         ["*", at("ai", "u", "k"), at("br", "k", "v")]]
        return summation("k", "d", body)
    for part in ("r", "i"):
        add("mul_"+part, ["ar", "ai", "br", "bi", "d", "u", "v"], product_part(part))
    expectation = summation("u", "d", summation("v", "d",
        ["-", ["*", at("rr", "u", "v"), at("ar", "v", "u")],
              ["*", at("ri", "u", "v"), at("ai", "v", "u")]]))
    add("expect_r", ["ar", "ai", "rr", "ri", "d"], expectation)
    qi_r = lam(["u", "v"], at("qr", "i", "a", "u", "v"))
    qi_i = lam(["u", "v"], at("qi", "i", "a", "u", "v"))
    qj_r = lam(["u", "v"], at("qr", "j", "b", "u", "v"))
    qj_i = lam(["u", "v"], at("qi", "j", "b", "u", "v"))
    common = ["qr", "qi", "rr", "ri", "d"]
    add("mean", common+["i", "a"], call("expect_r", qi_r, qi_i, "rr", "ri", "d"))
    products = [lam(["u", "v"], call("mul_"+part, qi_r, qi_i, qj_r, qj_i, "d", "u", "v"))
                for part in ("r", "i")]
    mi = call("mean", *common, "i", "a")
    mj = call("mean", *common, "j", "b")
    add("cov", common+["i", "j", "a", "b"],
        ["-", call("expect_r", *products, "rr", "ri", "d"), ["*", mi, mj]])
    for part in ("r", "i"):
        q = "qr" if part == "r" else "qi"
        bracket = ["-", ["-", call("mul_"+part, qi_r, qi_i, qj_r, qj_i, "d", "u", "v"),
            ["*", mj, at(q, "i", "a", "u", "v")]],
            ["*", mi, at(q, "j", "b", "u", "v")]]
        piece = ["*", call("cov", *common, "i", "j", "a", "b"), bracket]
        interaction = summation("i", "n", summation("j", "n",
            ["ite", ["<", "i", "j"], summation("a", ["nat", 3],
              summation("b", ["nat", 3], piece)), "0"]))
        local = summation("i", "n", summation("a", ["nat", 3],
            ["*", at("alpha", "i", "a"), at(q, "i", "a", "u", "v")]))
        add("H_"+part, common+["n", "kappa", "alpha", "u", "v"],
            ["+", ["*", "kappa", interaction], local])
    h = [lam(["u", "v"], call("H_"+part, *common, "n", "kappa", "alpha", "u", "v"))
         for part in ("r", "i")]
    comm = []
    for part in ("r", "i"):
        comm.append(["-", call("mul_"+part, *h, "rr", "ri", "d", "u", "v"),
                          call("mul_"+part, "rr", "ri", *h, "d", "u", "v")])
    equation = ["and", ["=", at("drr", "u", "v"), comm[1]],
        ["=", at("dri", "u", "v"), ["-", "0", comm[0]]]]
    predicate = ["forall", "u", ["forall", "v", ["->",
        ["and", ["in", "u", ["range", "d"]], ["in", "v", ["range", "d"]]], equation]]]
    expr = lam(common+["n", "kappa", "alpha", "drr", "dri"], predicate)
    return {"expression": expr, "definitions": defs,
        "scope": "Explicit finite matrix derivative equations only; quantum/time/domain imports priced separately",
        "operator_input": "qr/qi are curried Pauli tensor arrays with supplied elementary factorization",
        "decision_status": "No complete continuum/chain decision certificate"}


if __name__ == "__main__":
    (ROOT/"L0_CORE.json").write_text(json.dumps(build(), separators=(",", ":"))+"\n")
