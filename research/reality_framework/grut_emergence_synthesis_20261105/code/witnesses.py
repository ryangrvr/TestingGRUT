"""Finite illustrations with known analytic answers; not cosmological tests."""
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
checks = []


def check(name, value, expected, atol=1e-11):
    error = float(np.max(np.abs(np.asarray(value) - np.asarray(expected))))
    checks.append({"name": name, "max_absolute_error": error,
                   "tolerance": atol, "pass": error <= atol})


def rho(psi):
    return np.outer(psi, psi.conj())


def reduce(r, dims, keep):
    a = r.reshape(tuple(dims) * 2)
    live = list(dims)
    for i in reversed(range(len(dims))):
        if i not in keep:
            a = np.trace(a, axis1=i, axis2=i + len(live))
            live.pop(i)
    n = int(np.prod(live))
    return a.reshape(n, n)


def entropy(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


def distance(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b))) / 2)


def cnot(control, target, n=3):
    u = np.zeros((2**n, 2**n), complex)
    for j in range(2**n):
        k = j ^ (1 << (n-1-target)) if (j >> (n-1-control)) & 1 else j
        u[k, j] = 1
    return u


def whole_parts():
    h = np.kron(X, X)
    plus = np.array([1, 0, 0, 1j]) / np.sqrt(2)
    minus = plus.conj()
    a, b = rho(plus), rho(minus)
    for i in [0, 1]:
        check(f"same initial marginal {i}", reduce(a, [2,2], [i]), I/2)
        check(f"same other marginal {i}", reduce(b, [2,2], [i]), I/2)
    check("same mean energy", np.trace(h @ a), np.trace(h @ b))
    t = np.pi/8
    u = expm(-1j*h*t)
    at, bt = u@a@u.conj().T, u@b@u.conj().T
    d = distance(reduce(at,[2,2],[0]), reduce(bt,[2,2],[0]))
    check("different local futures", d, np.sin(2*t))
    check("unitarity", u.conj().T @ u, np.eye(4))
    check("global purity", np.trace(at@at), 1)
    check("reverse recovers preparation", u.conj().T@at@u, a)
    shifted = expm(-1j*(h + 3.75*np.eye(4))*t)
    check("constant energy shift invisible here", shifted@a@shifted.conj().T, at)
    p = rho(np.array([1,0,0,0], complex))
    ent = [entropy(reduce(expm(-1j*h*s)@p@expm(1j*h*s),[2,2],[0]))
           for s in [0, np.pi/4, np.pi/2]]
    check("local entropy can rise and fall", ent, [0,np.log(2),0])
    return {"local_trace_distance":d, "local_entropy_cycle_nats":ent}


def records():
    d, nt = 8, 8
    gates = [cnot(0,1), cnot(0,2), np.eye(d), np.eye(d),
             cnot(0,2), cnot(0,1), np.eye(d)]
    ws = [np.eye(d, dtype=complex)]
    for u in gates:
        ws.append(u @ ws[-1])
    v = np.vstack(ws) / np.sqrt(nt)
    c = np.zeros((d*nt,d*nt), complex)
    for t,u in enumerate(gates):
        sl, sr = slice(t*d,(t+1)*d), slice((t+1)*d,(t+2)*d)
        c[sl,sl] += np.eye(d)/2
        c[sr,sr] += np.eye(d)/2
        c[sr,sl] -= u/2
        c[sl,sr] -= u.conj().T/2
    check("history isometry", v.conj().T@v, np.eye(d))
    check("constraint annihilates all histories", c@v, np.zeros_like(v))
    check("physical-space dimension", np.sum(np.abs(np.linalg.eigvalsh(c))<1e-10), d)
    check("constraint positive", min(np.linalg.eigvalsh(c)), 0)
    psi0 = np.zeros(d,complex)
    psi0[0],psi0[4] = np.sqrt(.3),np.sqrt(.7)
    hist = v@psi0
    check("history norm", np.vdot(hist,hist), 1)
    maximix = v@v.conj().T/d
    rr, mi = [], []
    for t,w in enumerate(ws):
        s = slice(t*d,(t+1)*d)
        cond = rho(hist[s])*nt
        check(f"conditional history {t}", cond, rho(w@psi0))
        check(f"maximally mixed conditional {t}", maximix[s,s]*nt, np.eye(d)/d)
        rs = reduce(cond,[2,2,2],[0])
        rec = reduce(cond,[2,2,2],[1])
        joint = reduce(cond,[2,2,2],[0,1])
        rr.append(entropy(rec))
        mi.append(entropy(rs)+entropy(rec)-entropy(joint))
    hb = -.3*np.log(.3)-.7*np.log(.7)
    check("blank-record entropy", rr[0], 0)
    check("written-record entropy", rr[1], hb)
    check("copied-record quantum MI", mi[2], hb)
    check("coherent pre-copy quantum MI", mi[1], 2*hb)
    check("erased record", rr[7], 0)
    return {"record_entropy_nats":rr, "system_record_quantum_MI_nats":mi,
            "physical_dimension":d, "clock_readings":nt}


def entropy_identity():
    beta = 1.7
    hb = np.diag([0.,1.3])
    tau = expm(-beta*hb); tau /= np.trace(tau)
    ra = np.diag([.12,.88])
    initial = np.kron(ra,tau)
    h = np.kron(.45*Z,I)+np.kron(I,hb)+.37*(np.kron(X,X)+.2*np.kron(Z,Y))
    rows = []
    for t in [.7,2.1]:
        u = expm(-1j*h*t); final = u@initial@u.conj().T
        a,b = reduce(final,[2,2],[0]),reduce(final,[2,2],[1])
        ds = entropy(a)-entropy(ra)
        de = float(np.trace(hb@(b-tau)).real)
        mi = entropy(a)+entropy(b)-entropy(final)
        rel = -entropy(b)-float(np.trace(b@np.diag(np.log(np.diag(tau)))).real)
        check(f"entropy identity t={t}", ds+beta*de, mi+rel)
        check(f"total entropy t={t}", entropy(final), entropy(initial))
        check(f"total Hamiltonian energy t={t}", np.trace(h@final), np.trace(h@initial))
        rows.append({"t":t,"delta_S_A":ds,"beta_delta_E_B":beta*de,
                     "mutual_information":mi,"relative_entropy":rel})
    return rows


def agency():
    joint = .5*np.array([[.9,.1],[.1,.9]])
    scrambled = np.outer(joint.sum(1),joint.sum(0))
    success = float(np.trace(joint))
    null = float(np.trace(scrambled))
    passive = float(joint[0].sum())
    mi = float(sum(joint[e,m]*np.log2(joint[e,m]/scrambled[e,m])
                   for e in range(2) for m in range(2)))
    check("informative controller success", success,.9)
    check("scrambling destroys advantage", null,.5)
    check("passive control",passive,.5)
    check("memory marginals preserved",joint.sum(0),scrambled.sum(0))
    check("environment marginals preserved",joint.sum(1),scrambled.sum(1))
    return {"success_with_memory":success,"success_scrambled":null,
            "passive_success":passive,"mutual_information_bits":mi,
            "supplied_inputs":["sensor error .1","policy action=memory",
                               "uniform environment","success function action=environment"],
            "not_demonstrated":["origin of controller","life","conscious experience",
                                "zero energetic cost of preparing/scrambling memory"]}


def main():
    out={"scope":"same-author finite mathematical illustrations; analytic answers known",
         "whole_parts":whole_parts(), "records":records(),
         "entropy_identity":entropy_identity(), "agency":agency()}
    out.update(checks=checks,check_count=len(checks),
               all_checks_pass=all(c["pass"] for c in checks),
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               runtime={"numpy":np.__version__,"scipy":scipy.__version__})
    (ROOT/"results").mkdir(exist_ok=True)
    (ROOT/"results"/"witnesses.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["check_count","all_checks_pass","runtime"]}))
    if not out["all_checks_pass"]: raise SystemExit(1)


if __name__ == "__main__": main()
