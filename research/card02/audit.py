"""Terminal-failure documentation for one frozen law; no rescue or law search.

Exact Gaussian-rational controls are independent of the numpy law arithmetic.
Public consistency fixtures are not fresh charter holdouts or external review.
"""
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np
from scipy.integrate import solve_ivp
import law

ROOT = Path(__file__).resolve().parent
FREEZE_SHA = "b871073fc2c7e02097dd82657480c6a667a012c4"


@dataclass(frozen=True)
class G:
    r: F = F(0)
    i: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "r", F(self.r))
        object.__setattr__(self, "i", F(self.i))

    def __add__(self, other):
        other = asg(other)
        return G(self.r+other.r, self.i+other.i)

    __radd__ = __add__

    def __neg__(self):
        return G(-self.r, -self.i)

    def __sub__(self, other):
        return self + (-asg(other))

    def __rsub__(self, other):
        return asg(other)-self

    def __mul__(self, other):
        other = asg(other)
        return G(self.r*other.r-self.i*other.i, self.r*other.i+self.i*other.r)

    __rmul__ = __mul__

    def conjugate(self):
        return G(self.r, -self.i)


def asg(x):
    return x if isinstance(x, G) else G(x)


def mat(rows):
    return [[asg(x) for x in row] for row in rows]


def zero(d):
    return mat([[0]*d for _ in range(d)])


def eye(d):
    return mat([[int(i == j) for j in range(d)] for i in range(d)])


def add(a, b):
    return [[x+y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(s, a):
    return [[asg(s)*x for x in row] for row in a]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), G())
             for j in range(len(b[0]))] for i in range(len(a))]


def adj(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def comm(a, b):
    return add(mul(a, b), scale(-1, mul(b, a)))


def tr(a):
    return sum((a[i][i] for i in range(len(a))), G())


def kron(a, b):
    return [[a[i][j]*b[k][l] for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]


P = [mat([[0, 1], [1, 0]]), mat([[0, G(0,-1)], [G(0,1), 0]]), mat([[1,0],[0,-1]])]
R0 = mat([[1,0],[0,0]])
R1 = mat([[0,0],[0,1]])
RP = scale(F(1,2), mat([[1,1],[1,1]]))


def qops(n):
    out = []
    for site in range(n):
        row = []
        for p in P:
            q = mat([[1]])
            for i in range(n):
                q = kron(q, p if i == site else eye(2))
            row.append(q)
        out.append(row)
    return out


def exact(rho, ops):
    means = [[tr(mul(rho, q)).r for q in row] for row in ops]
    cov, h = {}, zero(len(rho))
    for i, j in combinations(range(len(ops)), 2):
        cov[i,j] = []
        for a in range(3):
            row = []
            for b in range(3):
                c = tr(mul(rho, mul(ops[i][a],ops[j][b]))).r-means[i][a]*means[j][b]
                row.append(c)
                bracket = add(add(mul(ops[i][a],ops[j][b]), scale(-means[j][b],ops[i][a])),
                              scale(-means[i][a],ops[j][b]))
                h = add(h,scale(c,bracket))
            cov[i,j].append(row)
    deriv = scale(G(0,-1),comm(h,rho))
    return means,cov,h,deriv


def arr(a):
    return np.array([[complex(float(x.r),float(x.i)) for x in row] for row in a])


def serial(a):
    return [[[str(x.r),str(x.i)] for x in row] for row in a]


def save_json(path, value):
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")


def register(name, n, rho, controls=None):
    """Log every candidate fixture before executing its generator/energy/extractor.

    Replay de-duplicates the identical hashed instance, not a different input.
    """
    payload = {"id":name,"n":n,"kappa":"1","controls":controls or "zero",
               "rho":serial(rho) if isinstance(rho,list) else
               [[[float(x.real),float(x.imag)] for x in row] for row in rho]}
    digest = hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    path=ROOT/"ATTEMPT_LEDGER.json"
    ledger=json.loads(path.read_text())
    rows=[r for r in ledger["entries"] if r.get("instance_sha256")==digest]
    if not rows:
        if any(r["id"]==name for r in ledger["entries"]):
            raise RuntimeError("changed instance under existing ID")
        ledger["entries"].append({"id":name,"type":"public_candidate_instance_conservatively_taxed",
            "counts_for_IP12":True,"instance_sha256":digest,
            "logged_at_utc":datetime.now(timezone.utc).isoformat(),"instance":payload})
        save_json(path,ledger)
    return digest


def main():
    freeze=json.loads((ROOT/"FREEZE.json").read_text())
    for name,digest in freeze["files"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,("frozen file changed",name)
    receipt=json.loads((ROOT/"FREEZE_RECEIPT.json").read_text())
    assert receipt["remote_freeze_commit"]==FREEZE_SHA
    ops=qops(2)
    bell=zero(4)
    for i,j in product([0,3],repeat=2):bell[i][j]=G(F(1,2))
    schmidt=zero(4)
    schmidt[0][0],schmidt[3][3]=G(F(16,25)),G(F(9,25))
    schmidt[0][3]=schmidt[3][0]=G(F(12,25))
    mixture=scale(F(1,2),add(kron(R0,R0),kron(RP,R1)))
    parity=zero(8)
    for i in range(8):
        if i.bit_count()%2==0:parity[i][i]=G(F(1,4))
    classical=scale(F(1,2),add(kron(R0,R0),kron(R1,R1)))
    states=[("product-plus-plus",2,kron(RP,RP)),("product-zero-zero",2,kron(R0,R0)),
        ("product-plus-one",2,kron(RP,R1)),("correlated-classical-mixture",2,mixture),
        ("rational-Schmidt-16-25",2,schmidt),("maximally-mixed",2,scale(F(1,4),eye(4))),
        ("three-body-parity",3,parity),("two-correlated-blocks",4,kron(bell,classical)),
        ("mixture-times-product-third",3,kron(mixture,R0))]
    rows=[]; max_rhs=0.0; noether_cells=0; exact_controls=0; local_drive_numeric_controls=0
    data={}
    for name,n,rho in states:
        digest=register(name,n,rho)
        oo=qops(n); means,cov,h,deriv=exact(rho,oo)
        assert h==adj(h) and deriv==adj(deriv)
        assert tr(deriv)==G()
        local=[tr(mul(deriv,q)) for row in oo for q in row]
        assert all(x==G() for x in local)
        noether_cells+=len(local); exact_controls+=3+len(local)
        numeric=law.generator(arr(rho),law.operators(n))
        err=float(np.max(np.abs(numeric-arr(deriv)))); max_rhs=max(max_rhs,err)
        assert err<2e-12
        cov_zero=all(x==0 for c in cov.values() for row in c for x in row)
        rows.append({"id":name,"n":n,"instance_sha256":digest,"pair_covariance_zero":cov_zero,
            "exact_derivative_zero":deriv==zero(len(rho)),"trace_derivative":"0",
            "local_marginal_derivatives":["0"]*len(local),"numeric_exact_discrepancy":err,
            "attempted_chain":law.attempted_partition(arr(rho),law.operators(n))})
        data[name]=(rho,oo,means,cov,h,deriv)
    # Exact failure witnesses, not repairs or additional candidate laws.
    for name in ["product-plus-plus","product-zero-zero","product-plus-one","three-body-parity","maximally-mixed"]:
        assert data[name][5]==zero(len(data[name][0])); exact_controls+=1
    assert tr(mul(parity,kron(kron(P[2],P[2]),P[2])))==G(1); exact_controls+=1
    _,_,_,_,hm,dm=data["correlated-classical-mixture"]
    yz=kron(P[1],P[2])
    assert tr(mul(dm,yz))==G(1); exact_controls+=1
    assert data["product-zero-zero"][5]==zero(4) and data["product-plus-one"][5]==zero(4)
    # Rational Schmidt derivative agrees with the frozen phase equation.
    z=F(7,25); omega=-4*z*(2-z*z)
    expected=zero(4)
    expected[0][3]=G(0,-omega*F(12,25)); expected[3][0]=expected[0][3].conjugate()
    assert data["rational-Schmidt-16-25"][5]==expected; exact_controls+=1
    # Exact first variation of F at the mixed-state witness.
    delta=scale(F(1,32),kron(P[0],P[2]))
    means,cov=data["correlated-classical-mixture"][2:4]
    dmean=[[tr(mul(delta,q)).r for q in row] for row in ops]
    directional=F(0)
    for a,b in product(range(3),repeat=2):
        dc=tr(mul(delta,mul(ops[0][a],ops[1][b]))).r-dmean[0][a]*means[1][b]-means[0][a]*dmean[1][b]
        directional+=cov[0,1][a][b]*dc
    assert G(directional)==tr(mul(hm,delta)); exact_controls+=1
    assert directional==F(-1,16); exact_controls+=1
    # Exact product-block closure on the three-body witness.
    assert data["mixture-times-product-third"][5]==kron(dm,R0); exact_controls+=1
    # Basis covariance: Hadamard rotation, evaluated with rational matrices.
    v=kron(mat([[1,1],[1,-1]]),eye(2))
    rotated=scale(F(1,2),mul(mul(v,schmidt),adj(v)))
    digest=register("local-Hadamard-Schmidt",2,rotated)
    _,_,hr,dr=exact(rotated,ops)
    assert dr==scale(F(1,2),mul(mul(v,expected),adj(v))); exact_controls+=1
    rows.append({"id":"local-Hadamard-Schmidt","n":2,"instance_sha256":digest,"exact_local_basis_covariance":True})
    # A local drive can rotate a product state but cannot provide an interaction.
    controls=[[F(1,3),0,0],[0,0,0]]
    rho=kron(R0,R0); digest=register("product-local-drive",2,rho,controls=[[str(x) for x in row] for row in controls])
    driven=law.generator(arr(rho),law.operators(2),local_controls=np.array(controls,float))
    expected_drive=scale(G(0,-1),comm(scale(F(1,3),ops[0][0]),rho))
    assert np.max(np.abs(driven-arr(expected_drive)))<1e-14; local_drive_numeric_controls+=1
    rows.append({"id":"product-local-drive","n":2,"instance_sha256":digest,"intrinsic_hamiltonian_zero":True,
                 "local_drive_matches_product_rotation":True})
    # Exact conventional record-writing comparator at g*t=pi/4.
    initial=kron(RP,RP)
    w=zero(4)
    for i,sign in enumerate([1,-1,-1,1]):w[i][i]=G(1,-sign)
    written=scale(F(1,2),mul(mul(w,initial),adj(w)))
    system=[[written[i*2][j*2]+written[i*2+1][j*2+1] for j in range(2)] for i in range(2)]
    e0=scale(2,[[written[i][j] for j in range(2)] for i in range(2)])
    e1=scale(2,[[written[i+2][j+2] for j in range(2)] for i in range(2)])
    assert system==scale(F(1,2),eye(2)) and tr(mul(e0,e1))==G(); exact_controls+=2
    # Fixed exchange comparator: complete transfer at g*t=pi/2.
    exchange_h=scale(F(1,2),add(kron(P[0],P[0]),kron(P[1],P[1])))
    exchange_projector=mul(exchange_h,exchange_h)
    exchange_u=add(add(eye(4),scale(-1,exchange_projector)),scale(G(0,-1),exchange_h))
    assert mul(adj(exchange_u),exchange_u)==eye(4); exact_controls+=1
    exchange_final=mul(mul(exchange_u,kron(R1,R0)),adj(exchange_u))
    assert exchange_final==kron(R0,R1); exact_controls+=1
    assert tr(mul(exchange_final,kron(R1,eye(2))))==G()
    assert tr(mul(exchange_final,kron(eye(2),R1)))==G(1); exact_controls+=2
    # Numerical trajectory checks validate the frozen nonlinear phase forecast.
    trajectories=[]; max_trajectory=0.0
    times=np.linspace(0,1,9)
    for p in [F(3,4),F(5,8)]:
        psi=np.array([math.sqrt(float(p)),0,0,math.sqrt(float(1-p))],complex)
        rho=np.outer(psi,psi.conj()); nops=law.operators(2)
        digest=register("phase-p-"+str(p).replace('/','-'),2,rho)
        sol=solve_ivp(lambda t,y:law.generator(y.reshape(4,4),nops).reshape(-1),
            [0,1],rho.reshape(-1),t_eval=times,method="DOP853",rtol=1e-11,atol=1e-12)
        assert sol.success
        rate=-4*(2*p-1)*(2-(2*p-1)**2)
        outputs=[]
        for t,y in zip(sol.t,sol.y.T):
            state=y.reshape(4,4)
            ideal=psi.copy(); ideal[3]*=np.exp(1j*float(rate)*t)
            ideal=np.outer(ideal,ideal.conj())
            err=float(np.max(np.abs(state-ideal))); max_trajectory=max(max_trajectory,err)
            assert err<2e-9
            outputs.append({"time":float(t),"matrix_discrepancy":err,"energy":law.energy(state,nops),
                "purity":float(np.trace(state@state).real),"min_eigenvalue":float(np.min(np.linalg.eigvalsh(state))),
                "xx":float(np.trace(state@nops[0][0]@nops[1][0]).real),
                "yx":float(np.trace(state@nops[0][1]@nops[1][0]).real)})
        trajectories.append({"p":str(p),"omega":str(rate),"instance_sha256":digest,"outputs":outputs})
    assert F(31,56)==F(trajectories[1]["omega"])/F(trajectories[0]["omega"])
    # Price only the explicit arithmetic; do not fabricate full Price.
    sys.path.insert(0,str(ROOT.parents[1]))
    from development.final_precard.price_oracle import raw_price,table_price
    ast=json.loads((ROOT/"L0_CORE.json").read_text())
    prices={str(p):raw_price(ast["expression"],p,definitions=ast["definitions"]) for p in [6,10,16]}
    seed_literals=[["rat",x.r.numerator,x.r.denominator] for q in P for row in q for x in row]
    seed_literals += [["rat",x.i.numerator,x.i.denominator] for q in P for row in q for x in row]
    seed_bits=table_price(seed_literals)
    ledger=json.loads((ROOT/"ATTEMPT_LEDGER.json").read_text())
    current=sum(bool(e["counts_for_IP12"]) for e in ledger["entries"])
    local=current+ledger["prior_card1_disclosed_entries"]
    price={"pointwise_arithmetic":prices,"Pauli_seed_literal_table_bits":seed_bits,
        "card2_taxed_entries":current,"prior_card1_disclosed_entries":7,
        "local_taxed_entries":local,"selection_tax_local_lower_bound_bits":math.log2(1+local),
        "core_plus_seed_plus_tax_lower_bound_bits":prices["10"]["bits"]+seed_bits+math.log2(1+local),
        "full_price_certified":False,"oracle_external_clearance":False,"awarded_b_J":0,
        "normal_form_and_import_unfolding_complete":False,
        "unpriced_obligations":["quantum/time/domain axioms","input/extractor/embedding statements",
            "full constant range/pass windows","selection maxima","IP13 exclusions","full vocabulary/checklists"]}
    save_json(ROOT/"price_report.json",price)
    result={"status":"FROZEN_CARD_2_AUTHOR_KILL_PENDING_REVIEW","formal_first_terminal":"S0_CARD_INCOMPLETE",
        "scientific_failure":"Local-unitary invariant correlation energy conserves all local marginals; product-block organization cannot arise across initially uncorrelated cuts",
        "freeze_commit":FREEZE_SHA,"candidate_law_changed":False,"candidate_fixtures":rows,
        "exact_controls":exact_controls,"numeric_local_drive_controls":local_drive_numeric_controls,
        "exact_local_marginal_cells":noether_cells,
        "max_reference_exact_rhs_discrepancy":max_rhs,"max_trajectory_discrepancy":max_trajectory,
        "phase_trajectories":trajectories,"predicted_phase_rate_ratio":"31/56","fixed_standard_phase_ratio":"1",
        "affinity_witness":{"rho":"half |00><00| + half |+1><+1|",
            "mixture_YA_ZB_derivative":"1","average_product_branch_derivative":"0",
            "linear_CPTP_extension_exists":False,"superluminal_signal_proved":False},
        "conventional_record_comparator":{"hamiltonian":"g Z_A Z_B","g":"1","duration":"pi/4",
            "preparation":"|+,+>","candidate_visibility":"1","candidate_record_D":"0",
            "standard_visibility":"0","standard_record_D":"1","norm_ceiling":"27 kappa"},
        "conventional_exchange_comparator":{"hamiltonian":"g(XX+YY)/2","g":"1","duration":"pi/2",
            "preparation":"|1,0>","candidate_B_excitation_probability":"0","standard_B_excitation_probability":"1"},
        "charter":{"generated_partition":"FAIL","common_carrier":"UNRESOLVED_FAIL",
            "generated_readout_interface":"UNRESOLVED_FAIL","epsilon":"NOT_COMPUTED_CHAIN_FAILED",
            "Q5":"FAIL_NO_GENERATED_JOINT_COORDINATES","LT1":"NOT_CERTIFIED","fresh_holdouts":"NOT_DRAWN_S0_FAILED",
            "full_batteries":"NOT_RUN_AFTER_S0_FAILURE","independent_recomputer":False,"external_kit_review":False},
        "price":price,"slots_consumed":2,"slots_remaining":1,"card3_authorized":False,
        "new_validated_physics":False,"experimental_measurements":False,
        "mathematical_controls_grade":"Author exact arithmetic and analytic proofs; trajectories are numerical evidence"}
    save_json(ROOT/"results.json",result)
    table=["# Generated Card-2 controls", "", "These are failure-documentation controls, not a passing charter score.", "",
        "| Quantity | Result |","|---|---|",
        f"| Logged Card-2 entries, including symbolic | {current} |",
        f"| Exact rational control assertions | {exact_controls} |",
        f"| Numerical local-drive comparison | {local_drive_numeric_controls} |",
        f"| Exact local-marginal derivative cells | {noether_cells} zero |",
        f"| Maximum reference/exact RHS discrepancy | {max_rhs:.6g} |",
        f"| Maximum finite-trajectory discrepancy | {max_trajectory:.6g} |",
        "| Frozen phase ratio / fixed standard comparator | 31/56 / 1 |",
        "| Standard record-writing outcome (V,D) | (0,1) |",
        "| Candidate product-input outcome (V,D) | (1,0) |",
        f"| Arithmetic + seed + local tax price lower bound (p=10) | {price['core_plus_seed_plus_tax_lower_bound_bits']:.6f} bits |",
        "| Complete price / independent review | Unestablished / pending |",
        "| Charter first terminal | S0 CARD-INCOMPLETE |",
        "| Scientific author verdict | KILL |"]
    (ROOT/"CONTROL_TABLES.md").write_text("\n".join(table)+"\n")
    print(result["status"])
    print(f"{exact_controls} exact controls; {local_drive_numeric_controls} numerical drive control; {noether_cells} local-marginal cells; {current} Card-2 entries")
    print("Numerical checks document the failed frozen law; no passing card or experiment.")


if __name__ == "__main__":
    main()
