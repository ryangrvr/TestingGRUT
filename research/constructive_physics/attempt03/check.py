#!/usr/bin/env python3
"""CP-03: exact algebra and independently integrated leading scalar bubble.

No laboratory data. No all-orders quantum or Lorentzian-reconstruction claim.
The bubble is generated from its momentum polynomial, not a stored K formula.
"""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
GROUPS = {}


def exact_array(a):
    return np.array([[F(x) for x in row] for row in a], dtype=object)


def eye(n):
    return exact_array(np.eye(n, dtype=int))


def inverse(a):
    n = len(a)
    r = np.concatenate((a.copy(), eye(n)), axis=1)
    for i in range(n):
        pivot = next(j for j in range(i, n) if r[j, i])
        r[[i, pivot]] = r[[pivot, i]]
        r[i] = r[i] / r[i, i]
        for j in range(n):
            if j != i:
                r[j] -= r[j, i] * r[i]
    return r[:, n:]


def rank(a):
    r = a.copy()
    nrow, ncol = r.shape
    i = 0
    for col in range(ncol):
        pivots = [j for j in range(i, nrow) if r[j, col]]
        if not pivots:
            continue
        j = pivots[0]
        r[[i, j]] = r[[j, i]]
        r[i] /= r[i, col]
        for j in range(nrow):
            if j != i:
                r[j] -= r[j, col] * r[i]
        i += 1
        if i == nrow:
            break
    return i


def tr(a):
    return sum(a[i, i] for i in range(len(a)))


def determinant(a):
    total = F(0)
    for perm in permutations(range(len(a))):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(len(a)) for j in range(i+1, len(a)))
        term = F(sign)
        for i, j in enumerate(perm):
            term *= a[i, j]
        total += term
    return total


def serial(a):
    if isinstance(a, np.ndarray):
        return serial(a.tolist())
    if isinstance(a, F):
        return str(a)
    if isinstance(a, (list, tuple)):
        return [serial(x) for x in a]
    if isinstance(a, dict):
        return {str(k): serial(v) for k, v in a.items()}
    if isinstance(a, np.generic):
        return a.item()
    return a


# Polynomial variables (k0,k1,k2,k3,t,x); t^2=Q_E^2.
ZERO = (0,)*6


def const(c):
    return {ZERO: F(c)} if c else {}


def variable(i):
    ex = list(ZERO); ex[i] = 1
    return {tuple(ex): F(1)}


def add(*ps):
    out = {}
    for p in ps:
        for ex, c in p.items():
            out[ex] = out.get(ex, F(0)) + c
    return {ex:c for ex,c in out.items() if c}


def mul(*ps):
    out = const(1)
    for p in ps:
        new = {}
        for a, ca in out.items():
            for b, cb in p.items():
                ex = tuple(x+y for x,y in zip(a,b))
                if ex[4] <= 2:  # Only the frozen small-momentum jet.
                    new[ex] = new.get(ex, F(0)) + ca*cb
        out = {ex:c for ex,c in new.items() if c}
    return out


def scale(c, p):
    return mul(const(c), p)


def integrate_bubble_polynomial(p, m2=F(1)):
    """Use isotropic DR integrals then integrate x exactly on [0,1].

    Common factor 1/(16*pi^2*epsilon) is stripped. J0=2,
    <k_i^2>=-Delta, <k_i^4>=3 Delta^2/4,
    <k_i^2 k_j^2>=Delta^2/4, Delta=m^2+x(1-x)t^2.
    Odd monomials vanish. These follow from radial Gamma integrals.
    """
    x, t = variable(5), variable(4)
    z = add(x, scale(-1, mul(x,x)))
    delta = add(const(m2), mul(z,t,t))
    out = [F(0), F(0), F(0)]
    for ex, coeff in p.items():
        powers = ex[:4]
        if any(j % 2 for j in powers):
            continue
        degree = sum(powers)
        if degree == 0:
            moment = const(2)
        elif degree == 2:
            moment = scale(-1, delta)
        elif degree == 4:
            moment = scale(F(3,4) if 4 in powers else F(1,4), mul(delta,delta))
        else:
            raise AssertionError(degree)
        base = list(ZERO); base[4],base[5] = ex[4],ex[5]
        for power, c in mul({tuple(base):coeff}, moment).items():
            out[power[4]] += c / (power[5]+1)
    assert out[1] == 0
    return out[0], out[2]


def derive_bubble():
    ks = [variable(i) for i in range(4)]
    t, x = variable(4), variable(5)
    b = add(x, const(F(-1,2)))
    z = add(x, scale(-1,mul(x,x)))
    k2 = add(*(mul(k,k) for k in ks))
    e = [[{} for _ in range(4)] for _ in range(4)]
    for mu,nu in product(range(4),repeat=2):
        diag = int(mu == nu)
        # Feynman-shifted Euclidean stress insertion, derived before execution.
        e[mu][nu] = add(
            scale(2,mul(ks[mu],ks[nu])),
            scale(-diag,add(k2,const(1))),
            scale(2,mul(b,t,add(scale(int(nu==3),ks[mu]),
                scale(int(mu==3),ks[nu]),scale(-diag,ks[3])))),
            scale(diag-2*int(mu==3 and nu==3),mul(z,t,t)))
    b0,b2 = np.empty((4,4,4,4),dtype=object),np.empty((4,4,4,4),dtype=object)
    for mu,nu,rho,sigma in product(range(4),repeat=4):
        b0[mu,nu,rho,sigma],b2[mu,nu,rho,sigma] = integrate_bubble_polynomial(mul(e[mu][nu],e[rho][sigma]))
    pi = np.empty((4,4,4,4),dtype=object)
    for mu,nu,rho,sigma in product(range(4),repeat=4):
        pi[mu,nu,rho,sigma] = F(int(mu==rho and nu==sigma)+int(mu==sigma and nu==rho)-int(mu==nu and rho==sigma),2)
    assert np.array_equal(b0,2*pi)
    spatial = []
    for i in range(3):
        for j in range(i,3):
            a=eye(4)*0;a[i,j]=a[j,i]=F(1);spatial.append(a)
    c0=(pi.reshape(16,16)@b0.reshape(16,16)).reshape(4,4,4,4)
    c2=(pi.reshape(16,16)@b2.reshape(16,16)).reshape(4,4,4,4)
    cells=0
    for source in spatial:
        r0=(c0.reshape(16,16)@source.reshape(16)).reshape(4,4)
        r2=(c2.reshape(16,16)@source.reshape(16)).reshape(4,4)
        assert np.array_equal(r0[:3,:3],2*source[:3,:3])
        assert np.array_equal(r2[:3,:3],F(2,3)*source[:3,:3])
        cells += 18
    GROUPS['independent_scalar_bubble']={
        'method':'Momentum polynomial -> Feynman shift -> radial/angular DR integrals -> exact Feynman-parameter integral -> Pi contraction -> conserved-source quotient',
        'exact_transverse_cells':cells,'B0':b0,'B_Q2':b2,
        'common_factor':'1/(16*pi^2*epsilon)',
        'contracted_coefficients':['2 m^4','2 m^2 Q_E^2/3'],
        'K_Lorentz':'-N*m^4/(32*pi^2*V0*epsilon)*(1-q^2/(3*m^2))',
        'finite_parts_or_full_q_dependence_computed':False}


def classical():
    fixtures=[]
    for n in (4,6,8):
        g=eye(4)
        if n>4:
            extra=exact_array([[(i+1)*(j+1)%7-3 for j in range(n-4)] for i in range(4)])
            g=np.concatenate((g,extra),axis=1)
        q=g@g.T;qi=inverse(q);v=F(7)
        metric=q/v;mi=inverse(metric)
        assert np.array_equal(q-metric*(tr(mi@q)/2-v),eye(4)*0)
        assert determinant(q)>0
        assert all(determinant(q[:k,:k])>0 for k in range(1,5))
        l=math.sqrt(float(determinant(q)))/float(v)
        gf=np.array(g,dtype=float);qif=np.array(qi,dtype=float)
        p=l*(qif[0]@gf)
        assert np.max(abs(p@gf.T-np.array([l,0,0,0])))<1e-11
        normal=eye(n)-g.T@qi@g
        assert rank(normal)==n-4 and np.array_equal(normal@g.T,np.zeros((n,4),dtype=object))
        h=l*qif[0,0]*np.array(normal,dtype=float)
        # Independent numerical derivative of the canonical momentum.
        def momentum(v0):
            a=gf.copy();a[0]=v0;qa=a@a.T
            la=math.sqrt(float(np.linalg.det(qa)))/7
            return la*(np.linalg.inv(qa)[0]@a)
        step=1e-5
        actual=np.column_stack([(momentum(gf[0]+step*np.eye(n)[i])-momentum(gf[0]-step*np.eye(n)[i]))/(2*step) for i in range(n)])
        err=float(np.max(abs(actual-h)))
        assert err<1e-7
        fixtures.append({'N':n,'gradient_matrix':g,'Gram_matrix':q,'all_leading_minors_positive':True,'Hessian_rank_exact':n-4,'Hessian_derivative_error':err})
    GROUPS['classical_metric_and_signature']={'fixtures':fixtures,
        'theorem':'v^mu Q_mu_nu v^nu=sum_a(v^mu partial_mu phi^a)^2 >=0; at rank 4 Q is positive definite. Multiplying by a nonzero scalar cannot give Lorentzian inertia.',
        'real_Lorentzian_classical_branch':False,
        'quantum_Lorentzian_continuation_refuted':False}


def canonical():
    # Fixed local spatial embedding, arbitrary off-constraint momentum.
    g=exact_array([[1,0,0,1,2,0],[0,1,0,0,1,2],[0,0,1,2,0,1]])
    q=g@g.T;qi=inverse(q);d=determinant(q);u=F(1,7)
    p=np.array([F(1),F(2),F(3),F(4),F(5),F(6)],dtype=object)
    hi=g@p
    dh_dg=-u*u*d*(qi@g)
    bracket_coefficient=dh_dg@p
    expected=-u*u*d*(qi@hi)
    assert np.array_equal(bracket_coefficient,expected)
    GROUPS['constraint_algebra']={'spatial_gradient':g,'q_inverse':qi,
        'H_i':hi,'bracket_coefficient':bracket_coefficient,
        'H0':'(p^2-U(phi)^2 det(q_ij))/2',
        'exact_smeared_bracket':'{H0[n],H0[m]} = D[-U^2 det(q) q^ij (n partial_j m-m partial_j n)]',
        'normalized_weak_bracket':'{H[N],H[M]} approximately -D[q^ij (N partial_j M-M partial_j N)]',
        'H':'H0/(U sqrt(det q))',
        'signature':'Euclidean normal-deformation sign for positive target metric',
        'ADM_metric_and_quantum_constraint_algebra_derived':False,
        'constraints_inherited_from_postulated_reparametrization_invariance':True}


def quartic():
    eta=exact_array(np.diag([1,-1,-1,-1]));v=F(7)
    bs=[eye(4),exact_array([[1,2,0,1],[2,3,1,0],[0,1,2,1],[1,0,1,4]]),
        exact_array([[2,0,1,0],[0,1,0,2],[1,0,3,0],[0,2,0,1]])]
    fixtures=[]
    for b in bs:
        a=eta@b/v
        # Direct permutation determinant coefficients through quadratic order.
        d=[F(0)]*3
        for perm in permutations(range(4)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
            coeff=[F(sign),F(0),F(0)]
            for i,j in enumerate(perm):
                c0,c1=F(i==j),a[i,j]
                coeff=[coeff[0]*c0,coeff[1]*c0+coeff[0]*c1,coeff[2]*c0+coeff[1]*c1]
            d=[x+y for x,y in zip(d,coeff)]
        assert d[0]==1
        for potential in (F(0),F(2,5)):
            interaction=v*(d[2]/2-d[1]*d[1]/8-potential*d[1]/(2*v)+potential**2/v**2)
            stress=b-eta*(tr(eta@b)/2-potential)
            stress_route=-(tr((eta@stress)@(eta@stress))-tr(eta@stress)**2/2)/(4*v)
            assert interaction==stress_route
            fixtures.append({'B':b,'potential_delta':potential,'quartic_coefficient':interaction})
    GROUPS['determinant_to_total_stress_interaction']={'exact_fixtures':fixtures,'L_int':'-(T_mu_nu T^mu_nu-T^2/2)/(4 V0)'}


def helicities():
    eta=exact_array(np.diag([1,-1,-1,-1]))
    basis=[]
    for j in range(6):
        a,b,c,d,e,f=[F(i==j) for i in range(6)]
        basis.append(exact_array([[a,b,c,a],[b,d,f,b],[c,f,e,c],[a,b,c,a]]))
    q=np.array([F(1),F(0),F(0),F(1)],dtype=object)
    for source in basis:
        assert all(x==0 for x in (q@eta)@source)
    def residue(a,b):
        return tr(a@eta@b@eta)-tr(eta@a)*tr(eta@b)/2
    gram=np.array([[residue(a,b) for b in basis] for a in basis],dtype=object)
    assert rank(gram)==2
    eig=np.linalg.eigvalsh(np.array(gram,dtype=float))
    assert min(eig)>-1e-13
    for vals in ((1,2,3,4,5,6),(0,1,0,2,2,0),(1,0,0,0,0,1)):
        source=sum((F(c)*b for c,b in zip(vals,basis)),start=eye(4)*0)
        assert residue(source,source)==F(vals[3]-vals[4])**2/2+2*F(vals[5])**2
    timelike=exact_array(np.diag([0,1,1,1]))
    assert residue(timelike,timelike)==F(-3,2)
    GROUPS['conserved_source_helicity_residue']={'null_q':[1,0,0,1],
        'residue_Gram':gram,'rank':2,'eigenvalues':eig,
        'exact_positive_form':'(T11-T22)^2/2+2 T12^2',
        'timelike_trace_control':residue(timelike,timelike),
        'interpretation':'Two positive radiative helicities in the conditional massless exchange. Timelike negative trace is not proof of a physical ghost at that point; a massive extension needs a separate constraint/positivity analysis.'}


def tuning_and_universality():
    n,m,eps=200,1.0,-0.01
    v=-n*m**4/(32*math.pi**2*eps)
    mp2=-n*m*m/(48*math.pi**2*eps)
    lam=n*m**4/(32*math.pi**2*v*eps)
    assert abs(lam+1)<1e-14 and abs(2*v/(3*m*m)-mp2)<1e-13
    rows=[]
    for qe2 in (0.01,0.001,0.0001):
        kernel=-lam*(1+qe2/(3*m*m))
        a=-1/(2*v)/(1-kernel)
        pole=1/(mp2*qe2)
        err=abs(a/pole-1)
        assert err<1e-10
        rows.append({'Q_E2_over_m2':qe2,'resummed_scalar_prefactor':a,'pole_prefactor':pole,'relative_error':err})
    sensitivities=[]
    for delta in (-0.01,0.01):
        lam_new=-1/(1+delta)
        # Diagnostic root in the low-q denominator only; background not certified.
        q2=3*m*m*(1+lam_new)/lam_new
        sensitivities.append({'V0_fractional_shift':delta,'lambda':lam_new,'diagnostic_q2_pole':q2,'complete_flat_vacuum_model':False})
    s4=F(120*1**4+80*2**4);s2=F(120*1**2+80*2**2)
    b=s2/(3*s4)
    vp=float(s4)/(-32*math.pi**2*eps)
    mp2p=float(s2)/(-48*math.pi**2*eps)
    assert abs(1/(2*vp*float(b))-1/mp2p)<1e-14
    # One total stress square forces equal cross/mixed couplings.
    species_prefactors=[-1/(2*vp*float(b)) for _ in range(3)]
    assert max(species_prefactors)==min(species_prefactors)
    d=4-eps
    full_v=-n*(d/2-1)*math.gamma(-d/2)/(2*(4*math.pi)**(d/2))
    full_mp=n*math.gamma(1-d/2)/(6*(4*math.pi)**(d/2))
    assert full_v>0 and full_mp>0
    GROUPS['critical_pole_and_regulator_sensitivity']={'software_N':n,'software_m':m,'epsilon':eps,
        'leading_V0':v,'leading_M_P2':mp2,'lambda':lam,'rows':rows,'detuning_diagnostics':sensitivities,
        'full_D_V0_over_leading':full_v/v,'full_D_MP_power_over_leading_MP2':full_mp/mp2,
        'physical_regulator_or_tuning_generated':False}
    GROUPS['two_species_universal_scalar_coupling']={'species':[[120,1],[80,2]],'S4':s4,'S2':s2,
        'small_momentum_coefficient':b,'same_pole_prefactor_aa_ab_bb':species_prefactors,
        'generated_observed_matter_inventory':False,'coupling_to_entire_gravitational_energy_independently_proved':False}
    # Independent heat-kernel sign/coefficient route: Gamma(-2+eps/2)=1/eps,
    # Gamma(-1+eps/2)=-2/eps; a1=R/6 for a minimally coupled scalar.
    gamma1_pole=F(-2)
    coeff_R=-F(n,2)*F(1,6)*gamma1_pole/F(16)  # pi^2 epsilon stripped.
    mp2_pole=-2*coeff_R
    assert mp2_pole==F(-n,48)
    GROUPS['heat_kernel_normalization']={'R_coefficient_stripped_pi2_epsilon':coeff_R,
        'M_P2_coefficient_stripped_pi2_epsilon':mp2_pole,
        'quantum_equivalence_of_auxiliary_and_constrained_measure':False}


def newton_limit():
    eta=exact_array(np.diag([1,-1,-1,-1]))
    sources=[]
    for mass in (F(1),F(2)):
        p=np.array([mass,F(0),F(0),F(0)],dtype=object)
        pp=-p
        vertex=-(np.outer(p,pp)+np.outer(pp,p))+eta*((p@eta@pp)+mass*mass)
        expected=eye(4)*0;expected[0,0]=2*mass*mass
        assert np.array_equal(vertex,expected)
        sources.append(vertex)
    numerator=tr(sources[0]@eta@sources[1]@eta)-tr(eta@sources[0])*tr(eta@sources[1])/2
    assert numerator==2*F(1)**2*F(2)**2
    # Dividing by the four external nonrelativistic normalization factors and
    # Fourier transforming 1/k^2 gives -M1*M2/(8*pi*M_P^2*r).
    born_coefficient=-numerator/(4*F(1)*F(2))
    assert born_coefficient==F(-1)*F(1)*F(2)/2
    GROUPS['conditional_Newton_limit']={'source_vertices':sources,
        'static_exchange_numerator':numerator,'Born_potential_coefficient_without_MP2_k2':born_coefficient,
        'Fourier_identity':'integral d^3k exp(i k.r)/[(2*pi)^3 k^2]=1/(4*pi*r), r>0',
        'potential':'-M1*M2/(8*pi*M_P^2*r)','G_N':'1/(8*pi*M_P^2)',
        'Lorentzian_quantum_theory_and_asymptotic_states_assumed_for_this_limit':True,
        'new_GRUT_prediction':False}


def invariants():
    freeze=json.loads((ROOT/'FREEZE.json').read_text())
    for name,expected in freeze['sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected
    baseline=json.loads((ROOT/'BASELINE_INTEGRITY.json').read_text())
    repo=ROOT.parents[2]
    for name,expected in {**baseline['sha256'],**baseline['extra_frozen_files']}.items():
        assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==expected,name
    GROUPS['frozen_record_preservation']={'unchanged_manifested_files':baseline['manifested_count'],
        'additional_frozen_program_files':len(baseline['extra_frozen_files']),
        'CP03_freeze_unchanged':True}


def main():
    invariants();classical();canonical();quartic();derive_bubble();helicities();tuning_and_universality();newton_limit()
    result={'status':'CP03_CONDITIONAL_COMPOSITE_POLE_CHECKS_PASS',
        'scientific_verdict':'NOT_YET_ESTABLISHED_AS_PHYSICAL_GRAVITY_OR_NEW_GRUT_LAW',
        'real_classical_Lorentzian_branch':'EXCLUDED_FOR_FROZEN_POSITIVE_GRAM_ACTION',
        'quantum_Lorentzian_reconstruction':'UNPROVED_NOT_REFUTED_BY_CLASSICAL_GRAM_TEST',
        'independent_review':False,'experiment':False,'Stage3_Card3_opened':False,
        'software':{'python':sys.version.split()[0],'numpy':np.__version__},'groups':GROUPS}
    (ROOT/'results.json').write_text(json.dumps(serial(result),indent=2)+'\n')
    print(result['status']);print('groups='+str(len(GROUPS)))
    print(result['scientific_verdict'])
    print(result['real_classical_Lorentzian_branch'])


if __name__=='__main__':
    main()
