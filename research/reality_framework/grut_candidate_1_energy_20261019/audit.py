#!/usr/bin/env python3
"""Root audit: fresh stricter continuum replay and separate FD determinants.

Executed calculations may write their own machine-readable records. No action,
matching or reference-profile parameter is changed by this audit.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_banded
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import splu

import graphs
import spectral as sp

HERE = Path(__file__).resolve().parent


def save(name, data):
    (HERE / name).write_text(json.dumps(data, indent=2) + '\n')


def product(a, b):
    q = [a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]]
    return min(q), max(q)


def square(a):
    return (0. if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2), max(a[0]**2, a[1]**2))


def interval(a, factor, b, factor2=1., offset=0.):
    ends = [factor*x + factor2*y + offset for x in a for y in b]
    return min(ends), max(ends)


def enclose(R):
    """Bernstein convex-hull range on every cell, including a split endpoint."""
    edges = np.r_[sp.SPL.x[(sp.SPL.x >= 0) & (sp.SPL.x < R)], R]
    vmax = i2 = ir3 = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        h = hi - lo
        d, c, b, a = [sp.SPL(lo, nu=k)[[0, 2]] / np.math.factorial(k) for k in range(4)] if False else [sp.SPL(lo, nu=0)[[0, 2]], sp.SPL(lo, nu=1)[[0, 2]], sp.SPL(lo, nu=2)[[0, 2]]/2, sp.SPL(lo, nu=3)[[0, 2]]/6]
        bern = np.array([d, d+c*h/3, d+2*c*h/3+b*h*h/3, d+c*h+b*h*h+a*h**3])
        F, S = [(float(bern[:,j].min())-1e-10, float(bern[:,j].max())+1e-10) for j in range(2)]
        f2, s2 = square(F), square(S)
        A = interval(s2, 3, f2, offset=-3)
        B = interval(s2, 1, f2, .2, -1)
        C = (f2[0]*.1, f2[1]*.1)
        D = tuple(2*x for x in product(F,S))
        E = (.01*(s2[0]-1), .01*(s2[1]-1))
        diag = [A, interval(B,1,C), interval(B,1,C,-1), E]
        mid = np.diag([(v[0]+v[1])/2 for v in diag])
        mid[0,1] = mid[1,0] = (D[0]+D[1])/2
        radii = np.diag([(v[1]-v[0])/2 for v in diag])
        radii[0,1] = radii[1,0] = (D[1]-D[0])/2
        vmax = max(vmax, np.max(abs(np.linalg.eigvalsh(mid))) + np.max(np.sum(radii, axis=1)))
        tr2 = sum(square(v)[1] for v in diag) + 2*square(D)[1]
        i2 += 4*np.pi*tr2*(hi**3-lo**3)/3
        ir3 += tr2*(hi**4-lo**4)/4
    return {'R':R, 'Vmax_upper_spline':float(vmax), 'I2_upper_spline':float(i2), 'Irad3_upper_spline':float(ir3), 'scope':'Padded floating-point Bernstein enclosures of supplied cubic spline, not exact BVP interval certification'}


def tail_bounds(R, V=2., I2=1800., Ir3=600., L=40, Z=10.):
    c = 1-sp.OMEGA**2
    uv = V*I2/(48*np.pi**2)*quad(lambda z: 1/((z*z+c-V)*np.sqrt(z*z+c)), Z, np.inf, epsabs=1e-12)[0]
    ell = np.arange(L+1, 100001, dtype=float)
    ratio = V*R*R/((ell-1)*(ell+2))
    if np.max(ratio) >= 1: raise ValueError('Angular expansion norm not below one')
    angular = Z*Ir3/(3*np.pi)*np.sum(ratio/((1-ratio)*(2*ell+3)*(2*ell-1)))
    # For ell>=M+1, (ell-1)(ell+2)>=ell^2 and (2ell+3)(2ell-1)>=3ell^2.
    M=100000.; rr=V*R*R/M**2
    angular_tail = Z*Ir3/(3*np.pi)*V*R*R/(9*(1-rr)*M**3)
    return {'R':R,'Vmax_used':V,'I2_used':I2,'Irad3_used':Ir3,'UV_bound':float(uv),'angular_bound':float(angular+angular_tail),'angular_series_tail_bound':float(angular_tail),'sum_bound':float(uv+angular+angular_tail)}


def independent_fd(ell, zeta, n, R=25.):
    h = R/(n+1); r = np.arange(1,n+1)*h
    k2 = np.array([1+(zeta+1j*sp.OMEGA)**2,1+(zeta-1j*sp.OMEGA)**2,2+zeta*zeta])
    V = np.array([sp.potential(x)[0] for x in r])
    G=[]; base=[]; ab=np.zeros((3,n),complex); ab[0,1:]=-1/h**2; ab[2,:-1]=-1/h**2
    for k in k2:
        diagonal = 2/h**2+ell*(ell+1)/r**2+k
        ab[1]=diagonal
        G.append(solve_banded((1,1),ab,np.eye(n),check_finite=False))
        base.append(diags([-np.ones(n-1)/h**2,diagonal,-np.ones(n-1)/h**2],[-1,0,1],format='csc'))
    blocks = [[base[i]+diags(V[:,i,i]) if i==j else diags(V[:,i,j]) for j in range(3)] for i in range(3)]
    logfull=float(np.log(abs(splu(bmat(blocks,format='csc')).U.diagonal())).sum())
    logvac=sum(float(np.log(abs(splu(x).U.diagonal())).sum()) for x in base)
    b1=sum(np.dot(np.diag(G[i]),V[:,i,i]) for i in range(3))
    b2=-.5*sum(np.einsum('mn,mn,m,n->',G[i],G[j],V[:,j,i],V[:,i,j],optimize=True) for i in range(3) for j in range(3))
    return {'ell':ell,'zeta':zeta,'n':n,'R':R,'h':h,'logdet':logfull-logvac,'Born1':float(b1.real),'Born2':float(b2.real),'remainder':float(logfull-logvac-b1.real-b2.real),'max_Born_imag':float(max(abs(b1.imag),abs(b2.imag)))}


def support_graph(R):
    r=np.linspace(0,R,int(R*200)+1); F,_,S,_=sp.SPL(r); V=graphs.potentials(F,S)
    integ=4*np.pi*graphs.simpson(r[:,None]**2*V,x=r,axis=0)
    k=np.linspace(0,40,1601); vh=graphs.transform(r,V,k); kernel,_=graphs.kernels(k,sp.OMEGA,128)
    return {'R':R,'first':float(np.sum(graphs.first_graph(integ))), 'second':float(np.sum(graphs.second_graph(k,vh,kernel)))}


def main():
    start=time.perf_counter(); raw=(HERE/'audit_protocol.json').read_bytes(); protocol=json.loads(raw)
    result={'protocol_sha256':hashlib.sha256(raw).hexdigest(),'profile_sha256':hashlib.sha256(sp.REF.read_bytes()).hexdigest(),'enclosures':[enclose(R) for R in [20.,25.]],'tails':[tail_bounds(R) for R in [20.,25.]],'independent_fd':[],'replay_points':[],'replay_partial_integrals':[]}
    for ell,zeta in protocol['independent_discretization']['points']:
        continuum=sp.run_point(ell,zeta,R=25.,tol=1e-10,r0=.001,label='FD_reference')
        for n in protocol['independent_discretization']['interior_sizes']:
            x=independent_fd(ell,zeta,n);x['continuum_remainder']=continuum['coupled_remainder'];x['discrepancy']=x['remainder']-x['continuum_remainder'];result['independent_fd'].append(x);print(json.dumps(x),flush=True)
    result['support_graphs']=[support_graph(R) for R in [20.,25.]]
    save('audit_results.json',result)
    for ell in range(41):
        subtotal=dtotal=0.
        for zeta,weight in sp.quadrature(8):
            if time.perf_counter()-start>600: raise RuntimeError('Audit budget; no complete binding result')
            rec=sp.run_point(ell,float(zeta),R=25.,tol=1e-9,r0=.001,label='audit_replay')
            if not rec['success'] or abs(rec['logJost_phase_mod_2pi'])>1e-7: raise RuntimeError('Failed replay determinant/phase')
            result['replay_points'].append(rec)
            zero=np.log(zeta*zeta/(zeta*zeta+1)) if ell in [0,1] else 0.
            subtotal+=weight*(rec['coupled_remainder']-zero);dtotal+=weight*rec['D_remainder']
        if ell in [0,1]: subtotal+=10*np.log(100/101)-2*np.arctan(10.)
        entry={'ell':ell,'energy':float((2*ell+1)*(subtotal+dtotal)/(2*np.pi))}
        result['replay_partial_integrals'].append(entry);print(json.dumps(entry),flush=True)
        save('audit_results.json',result)
    result['energy_replay']=sum(x['energy'] for x in result['replay_partial_integrals'])
    result['original_energy']=json.loads((HERE/'spectral_results.json').read_text())['energy_n8']
    result['replay_difference']=result['energy_replay']-result['original_energy']
    result['elapsed_seconds']=time.perf_counter()-start;result['completed']=True
    result['checks']={'complete_replay':len(result['replay_points'])==41*32,
        'replay_change_below_allowance':abs(result['replay_difference'])<.02,
        'FD_finest_agreement':all(abs(x['discrepancy'])<.002 for x in result['independent_fd'] if x['n']==1000),
        'spline_bound_inputs':all(x['Vmax_upper_spline']<2 and x['I2_upper_spline']<1800 and x['Irad3_upper_spline']<600 for x in result['enclosures'])}
    result['all_checks_pass']=all(result['checks'].values());save('audit_results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['replay_points','replay_partial_integrals']},indent=2),flush=True)
    if not result['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
