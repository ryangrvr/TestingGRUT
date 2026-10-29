"""Recreate the parent 10% pulse and same-lattice stationary control."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy import sparse
from scipy.sparse.linalg import spsolve

ROOT=Path(__file__).resolve().parents[1]
M1=ROOT.parent/'grut_matched_candidate_20261027'
PROFILE=np.loadtxt(M1/'PROFILE_M1.csv',delimiter=',',skiprows=1)
OMEGA=json.loads((M1/'DEFORMATION_RESULTS.json').read_text())['primary']['Omega']
F=CubicHermiteSpline(PROFILE[:,0],PROFILE[:,1],PROFILE[:,2])
S=CubicHermiteSpline(PROFILE[:,0],PROFILE[:,3],PROFILE[:,4])

def geometry(h,L):
    r=(np.arange(round(L/h))+.5)*h;v=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    face=4*np.pi*(r+h/2)**2/h;face[-1]*=2;left=np.r_[0,face[:-1]]
    lap=sparse.diags([face[:-1]/v[:-1],-(face+left)/v,left[1:]/v[1:]],[1,0,-1],format='csr')
    return r,v,face,lap

def values(r):
    a=np.minimum(r,160.);return np.where(r>160,0,F(a)),np.where(r>160,1,S(a))

def stationary(h,L,q):
    r,v,face,lap=geometry(h,L);n=len(r);f,s=values(r);om=OMEGA
    def residual(f,s,om):
        rs=-lap@s+(3*(s*s-1)+f*f)*s-.6*(s-1)**2;rs[-1]-=face[-1]/v[-1]
        return np.r_[-lap@f+(s*s+.1*f*f-om*om)*f,rs,om*(v@f**2)/(.01*q)-1]
    for it in range(20):
        res=residual(f,s,om)
        if np.max(abs(res))<2e-10:break
        ff=-lap+sparse.diags(s*s+.3*f*f-om*om);ss=-lap+sparse.diags(9*s*s-3+f*f-1.2*(s-1));fs=sparse.diags(2*f*s)
        jac=sparse.bmat([[ff,fs,sparse.csr_matrix((-2*om*f)[:,None])],[fs,ss,sparse.csr_matrix((n,1))],
                         [sparse.csr_matrix((2*om*v*f/(.01*q))[None,:]),sparse.csr_matrix((1,n)),sparse.csr_matrix([[v@f**2/(.01*q)]])]],format='csr')
        delta=spsolve(jac,-res);alpha=1.
        for _ in range(15):
            nf=f+alpha*delta[:n];ns=s+alpha*delta[n:2*n];no=om+alpha*delta[-1]
            if np.linalg.norm(residual(nf,ns,no))<np.linalg.norm(res):break
            alpha*=.5
        else:raise RuntimeError('stationary Newton line search failed')
        f,s,om=nf,ns,no
    err=float(max(abs(residual(f,s,om))));assert err<2e-9,err
    return r,f,s,-om*f,{'omega':float(om),'max_residual':err,'Q':q,'iterations':it}

def main():
    p=argparse.ArgumentParser();p.add_argument('--h',type=float,required=True);p.add_argument('--L',type=float,default=200.)
    p.add_argument('--type',choices=['pulse','stationary'],default='pulse');p.add_argument('--tag',required=True);a=p.parse_args()
    r,v,_,_=geometry(a.h,a.L)
    if a.type=='pulse':
        f,s=values(r);shape=np.exp(-(r/12)**2);f*=1+.1*shape;s+=.1*(1-s)*shape
        py=-OMEGA*f;py*=2000/np.sum(-f*py*v);meta={'omega_reference':OMEGA,'amplitude':.1}
    else:r,f,s,py,meta=stationary(a.h,a.L,200000.)
    data=np.c_[r,f,0*r,s,0*r,py,0*r];path=ROOT/'.work'/(a.tag+'.csv')
    np.savetxt(path,data,delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
    meta.update({'h':a.h,'L':a.L,'type':a.type,'charge':float(np.sum(-f*py*v)/.01)})
    (ROOT/'.work'/(a.tag+'.json')).write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
