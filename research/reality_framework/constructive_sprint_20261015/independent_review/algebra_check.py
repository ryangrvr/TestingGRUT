"""Small independent determinant and square-free polynomial checks."""
from pathlib import Path
import importlib.util
from fractions import Fraction as F
from itertools import permutations
from math import factorial
import json
import numpy as np
from scipy.linalg import expm

SOURCE = Path(__file__).resolve().parent.parent/'charged_access'/'cycle_cover.py'
spec=importlib.util.spec_from_file_location('review_cover',SOURCE)
cc=importlib.util.module_from_spec(spec);spec.loader.exec_module(cc)

def add(a,b):
    c=a.copy()
    for k,v in b.items():c[k]=c.get(k,F(0))+v
    return {k:v for k,v in c.items() if v}
def mul(a,b):
    c={}
    for i,x in a.items():
        for j,y in b.items():
            if i&j:continue
            c[i|j]=c.get(i|j,F(0))+x*y
    return {k:v for k,v in c.items() if v}
def scale(a,c):return {k:c*v for k,v in a.items() if c*v}
def mm(a,b):
    n=len(a);return [[sum_poly(mul(a[i][k],b[k][j]) for k in range(n)) for j in range(n)] for i in range(n)]
def sum_poly(xs):
    ans={}
    for x in xs:ans=add(ans,x)
    return ans
def det(a):
    n=len(a);ans={}
    for p in permutations(range(n)):
        term={0:F((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))}
        for i in range(n):term=mul(term,a[i][p[i]])
        ans=add(ans,term)
    return ans

def squarefree(n,arcs,S):
    q=len(arcs);A=[[{} for j in range(n)] for i in range(n)]
    for e,(i,j) in enumerate(arcs):A[i][j]={1<<e:F(1)}
    power=[[{0:F(int(i==j))} if i==j else {} for j in range(n)] for i in range(n)]
    cosh=[[{} for j in range(n)] for i in range(n)]
    tanh=[[{} for j in range(n)] for i in range(n)]
    # Independently divide sinh(sqrt(x))/sqrt(x) by cosh(sqrt(x));
    # do not reuse the cover implementation's logarithmic coefficients.
    u=[]
    for r in range(q+1):
        u.append(F(1,factorial(2*r+1))-
                 sum(F(1,factorial(2*k))*u[r-k] for k in range(1,r+1)))
    for r in range(q+1):
        for i in range(n):
            for j in range(n):
                cosh[i][j]=add(cosh[i][j],scale(power[i][j],F(1,factorial(2*r))))
                tanh[i][j]=add(tanh[i][j],scale(power[i][j],u[r]))
        power=mm(power,A)
    minor=[[tanh[i][j] for j in S] for i in S]
    result=mul(det(cosh),det(minor))
    return result.get((1<<q)-1,F(0))

def run():
    rng=np.random.default_rng(761);identity=[]
    n=4;es=[(0,1),(1,2),(2,3),(0,3)];q=len(es)
    for k in (1,2,3):
        B=np.zeros((n,q),complex)
        for e,(i,j) in enumerate(es):
            B[i,e]=np.exp(1j*rng.uniform(-3,3));B[j,e]=-np.exp(1j*rng.uniform(-3,3))
        h=np.block([[np.zeros((n,n)),B],[B.conj().T,np.zeros((q,q))]])
        T=.173;U=expm(-1j*T*h);C=B@B.conj().T
        vals,V=np.linalg.eigh(C);root=np.sqrt(np.maximum(vals,0))
        cos=V@np.diag(np.cos(T*root))@V.conj().T
        f=np.array([np.tan(T*x)/x if x>1e-12 else T for x in root])
        mat=(V@np.diag(f)@V.conj().T)@B
        S=list(range(k));H=list(range(k));rows=S+[n+e for e in range(q) if e not in H]
        direct=np.linalg.det(U[np.ix_(rows,range(n,n+q))])
        formula=(-1j)**k*np.linalg.det(cos)*np.linalg.det(mat[np.ix_(S,H)])
        error=abs(abs(direct)-abs(formula));assert error<1e-13
        identity.append({'holes':k,'absolute_magnitude_difference':float(error)})
    covers=[]
    for n,arcs in [(3,[(0,1),(1,2),(2,0)]),(4,[(0,1),(1,2),(2,3),(3,0)])]:
        for S in [(),(0,),(0,1),(0,2)]:
            poly=squarefree(n,arcs,S);ports=cc.cover_coefficient(n,arcs,S)
            assert poly==ports
            covers.append({'N':n,'marked':S,'coefficient':str(poly)})
    result={'status':'PASS','slater_row_replacement_identity':identity,'independent_squarefree_determinants':covers,
      'scope':'Independent small algebra checks; not a re-execution of full charged graph-orbit or Hamiltonian evolution runs.'}
    (Path(__file__).resolve().parent/'ALGEBRA_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
