"""Independent finite-support proof controls for the partial HB kit.

If Q=A#P+b with A in GL(k), it is in the same affine orbit. For nonsingular
covariances the whitening matrices give an orthogonal residual:
O=Sigma_Q^(-1/2) A Sigma_P^(1/2), OO^T=I. Thus d_q^BL=0, and one common
reference makes epsilon exactly zero. This is a theorem, not a small residual.
Singular A is not in GL(k), even if it happens to be injective on a finite support.
These finite fixtures do not replace CR-3's unbuilt HB families or card gates.
"""
from fractions import Fraction as F


def rational(x):
    try: return F(x)
    except (ValueError,TypeError,OverflowError,ZeroDivisionError) as exc:
        raise ValueError('finite rational data required') from exc


def determinant(matrix):
    a=[list(map(rational,row)) for row in matrix]
    n=len(a)
    if n==0 or any(len(row)!=n for row in a): raise ValueError('nonempty square matrix required')
    sign=1; value=F(1)
    for col in range(n):
        pivot=next((i for i in range(col,n) if a[i][col]),None)
        if pivot is None: return F(0)
        if pivot!=col: a[col],a[pivot]=a[pivot],a[col]; sign=-sign
        p=a[col][col]; value*=p
        for i in range(col+1,n):
            factor=a[i][col]/p
            for j in range(col,n): a[i][j]-=factor*a[col][j]
    return sign*value


def law(points,weights):
    if len(points)!=len(weights) or not points: raise ValueError('support/weight mismatch or empty law')
    atoms=[tuple(map(rational,p)) for p in points]
    k=len(atoms[0]); ws=tuple(map(rational,weights))
    if not k or any(len(p)!=k for p in atoms): raise ValueError('ragged/empty record chart')
    if any(w<0 for w in ws) or sum(ws)<=0: raise ValueError('invalid probability measure')
    total=sum(ws); out={}
    for p,w in zip(atoms,ws):
        if w: out[p]=out.get(p,F(0))+w/total
    return out


def covariance(p):
    k=len(next(iter(p))); mean=[sum(w*x[i] for x,w in p.items()) for i in range(k)]
    return [[sum(w*(x[i]-mean[i])*(x[j]-mean[j]) for x,w in p.items())
             for j in range(k)] for i in range(k)]


def push(p,A,b):
    k=len(next(iter(p)))
    if len(A)!=k or len(b)!=k or any(len(row)!=k for row in A): raise ValueError('map/chart mismatch')
    A=[list(map(rational,row)) for row in A]; b=list(map(rational,b))
    out={}
    for x,w in p.items():
        y=tuple(sum(A[i][j]*x[j] for j in range(k))+b[i] for i in range(k))
        out[y]=out.get(y,F(0))+w
    return out


def affine_certificate(p,q,A,b,source_domain=None):
    if source_domain is not None and not set(p)<=set(source_domain):
        raise ValueError('calibration domain misses source support')
    if determinant(A)==0: raise ValueError('singular interface is not a GL equivalence')
    if determinant(covariance(p))<=0 or determinant(covariance(q))<=0:
        raise ValueError('finite nonsingular covariance required')
    if push(p,A,b)!=q: return False
    return True


def skew_square(p):
    if any(len(x)!=1 for x in p): raise ValueError('scalar witness only')
    mu=sum(w*x[0] for x,w in p.items())
    m2=sum(w*(x[0]-mu)**2 for x,w in p.items())
    m3=sum(w*(x[0]-mu)**3 for x,w in p.items())
    if m2<=0: raise ValueError('zero variance')
    return m3*m3/(m2*m2*m2)


def scalar_control(laws):
    # Unequal squares certify lack of any scalar affine/reflection equivalence.
    # Equal moments alone do NOT certify orbit equality.
    values=[skew_square(p) for p in laws]
    return {'gamma1_squared':values,'nonzero_witness':len(set(values))>1}
