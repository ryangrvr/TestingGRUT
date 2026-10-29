"""Check the quadrature of the CLASS source distribution and its derivative."""
from pathlib import Path
import json,numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parents[1]
def audit(name):
 q,f=np.loadtxt(ROOT/name,unpack=True);F=CubicSpline(q,f);rows=[]
 for N in (30,60,120,200,400,800):
  x=4*np.arange(1,N+1)/N;w=np.full(N,4/N);w[-1]*=.5
  fq=F(x);norm=w@(x*x*fq);deriv=-w@(x**3*F(x,1))/(3*norm);q2=(w@(x**4*fq))/norm
  rows.append({'bins':N,'qmax':4,'number_integral':float(norm),'adiabatic_derivative_identity':float(deriv),'q2_mean':float(q2),'min_interpolated_f':float(min(fq))})
 return rows
out={'discrete_smoothed_96':audit('psd_96.dat'),'moment_matched_maxwell':audit('psd_maxwell96.dat'),
 'interpretation':'For a sufficiently decaying smooth f(q), integration by parts gives -integral(q^3 f_prime)/(3 integral(q^2 f))=1. This is a diagnostic of quadrature consistency, not a proof of CLASS accuracy; finite qmax produces a boundary term.'}
(ROOT/'psd_audit.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
