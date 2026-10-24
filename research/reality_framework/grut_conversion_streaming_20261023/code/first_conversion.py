#!/usr/bin/env python3
"""First-conversion survival along the exact unconverted redshift trajectory.

Inverse reactions do not enter a first-event survival probability. They do
enter the separate active-population/last-escape calculation.
"""
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from streaming import ROOT,PARENT,H0,OR,OM,OL,hubble,YQ,S0,QOBJ,T0,GS0


def calculate(label,af):
    rows=json.loads((PARENT/'RATE_RESULTS.json').read_text())['rows']
    u=np.array([r['p_over_m'] for r in rows]);sv=np.array([r['sigma_v_conversion_cm3_s'] for r in rows])
    fit=PchipInterpolator(np.log(u),np.log(sv),extrapolate=False)
    u0=np.sqrt(.5);n0=S0*YQ/QOBJ;span=np.log(1/af)
    quadrature={n:np.polynomial.legendre.leggauss(n) for n in [32,64]}
    knots=np.sort(np.log(u0/u));knots=knots[knots>0]
    def rate(y):
        a=af*np.exp(y);p=u0*np.exp(-y)
        val=np.exp(fit(np.log(np.clip(p,u[0],u[-1]))))
        return n0/a**3*val/hubble(a)
    def depth(y,order=32):
        if y==0:return 0.
        edges=np.r_[0.,knots[knots<y],y];left=edges[:-1];right=edges[1:]
        x,w=quadrature[order]
        points=(left[:,None]+right[:,None])/2+(right-left)[:,None]/2*x
        return float(np.sum((right-left)/2*np.sum(rate(points)*w,axis=1)))
    total=depth(span)
    quantiles={}
    for fraction in [.1,.5,.9,.99]:
        target=-np.log1p(-fraction)
        y=brentq(lambda y:depth(y)-target,0,span,xtol=1e-13) if total>=target else None
        quantiles[str(fraction)]=dict(delta_log_a=y,z=float(1/(af*np.exp(y))-1)) if y is not None else None
    # Dense first-event region even when it lies far inside a transport step.
    ys=np.unique(np.r_[0.,np.geomspace(1e-10,.01,85),np.linspace(.01,span,160)])
    optical=np.array([depth(y) for y in ys])
    probes=[1e-7,.001,.1,1.,span]
    error=max(abs(depth(y)/depth(y,64)-1) for y in probes)
    return dict(label=label,a_form=af,initial_conversion_rate_over_H=float(rate(0)),
        optical_quadrature_32_vs_64_max_relative=error,
        first_event_quantiles=quantiles,total_optical_depth_to_today=total,
        first_conversion_probability_today=float(-np.expm1(-total)),
        history=[[float(1/(af*np.exp(y))-1),float(-np.expm1(-d)),float(d)] for y,d in zip(ys,optical)])


def main():
    cases=[('T20MeV',T0/.02*(GS0/10.75)**(1/3)),('T1keV',T0/1e-6),('z1e6',1/1000001)]
    rows=[calculate(n,a) for n,a in cases]
    (ROOT/'FIRST_CONVERSION_RESULTS.json').write_text(json.dumps(dict(
        scope='First conversion of each initial anti-X against fixed cold objects; imposed instantaneous source. Not a permanent-accretion fraction.',rows=rows),indent=2,allow_nan=False)+'\n')
    for r in rows:print(r['label'],'Gamma/H',r['initial_conversion_rate_over_H'],'first median',r['first_event_quantiles']['0.5'])


if __name__=='__main__':main()
