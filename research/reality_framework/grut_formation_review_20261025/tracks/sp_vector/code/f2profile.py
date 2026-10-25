#!/usr/bin/env python3
"""Track sp_vector (own code): window-averaged Higgs norm F_g^2(a) from a vrel history file -> JSON for proca P4.
Usage: python3 f2profile.py hist.txt out.json"""
import json, sys
import numpy as np
HI = 1.1211167304108865e13
h = 1e11 / HI
a, F2, F2X = [], [], []
for line in open(sys.argv[1]):
    if line.startswith('F '):
        _, x, f, fx = line.split(); a.append(float(x)); F2.append(float(f)); F2X.append(float(fx))
a = [1.0] + a; F2 = [1212.98] + F2; F2X = [2 * 23.826064298969072**2] + F2X
json.dump(dict(a=a, F2=F2, F2X=F2X, F2vac=36 * h * h, F2X_end=F2X[-1],
               note='F_g^2 = 2(|X|^2 + 9|Y|^2 + 36|S|^2), physical, H_i units, averaged over log bins (da-weighted)'),
          open(sys.argv[2], 'w'), indent=1)
print(list(zip(a[::4], F2[::4], F2X[::4])))
