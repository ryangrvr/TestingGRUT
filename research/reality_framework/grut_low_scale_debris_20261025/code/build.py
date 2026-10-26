"""Generate the published DOP853 coefficients and compile the C integrator.

Coefficient tables: scipy.integrate._ivp.dop853_coefficients, SciPy 1.17.0,
following Hairer's DOP853 method. This is an original C implementation of
that numerical method; no reviewer integrator source is copied.
"""
import os
import subprocess
from pathlib import Path
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as dc

HERE=Path(__file__).resolve().parent
WORK=Path(os.environ.get('GRUT_LOW_DEBRIS_WORKDIR','/tmp/grut_low_scale_debris_20261025'))
BUILD=WORK/'build'


def build():
    BUILD.mkdir(parents=True,exist_ok=True)
    def vector(name, a):
        a=np.asarray(a)
        if a.ndim==1:
            return 'static const double '+name+'['+str(len(a))+'] = {'+','.join(format(x,'.17g') for x in a)+'};\n'
        return 'static const double '+name+'['+str(len(a))+']['+str(a.shape[1])+'] = {\n'+',\n'.join('{'+','.join(format(x,'.17g') for x in row)+'}' for row in a)+'};\n'
    header='/* DOP853 coefficients, generated from SciPy; see build.py. */\n'
    for name,a in [('AA',dc.A[:12,:12]),('CC',dc.C[:12]),('BB',dc.B),('EE3',dc.E3),('EE5',dc.E5)]:
        header+=vector(name,a)
    (BUILD/'dop853.h').write_text(header)
    exe=BUILD/'release'
    subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(BUILD),str(HERE/'release.c'),'-lm','-o',str(exe)],check=True)
    return exe


if __name__=='__main__':print(build())
