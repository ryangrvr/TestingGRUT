"""Leading Wilson channels on supplied torus competitors, not their free energy."""
from pathlib import Path
from fractions import Fraction as F
import json
from loop_paths import coefficient
rows=[]
for n in [324,729,2916]:
 M=n*(n-1)//2;q=2*n;a=F(M-q,10*M)
 square=coefficient(4,a,exact=True);wind=coefficient(9,a,exact=False)
 rows.append(dict(N=n,Q=q,a=str(a),square_coefficient=square['float_coefficient'],square_channels_each_graph=n,thin_short_winding_cycles=n//9,thin_winding_coefficient=wind['float_coefficient'],thin_winding_matrix_element_order=18,balanced_short_winding_length=int(n**.5),balanced_winding_matrix_element_order=2*int(n**.5),scope='One-flux Feshbach vertices at neutral E0. Not total graph energy or Gibbs weights.'))
Path(__file__).with_name('geometry_channels_results.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
