from pathlib import Path
import itertools,json
from check import run
rows=[]
for n in [3,4,5,6,7,8]:
 edges=list(itertools.combinations(range(n),2));ei={e:i for i,e in enumerate(edges)}
 mask=sum(1<<ei[tuple(sorted((i,(i+1)%n)))] for i in range(n))
 row=run(n,n,2,mask);row['energy_per_vertex']=row['energies'][0]/n
 rows.append(row)
Path(__file__).with_name('cycle_results.json').write_text(json.dumps(rows,indent=2)+'\n')
print('CP12_CONDITIONAL_CYCLE_CONTROLS_COMPLETE')
