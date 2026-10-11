"""Rank control: no new dynamics; same link prescription as supplied probe."""
import json
import numpy as np
from supplied_probe import analyze, four_dimensional_diamond, SEED

n=180
geom=analyze(four_dimensional_diamond(n,SEED),'4D-link')
relations=int(geom['C'].sum())
rng=np.random.default_rng(71011)
block=np.zeros((90,90))
block.flat[rng.choice(block.size,relations,replace=False)]=1
C=np.zeros((n,n)); C[90:,:90]=block
fake=analyze(C,'4D-link')
singular=np.linalg.svd(block,compute_uv=False)
rank=int(np.sum(singular>1e-8))
assert fake['height']==2
assert fake['pair_fraction']==geom['pair_fraction']
assert fake['nonzero_commutator_eigenvalues']==2*rank
assert rank==90
assert fake['nonzero_commutator_eigenvalues']>=geom['nonzero_commutator_eigenvalues']
# Any lower/upper triangular strictly retarded finite matrix has determinant zero.
assert np.linalg.matrix_rank(geom['G'])<n
output={
 'relations':relations,
 'pair_fraction':fake['pair_fraction'],
 'sprinkled_height':geom['height'],
 'sprinkled_commutator_rank':geom['nonzero_commutator_eigenvalues'],
 'random_layered_height':fake['height'],
 'random_layered_block_rank':rank,
 'random_layered_commutator_rank':fake['nonzero_commutator_eigenvalues'],
 'random_layered_positive_modes':fake['sj_positive_modes'],
 'minimum_block_singular_value':float(singular.min()),
 'verdict':'PASS: high mode rank does not select manifoldlike orders'
}
print(json.dumps(output,indent=2))
