#!/usr/bin/env python3
"""Check every compiled source generator against an independent explicit average.
The expected coefficients use literal S4 composition, transpose of the seed,
global bosonic sandwich, and all six ordered physical-site permutations.
"""
import sys,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
root=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/data/werner_review/werner-problem6-proof')
sys.path.insert(0,str(root/'certificate'))
from basis_values import compile_basis
from exact_coordinates import CATALOG
P=list(itertools.permutations(range(4)));index={p:i for i,p in enumerate(P)}
perm3=list(itertools.permutations(range(3)))
def compose(g,p,h):return tuple(g[p[h[i]]] for i in range(4))
def inv(p):return tuple(p.index(i) for i in range(4))
counts={}
with np.load(root/'certificate/data/basis_data.npz',allow_pickle=False) as z:
 for f in ('31','22','13'):
  ix=z['basis'+f+'_indices'];compiled=compile_basis(f,ix)
  group=[p+(3,) for p in perm3] if f!='22' else [p+q for p in ((0,1),(1,0)) for q in ((2,3),(3,2))]
  den=12*len(group)**2
  for k,j in enumerate(ix):
   seed=tuple(P[x] for x in CATALOG[int(j)]);expected=Counter()
   for transpose in (False,True):
    words=tuple(inv(p) for p in seed) if transpose else seed
    for g in group:
     for h in group:
      a=tuple(index[compose(g,p,h)] for p in words)
      for s in perm3:
       b=tuple(a[i] for i in s)+(a[3],)
       flat=((b[0]*24+b[1])*24+b[2])*24+b[3];expected[flat]+=1
   lo=compiled.row_starts[k];hi=compiled.row_starts[k+1] if k+1<len(ix) else len(compiled.table_indices)
   actual=Counter()
   for x,c in zip(compiled.table_indices[lo:hi],compiled.coefficients[lo:hi]):actual[int(x)]+=int(c)
   assert compiled.row_denominators[k]==den and expected==actual,(f,k)
  counts[f]=len(ix)
result={'status':'INDEPENDENT_LITERAL_PROJECTOR_EXPANSIONS_PASS','generators_checked':counts,'total':sum(counts.values()),
 'method':'Independent literal full group average, no use of projected_generator to construct expected values',
 'includes':['Hermitian transpose of seed','both global bosonic projectors','all ordered physical-site permutations','orbit-average multiplicities']}
(Path(__file__).parent/'projector_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
