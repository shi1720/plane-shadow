import sys,json
from pathlib import Path
from collections import defaultdict
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from plane_shadow import *
allrows=[]
for p in [7,13,19,31,37,43,61,67,73,79,97]:
    groups=defaultdict(list)
    for ns in class_representatives():
        cs=tuple(i for i,n in enumerate(ns) for _ in range(n))
        poly=middle_polynomial(p,cs)
        groups[tuple(poly)].append(ns)
        allrows.append({'p':p,'multiplicities':ns,'signature':signature(cs),'planes':plane_count(cs),'middle_polynomial':poly})
    assert len(groups)==9,(p,len(groups))
    assert next(v for v in groups.values() if len(v)>1)==[(1,1,4),(2,2,2)]
(R/'results/zeta_classification.json').write_text(json.dumps(allrows,indent=2)+'\n')
print('Exactly nine zeta functions in ten classes at all eleven tested primes.')
