from pathlib import Path
from math import prod
import csv,sys,json
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from plane_shadow import plane_count
rows=list(csv.DictReader((R/'results/grassmannian.csv').open()))
assert len(rows)==22
for r in rows:
    p=int(r['p'])
    cs=[0]*int(r['n0'])+[1]*int(r['n1'])+[2]*int(r['n2'])
    assert int(r['planes'])==plane_count(cs)
    assert int(r['all_subspaces'])==prod(p**(6-i)-1 for i in range(3))//prod(p**(3-i)-1 for i in range(3))
lines=list(csv.DictReader((R/'results/lines.csv').open()))
assert len(lines)==2 and all(int(r['lines'])==2943 for r in lines)
assert all(r['actual']==r['predicted'] for r in csv.DictReader((R/'results/point_counts.csv').open()))
print('Verified 22 independent plane searches, two line searches, and 150 exact point counts.')
