"""Deterministic arithmetic checks separate from the mathematical formulas."""
from pathlib import Path
import sys,random,json
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'src'))
from plane_shadow import Field
rng=random.Random(602)
checks=0
for p,mod in [(7,None),(7,[1,0,1]),(7,[2,0,0,1]),(13,[2,0,1]),(19,[1,0,1])]:
 f=Field(p,mod)
 for a in range(1,f.q):
  assert f.power(a,f.q-1)==1
  checks+=1
 for _ in range(500):
  a,b,c=[rng.randrange(f.q) for _ in range(3)]
  assert f.mul(a,f.add(b,c))==f.add(f.mul(a,b),f.mul(a,c))
  assert f.mul(f.mul(a,b),c)==f.mul(a,f.mul(b,c))
  assert f.add_table[a][b]==f.add(a,b)
  checks+=3
try:
 Field(7,[0,0,1])
 raise AssertionError('reducible modulus accepted')
except ValueError:checks+=1
(R/'results/field_checks.json').write_text(json.dumps({'seed':602,'assertions':checks,'failures':0},indent=2)+'\n')
print(checks,'field arithmetic checks passed')
