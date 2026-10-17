import sys,json,csv,time,hashlib,platform
from pathlib import Path
from itertools import product
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tests')]
from plane_shadow import *
from independent import projective_count,substitute,rref_key
protocol=json.loads((ROOT/'experiments/protocol.json').read_text())
started=time.perf_counter()
rows=[]; direct=[]; planes=[]
# No random sampling: every cubic-class vector is inspected.
labeled=[]
for cs in product(range(3),repeat=6):
    s,b=signature(cs)
    n=tuple(cs.count(j) for j in range(3))
    expected=0
    if all(k%2==0 for k in n):
        expected=27
        for k in n: expected*=prod(range(1,k,2))
    assert expected==plane_count(cs)
    labeled.append({'classes':cs,'signature':[s,b],'planes':expected})
for p in protocol['base_primes']:
    f=Field(p); g=primitive_root(p)
    for ns in class_representatives():
        cs=tuple(j for j,n in enumerate(ns) for _ in range(n))
        coeffs=[pow(g,c,p) for c in cs]
        t=time.perf_counter(); actual=count_points(f,coeffs)
        predicted=predicted_points(p,1,cs)
        assert actual==predicted,(p,ns,actual,predicted)
        rows.append(dict(p=p,m=1,q=p,modulus='prime',multiplicities=str(ns),actual=actual,predicted=predicted,seconds=time.perf_counter()-t))
        matrices=plane_matrices(coeffs,p)
        keys=[rref_key(A,p) for A in matrices]
        assert len(set(keys))==len(matrices)==plane_count(cs)
        assert all(substitute(coeffs,A,p) for A in matrices)
        planes.append(dict(p=p,multiplicities=ns,planes=len(matrices),substitutions=len(matrices),unique=True))
        if p in protocol['projective_enumeration_primes']:
            t=time.perf_counter(); naive=projective_count(coeffs,p)
            assert naive==actual
            direct.append(dict(p=p,multiplicities=ns,count=naive,seconds=time.perf_counter()-t))
    print('prime',p,'verified',flush=True)
for p,mod in protocol['extensions']:
    f=Field(p,mod); m=f.d; g=primitive_root(p)
    for ns in class_representatives():
        cs=tuple(j for j,n in enumerate(ns) for _ in range(n))
        coeffs=[pow(g,c,p) for c in cs]
        t=time.perf_counter(); actual=count_points(f,coeffs)
        pred=predicted_points(p,m,cs)
        assert actual==pred,(p,m,ns,actual,pred)
        rows.append(dict(p=p,m=m,q=f.q,modulus=str(mod),multiplicities=str(ns),actual=actual,predicted=pred,seconds=time.perf_counter()-t))
    print('extension',p,m,'verified',flush=True)
case=[]
for m in range(1,13):
    n=predicted_pair_points(7,m)
    assert n==predicted_points(7,m,(0,0,0,0,1,2))==predicted_points(7,m,(0,0,1,1,2,2))
    nn=predicted_pair_points(7,2*m); Q=7**m
    assert (n*n+nn-2*(1+Q**4)*n)%(2*Q**2)==0
    case.append({'m':m,'field_size':Q,'points_X':n,'points_Y':n,'planes_X':405 if m%3==0 else 0,'planes_Y':405 if m%3==0 else 27,'lines_both':(n*n+nn-2*(1+Q**4)*n)//(2*Q**2)})
for name,data in [('point_counts',rows),('direct_enumeration',direct)]:
    with (ROOT/f'results/{name}.csv').open('w') as out:
        w=csv.DictWriter(out,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
for name,data in [('labeled_classes',labeled),('plane_checks',planes),('case_study',case)]:
    (ROOT/f'results/{name}.json').write_text(json.dumps(data,indent=2)+'\n')
summary={'protocol_sha256':hashlib.sha256((ROOT/'experiments/protocol.json').read_bytes()).hexdigest(),'python':sys.version,'platform':platform.platform(),'class_vectors':len(labeled),'point_count_cases':len(rows),'projective_enumerations':len(direct),'plane_substitutions':sum(x['substitutions'] for x in planes),'failures':0,'elapsed_seconds':time.perf_counter()-started,'J7':jacobi_prime(7),'A7':quadratic_trace(7)}
(ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
