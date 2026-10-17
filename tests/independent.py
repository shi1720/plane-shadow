"""Independent oracles: projective enumeration and polynomial substitution."""
from itertools import product


def projective_count(coefficients,p):
    # Last nonzero coordinate equals one. Each projective point occurs once.
    cubes=[x**3%p for x in range(p)]
    total=0
    for last in range(6):
        for prefix in product(range(p),repeat=last):
            total += (coefficients[last]+sum(a*cubes[x] for a,x in zip(coefficients,prefix)))%p==0
    return total


def substitute(coefficients,matrix,p):
    # Ordered tensor expansion, deliberately not the pairing formula.
    coeff={}
    for a,row in zip(coefficients,matrix):
        for i,j,k in product(range(3),repeat=3):
            exponent=tuple([i,j,k].count(t) for t in range(3))
            coeff[exponent]=(coeff.get(exponent,0)+a*row[i]*row[j]*row[k])%p
    return all(v==0 for v in coeff.values())


def rref_key(matrix,p):
    A=[list(row) for row in zip(*matrix)]
    r=0
    for c in range(6):
        pivot=next((i for i in range(r,3) if A[i][c]),None)
        if pivot is None: continue
        A[r],A[pivot]=A[pivot],A[r]
        u=pow(A[r][c],-1,p)
        A[r]=[x*u%p for x in A[r]]
        for i in range(3):
            if i!=r:
                u=A[i][c]
                A[i]=[(x-u*y)%p for x,y in zip(A[i],A[r])]
        r+=1
        if r==3: break
    assert r==3
    return tuple(tuple(row) for row in A)
