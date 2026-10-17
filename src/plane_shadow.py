"""Exact diagonal-cubic invariants. Python standard library only."""
from itertools import combinations, product
from math import prod


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, int(p**0.5)+1))


def primitive_root(p):
    if not is_prime(p):
        raise ValueError('prime required')
    return next(g for g in range(1,p) if len({pow(g,k,p) for k in range(p-1)}) == p-1)


def signature(classes):
    if len(classes) != 6 or any(c not in (0,1,2) for c in classes):
        raise ValueError('six cubic classes required')
    s = sum(classes)%3
    b = tuple(sum((s+sum(classes[i] for i in I))%3 == j
                  for I in combinations(range(6),3)) for j in range(3))
    return s,b


def matchings(items):
    if not items:
        yield ()
        return
    i,*rest = items
    for j in rest:
        for tail in matchings([k for k in rest if k != j]):
            yield ((i,j),)+tail


def plane_count(classes):
    return 27*sum(all(classes[i]==classes[j] for i,j in M)
                  for M in matchings(list(range(6))))


def plane_matrices(coefficients,p):
    """Six-by-three column matrices for every rational plane, without duplicates."""
    out=[]
    for M in matchings(list(range(6))):
        roots=[[r for r in range(p) if (coefficients[i]+coefficients[j]*r**3)%p==0]
               for i,j in M]
        for rs in product(*roots):
            A=[[0]*3 for _ in range(6)]
            for k,((i,j),r) in enumerate(zip(M,rs)):
                A[i][k]=1; A[j][k]=r
            out.append(A)
    return out


def eisenstein_mul(a,b):
    # (a0+a1*z)(b0+b1*z), z^2+z+1=0
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]-a[1]*b[1]


def jacobi_prime(p,g=None):
    if p%3!=1 or not is_prime(p):
        raise ValueError('prime congruent to one modulo three required')
    g=primitive_root(p) if g is None else g
    logs={pow(g,i,p):i%3 for i in range(p-1)}
    h=[0,0,0]
    for x in range(2,p):
        h[(logs[x]+logs[(1-x)%p])%3]+=1
    J=(h[0]-h[2],h[1]-h[2])
    assert J[0]**2-J[0]*J[1]+J[1]**2==p
    return J


def quadratic_trace(p):
    J=jacobi_prime(p)
    K=eisenstein_mul(J,J)
    return p*(2*K[0]-K[1])


def predicted_pair_points(p,m):
    if m<1: raise ValueError('positive extension degree required')
    A=quadratic_trace(p)
    t0,t1=2,A
    for _ in range(2,m+1):
        t0,t1=t1,A*t1-p**4*t0
    Q=p**m
    return sum(Q**j for j in range(5))+(20 if m%3==0 else 2)*Q**2+t1


class Field:
    """Polynomial-basis F_(p^d), using a checked monic irreducible modulus.

    Encoding is little-endian base p. Addition and multiplication have no
    floating point operations. Degrees <=3 suffice for the executed protocol.
    """
    def __init__(self,p,modulus=None):
        if not is_prime(p): raise ValueError('prime characteristic required')
        self.p=p
        self.modulus=tuple(modulus or (0,1))
        self.d=len(self.modulus)-1
        if self.d not in (1,2,3) or self.modulus[-1]!=1:
            raise ValueError('monic modulus of degree one, two or three required')
        if self.d>1 and any(sum(c*pow(x,i,p) for i,c in enumerate(self.modulus))%p==0 for x in range(p)):
            raise ValueError('reducible modulus')
        self.q=p**self.d
        self.digits=[tuple((x//p**i)%p for i in range(self.d)) for x in range(self.q)]
        self.add_table=[[self.add(x,y) for y in range(self.q)] for x in range(self.q)]
    def enc(self,a): return sum(x*self.p**i for i,x in enumerate(a))
    def add(self,x,y):
        return self.enc([(a+b)%self.p for a,b in zip(self.digits[x],self.digits[y])])
    def mul(self,x,y):
        a,b=self.digits[x],self.digits[y]
        c=[0]*(2*self.d-1)
        for i,u in enumerate(a):
            for j,v in enumerate(b): c[i+j]+=u*v
        for k in range(2*self.d-2,self.d-1,-1):
            for j in range(self.d): c[k-self.d+j]-=c[k]*self.modulus[j]
        return self.enc([u%self.p for u in c[:self.d]])
    def power(self,x,n):
        ans=1
        while n:
            if n&1: ans=self.mul(ans,x)
            x=self.mul(x,x); n//=2
        return ans


def count_points(field,coefficients):
    if len(coefficients)!=6 or any(not 0<a<field.q for a in coefficients):
        raise ValueError('six nonzero field elements required')
    cubes=[field.power(x,3) for x in range(field.q)]
    counts=[1]+[0]*(field.q-1)
    for a in coefficients:
        hist=[0]*field.q
        for x in cubes: hist[field.mul(a,x)]+=1
        nonzero=[(y,n) for y,n in enumerate(hist) if n]
        nxt=[0]*field.q
        for x,n in enumerate(counts):
            if n:
                for y,k in nonzero: nxt[field.add_table[x][y]]+=n*k
        counts=nxt
    assert sum(counts)==field.q**6
    assert (counts[0]-1)%(field.q-1)==0
    return (counts[0]-1)//(field.q-1)


def class_representatives():
    seen=set()
    for n0 in range(7):
        for n1 in range(7-n0):
            ns=(n0,n1,6-n0-n1)
            canonical=min(ns,ns[1:]+ns[:1],ns[2:]+ns[:2])
            if canonical not in seen:
                seen.add(canonical)
                yield ns


def predicted_points(p,m,classes):
    if m<1: raise ValueError('positive degree required')
    s,b=signature(tuple(m*c%3 for c in classes))
    K=eisenstein_mul(jacobi_prime(p),jacobi_prime(p))
    alpha=(p*K[0],p*K[1])
    am=(1,0)
    for _ in range(m): am=eisenstein_mul(am,alpha)
    phase=[(1,0),(0,1),(-1,-1)][(-s)%3]
    u,v=eisenstein_mul(am,phase)
    Q=p**m
    assert b[1]==b[2]
    return sum(Q**j for j in range(5))+Q**2*(b[0]-b[1])+2*u-v


def poly_mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out


def middle_polynomial(p,classes):
    s,b=signature(classes)
    J=jacobi_prime(p)
    K=eisenstein_mul(J,J)
    u,v=eisenstein_mul((p*K[0],p*K[1]),[(1,0),(-1,-1),(0,1)][s])
    trace=2*u-v
    P=[1,-trace,p**4]
    for _ in range(b[0]+1): P=poly_mul(P,[1,-p**2])
    for _ in range(b[1]): P=poly_mul(P,[1,p**2,p**4])
    assert len(P)==24
    return P
