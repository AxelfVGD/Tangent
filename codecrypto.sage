from itertools import product


# GENERAL CODING STUFF

# C \cap D
def inter(C,D):
    return codes.LinearCode(matrix(list(C.parity_check_matrix()) + list(D.parity_check_matrix()))).dual_code()

# C * D
def code_product(C,D):
    GC = C.generator_matrix()
    GD = D.generator_matrix()
    return codes.LinearCode(matrix([GC[i].pairwise_product(GD[j]) for i,j in product(range(GC.nrows()), range(GD.nrows()))]))

# C * C
def square(C):
    G = C.generator_matrix()
    # C = codes.LinearCode(G)
    return codes.LinearCode(matrix([G[i].pairwise_product(G[j]) for i,j in product(range(G.nrows()), repeat=2) if i <= j]))

# X * C \subset D
def cond(C,D):
    return code_product(C,D.dual_code()).dual_code()

def shorten(C, I):
    G = C.generator_matrix()
    k = G.nrows()
    n = G.ncols()
    rows = [j for j in range(k) if not j in I]
    cols = [j for j in range(n) if not j in I]
    return codes.LinearCode(G[rows,cols])

# If F is a prime field, build the "trace matrix" easily
# if not, use "SubfieldSubcode" method
# NB: the prime field version works also in the nonprime field case but is harder 
# to code in sagemath
def trace_code(C, F):
    if F.is_prime_field():
        G = C.generator_matrix()
        m = G.base_ring().degree() / F.degree()
        k = G.nrows()
        n = G.ncols()

        TrG = matrix(Fq, m*k, n)
        count = 0
        for i in range(k):
            for j in range(n):
                c = list(G[i,j])
                for u in range(m):
                    TrG[m*i+u,j] = c[u]
                    count += 1
        return codes.LinearCode(TrG)

    else:
        return codes.SubfieldSubcode(C.dual_code(),F).dual_code()


# GENERALIZED REED SOLOMON CODES

# generate a random support for a GRS code 
def generate_support(n, F):
    if n > F.cardinality():
        exit('Error')
    return Permutations(F).random_element()[:n]

def random_multiplier(n, F):
    L = list(F)
    assert(L[0] == 0)
    return [L[randint(1,len(L)-1)] for _ in range(n)]

def Crel(C):
    G = C.generator_matrix()
    G2 = matrix([G[i].pairwise_product(G[j]) for i,j in product(range(G.nrows()), repeat=2) if i <= j])
    return G2.left_kernel().basis_matrix()

def formula_square_code(r, m, q):
    if r < q-1:
        return min(n, binomial(r*m+1,2)-m*binomial(r-1,2))
    e = ceil(log(r / (q-1)**2,q))+1
    print(e)
    # return min(n, binomial(r*m+1,2)-r*m/2*((2*e+1)*r-2*(q-1)*q**(e-1)-1))
    return r*m/2*((2*e+1)*r-2*(q-1)*q**(e-1)-1)

def hilbert_random(d, k, n):
    s = n-k
    return max([0,sum((-1)^i / (k+d-i-1)*binomial(s,i)*binomial(k+d-i-1,d-i+1)*binomial(k+d-i-1,d-i) for i in range(d+1))])

def relation_to_matrix(rel,k,F):
    assert(len(rel) == binomial(k+1,2))
    res = matrix(F,k,k)
    c = 0
    for i,j in product(range(k), repeat=2):
        if i < j:
            res[i,j] = rel[c]
            res[j,i] = rel[c]
            c += 1
        elif i == j:
            res[i,i] = 2*rel[c]
            c += 1
    return res

def vectorize(M):
    return vector(flatten([list(r) for r in M.rows()]))

# ALGEBRAIC GEOMETRY

# Given an ideal I and a point P in V(I), compute the dimension of the tangent space at P

def compute_local_dimension(F,P):
    # check that P is in V(I)
    for f in F:
        if f(*P) != 0:
            exit(f'Error: {P} is not in ({F})')
    X = F[0].parent().gens()
    J = jacobian(F,X)(*P)
    return len(X) - J.rank()

def rational_normal_curve(Points,BasePoints):

    '''
    Griffiths, P. and Harris, J. in Principles of Algebraic Geometry:
    "Through any n+3 points in general position in P^n there passes a unique rational normal curve."
    This algorithm implements there method for retrieving the parametrization of such curve
    '''

    def find_homography(a,b,c):
        # a,b and c must be points on the projective line
        return matrix([[0,a[1],0,-a[0]],[b[1],b[1],-b[0],-b[0]],[c[1],0,-c[0],0]]).right_kernel().basis_matrix()[0]

    Fq = Points[0].base_ring()
    P1.<u,v> = ProjectiveSpace(Fq,1)
    r = len(Points[0])-1
    P = ProjectiveSpace(Fq,r)

    # H[i] parametrizes the hyperplanes passing through Points[0],...,Points[i-1],Points[i+1],...,Points[r-1]
    H = []
    for i in range(r):
        M = matrix([list(Points[j]) for j in range(r) if j != i])
        plane = vector([u,v]) * M.right_kernel().basis_matrix()

        # we reparametrize so that H[i](0:1) = BasePoints[0], H[i](1:1) = BasePoints[1] and H[i](1:0) = BasePoints[2]
        lambda_j = []
        for p in BasePoints:
            f = vector(p) * plane
            lambda_j.append(P1([f.coefficient(v),-f.coefficient(u)]))
        # Find the homography that sens (0,1,\inf) onto (lambda_1,lambda_2,lambda_3)    
        a,b,c,d = find_homography(*lambda_j)
        H.append(plane(a*u+b*v,c*u+d*v))

    H = matrix(H)
    F = FractionField(H.base_ring())
    H = H.change_ring(F)
    Psi = H.right_kernel().basis_matrix()[0]
    return Psi

