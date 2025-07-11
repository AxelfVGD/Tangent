load('codecrypto.sage')
load('util.sage')
load('GH78.sage')

# Define the parameters

q = 11
r = 8
m = 3
Fq = GF(q)
Fqm.<a> = Fq.extension(m)

# A boolean variable to decide if the alternant code is Goppa or not
goppa = True
recover_grs = True

# Change this 
n = 250
support = generate_support(n, Fqm)
if goppa:
    Gamma = Fqm['X'].irreducible_element(r)
    multiplier = [1 / Gamma(x) for x in support]
else:
    multiplier = random_multiplier(n, Fqm)

GRS = codes.GeneralizedReedSolomonCode(support, r, multiplier)
C = trace_code(GRS, Fq)
dim = C.dimension()

# Check properness
assert(dim == r*m)

# Check whether the code is square distinguishable or bot
assert(square(C).dimension() < n)

# Create a random nonsigular rm x rm transition matrix
P = GL(dim, Fq).random_element()
D = codes.LinearCode(P * C.generator_matrix()) # D is equal to C

Hpub = D.generator_matrix()

# Compute the algebraic quadratic hull, or "code of relations"
C_rel = Crel(D)

print('Number of quadratic equations:', C_rel.nrows())

# Define the polynomial ring
S = PolynomialRing(Fq,dim,','.join(f'x{i}' for i in range(dim)))
x = S.gens()

Mons2 = vector([x[i]*x[j] for i,j in product(range(dim),repeat=2) if i <= j])

# F is a basis of I_2(Hpub)
F = list(C_rel * Mons2)
JacF = jacobian(F,x)

# complex structure
J = companion_matrix(Fqm.polynomial())
Jr = block_diagonal_matrix([J for _ in range(r)])

# Compute the stabilizing algebra
R = PolynomialRing(Fq,(dim)^2,','.join(f'u_{i}_{j}' for i,j in product(range(dim),repeat=2)))
u = R.gens()
U = matrix([[u[dim*i+j] for j in range(dim)] for i in range(dim)])

print('Computing tangent spaces...')
Eqs = []
i = 0
while len(Eqs) <= (dim)^2:
    g = Hpub.column(i)
    T = codes.LinearCode(JacF(*g).right_kernel())
    print('dim(T) =', T.dimension())
    M = T.parity_check_matrix() * U * T.generator_matrix().transpose()
    Eqs += M.coefficients()
    i += 1
print('done.')
print('Number of tangent spaces: ', i)

print('computing Linear system')
Mac,var = Sequence(Eqs).coefficient_matrix(sparse=False)
print('done.')
Solutions = Mac.right_kernel()

assert(Solutions.dimension() == m)

# compute A at random until we get a generator 
while True:
    A = vector_to_matrix(Solutions.random_element(),dim)
    if A.minimal_polynomial().degree() == m:
        break

R.<x> = Fqm[]
f = A.minimal_polynomial()
rho = f.roots(Fqm)[0][0]
X = rho.polynomial(x)(Jr)

A_Fqm = A.change_ring(Fqm)
X_Fqm = X.change_ring(Fqm)

b,Q = A_Fqm.is_similar(X_Fqm,transformation=True)

Gprime = (Q * Hpub).change_ring(Fq)
Vr = matrix(Fqm,r,n)
for i in range(r):
    for j in range(n):
        Vr[i,j] = Fqm(Gprime.column(j)[m*i:m*(i+1)])

if goppa:
    new_support,new_Gamma = GH78Goppa(Vr)
    new_multiplier = [1 / new_Gamma(x) for  x in new_support]
else:
    new_support, new_multiplier = GH78(Vr)

Test = codes.GeneralizedReedSolomonCode(new_support, r, new_multiplier)

print('FINAL TEST')
E = trace_code(codes.GeneralizedReedSolomonCode(new_support,r,new_multiplier),Fq)
print(f'Attack success: {E.is_subcode(D) and D.is_subcode(E)}')
