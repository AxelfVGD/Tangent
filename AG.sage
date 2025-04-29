load('codecrypto.sage')
load('util.sage')
load('GH78.sage')

# Define the parameters

q = 5
m = 4
k = 6

recover_ag = False

def vector_to_matrix(v,n):
    res = matrix(Fq,n)
    count = 0
    for i in range(n):
        for j in range(n):
            res[i,j] = v[count]
            count += 1
    return res

Fq = GF(q)
Fqm.<a> = Fq.extension(m)

P.<x,y,z> = ProjectiveSpace(Fqm,2)

# Change this
X = Curve(y^2*z - x^3 - x*z^2 + z^3)

assert(X.is_smooth())
g = X.genus()
print('genus(X) =', g)
r = k+1-g

F = X.function_field()
pls = F.places()

p_inf = X([0,1,0])
P_inf, = p_inf.places()
pls.remove(P_inf)

# Q = pls[randint(0,len(pls))]
# pls.remove(Q)
E = k*P_inf

AG = codes.EvaluationAGCode(pls, E)
n = AG.length()

multiplier = random_multiplier(n, Fqm)
AG = codes.LinearCode(AG.generator_matrix() * diagonal_matrix(multiplier))

AG_rel = Crel(AG)
print(AG_rel.nrows())

C = trace_code(AG, Fq)

print(C)
dim = C.dimension()
assert(dim == (k+1-g) * m)

P = GL(dim, Fq).random_element()
D = codes.LinearCode(P * C.generator_matrix())

G = D.generator_matrix()

C_rel = Crel(D)

N_rel = C_rel.nrows()
print(N_rel)
print(square(D))
assert(N_rel > binomial(dim+1,2)-n)

R = PolynomialRing(Fq,dim,','.join(f'x{i}' for i in range(dim)))
x = R.gens()
Mons2 = vector([x[i]*x[j] for i,j in product(range(dim),repeat=2) if i <= j])

F = list(C_rel * Mons2)
JacF = jacobian(F,x)

# complex structure
J = companion_matrix(Fqm.polynomial())
Jr = block_diagonal_matrix([J for _ in range(r)])

S = PolynomialRing(Fq,(dim)^2,','.join(f'u_{i}_{j}' for i,j in product(range(dim),repeat=2)))
u = S.gens()
U = matrix([[u[dim*i+j] for j in range(dim)] for i in range(dim)])

print('Computing tangent spaces...')
Eqs = []
i = 0
while len(Eqs) <= 2*(dim)^2:
    g = G.column(i)
    T = codes.LinearCode(JacF(*g).right_kernel())
    print('dim(T) =', T.dimension())
    M = T.generator_matrix() * U * T.parity_check_matrix().transpose()
    Eqs += M.coefficients()
    i += 1
print('done.')
print('Number of tangent spaces: ', i)

print('computing Linear system')
Mac,var = Sequence(Eqs).coefficient_matrix(sparse=False)
print('done.')
Solutions = Mac.right_kernel()
print('dimension =', Solutions.dimension())

while True:
    A = vector_to_matrix(Solutions.random_element(),dim).transpose()
    if A.minimal_polynomial().degree() == m:
        break


R.<x> = Fqm[]
f = A.minimal_polynomial()
rho = f.roots(Fqm)[0][0]
X = rho.polynomial(x)(Jr)

A_Fqm = A.change_ring(Fqm)
X_Fqm = X.change_ring(Fqm)

b,Q = A_Fqm.is_similar(X_Fqm,transformation=True)

Gprime = Q * G 
Vr = matrix(Fqm,r,n)
for i in range(r):
    for j in range(n):
        Vr[i,j] = sum(Gprime[m*i+k,j]*a^k for k in range(m))

E = codes.LinearCode(Vr)
print(E)

if recover_ag:
    for j in range(m):
        print('testing j =', j)
        b = True 
        for c in AG.generator_matrix():
            c = vector([ci^(q^j) for ci in c])
            b = b and (c in E)
        print(b)
        if b:
            break

else:
    Test = trace_code(E,Fq)
    print(f'Attack success: {Test == D}')