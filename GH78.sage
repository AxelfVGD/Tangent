load('codecrypto.sage')


# Griffiths, P. and Harris, J. in Principles of Algebraic Geometry:
# "Through any n+3 points in general position in P^n there passes a unique rational normal curve."
# Finding the parametrization enables one to recover a support and a multiplier of a given GRS code


# given 3 points a,b,c in Fq, find a homography s such that s(0) = a, s(1) = b and s(\infinity) = c
def find_homography(L):
    a,b,c = L
    M = matrix([[0,1,0,-a],[1,1,-b,-b],[1,0,-c,0]])
    return M.right_kernel().basis_matrix()[0]

# Given G a generator matrix of C = GRS_k(x,y), recovers x' and y' such that C = GRS_k(x',y')
def GH78(G):
    n = G.ncols()
    k = G.nrows() - 1  # the columns of G are identified with points in P^(k-1)
    Fq = G.base_ring()

    # ensure we can scale the points so as to have them all in the affine space A^(k-1) = P^(k-1) \ {x_k = 0}
    while 0 in G[-1]:
        G = GL(k+1,Fq).random_element() * G
    
    Points = []
    for i in range(k+3):
        g = G.column(i)
        Points.append([g[i] / g[-1] for i in range(k)])

    R = PolynomialRing(Fq, k+1, ','.join(f'X{i}' for i in range(1,k+1)) + ',t')
    X = list(R.gens())
    t = X.pop()


    # List the parametrization of each hyperplane
    # H[i] parametrizes the hyperplanes containing Points[0], ... , Points[n-1] but Points[i]
    # We also rearrange the parametrization so that Points[n] \in H[i](t = 0), Points[n+1] \in H[i](t = 1) and Points[n+2] \in H[i](t = \infinity)
    H = []

    for i in range(k):
        P = Points[(i+1)%k]
        M = matrix([list(Points[j]) + [1] for j in range(k) if j != i])

        B = M.right_kernel().basis_matrix()
        H.append(vector([t,1]) * B * vector([X[i] for i in range(k)] + [1]))

        lambda_j = []
        for j in range(3):
            var = list(Points[k+j]) + [t]
            fj = H[i](*var)
            lambda_j.append(Fq( - fj(*[0 for _ in range(k+1)]) / fj.coefficient(t)))

        a,b,c,d = find_homography(lambda_j)
        H[i] = H[i].subs(t=(a*t+b)/(c*t+d)).numerator()

    # The intersection between all hyperplanes defines a paramtrized curve of parameter t, so we solve for X_0, ... , X_{k-1}
    Mat = matrix([[H[i].coefficient(X[j]) for j in range(k)] for i in range(k)])

    A.<T> = Fq['T']

    Mat = Mat.subs(t=T)
    det = Mat.det()

    Target = []
    for i in range(k):
        var = [0 for _ in range(k)] + [t]
        Target.append(-H[i](*var).subs(t=T))

    # Psi is a parametrization of a rational normal curve
    Psi = Mat.inverse() * vector(Target)

    # Psi(t0) is the point at infinity
    t0 = det.roots()

    if len(det.roots()) == 1:
        t0 = det.roots()[0][0]
        Psi = [psi.subs(T=t0+1/T) for psi in Psi]

    new_support = []

    for j in range(n):
        g = G.column(j)
        P = [g[s] / g[-1] for s in range(k)]

        # find t such that Psi(t) = P
        fact = gcd([psi.numerator() - p*psi.denominator() for psi, p in zip(Psi, P)])
        if fact == 1:
            new_support.append(infinity)
        else:
            new_support.append(-fact.roots()[0][0])


    if infinity in new_support:
        # reparametrize so that 
        a = Set(Fq).difference(Set(new_support)).an_element()
        very_new_support = []
        for x in new_support:
            if x == infinity:
                very_new_support.append(0)
            else:
                very_new_support.append(1 / (x-a))

    D = codes.GeneralizedReedSolomonCode(very_new_support, k+1)
    very_new_multiplier = cond(D,C).generator_matrix()[0]

    return (very_new_support, very_new_multiplier)