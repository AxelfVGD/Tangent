load('codecrypto.sage')


# Griffiths, P. and Harris, J. in Principles of Algebraic Geometry:
# "Through any n+3 points in general position in P^n there passes a unique rational normal curve."
# Finding the parametrization enables one to recover a support and a multiplier of a given GRS code

# Given G a generator matrix of C = GRS_k(x,y), recovers x' and y' such that C = GRS_k(x',y')
def GH78(G):
    Fq = G.base_ring()
    r = G.nrows()-1
    P1.<u,v> = ProjectiveSpace(Fq,1)
    P = ProjectiveSpace(Fq,r)
    R.<t> = Fq['t']

    Points = [P(*g) for g in G.columns()]

    Psi = rational_normal_curve(Points,BasePoints=Points[r:r+3])
    Psi *= lcm(psi.denominator() for psi in Psi)

    # find the new support, possibly with infinity
    new_support = []
    for g in Points:
        if P(*Psi(1,0)) == g:
            new_support.append(P1(1,0))
        else:
            i = r
            while g[i] == 0:
                i -= 1
            f = gcd(psi(t,1).numerator() - gi/g[i] * psi(t,1).denominator() for gi,psi in zip(g,Psi/Psi[i]))
            new_support.append(P1(f.roots()[0][0],1))

    if P1(1,0) in new_support:
        a = Set(P1.rational_points()).difference(new_support).an_element()[0]
        new_support = [P1(p[1],p[0]-a*p[1]) for p in new_support]
        Phi = Psi(v+a*u,u)
        Phi *= lcm(phi.denominator() for phi in Phi)

    # now find the multiplier
    new_multiplier = []
    for g,x in zip(G.columns(),new_support):
        i = r
        while g[i] == 0:
            i -= 1
        new_multiplier.append(g[i] / Phi(*x)[i])

    return ([p[0] for p in new_support],new_multiplier)

# A specialized function when y = 1 / Gamma(x) : this function works in P^r instead of P^(r-1) and directly recovers a Goppa polynomial
def GH78Goppa(G):
    Fq = G.base_ring()
    r = G.nrows()
    P1.<u,v> = ProjectiveSpace(Fq,1)
    P = ProjectiveSpace(Fq,r)
    R.<t> = Fq['t']

    Points = [P(*(list(g)+[1])) for g in G.columns()]

    Psi = rational_normal_curve(Points,BasePoints=[Points[r],Points[r+1]]+[P([0 for _ in range(r)]+[1])])
    Psi *= lcm(psi.denominator() for psi in Psi)

    # find the new support. Since the third base points was (0:...:0:1), the point at infinity (1:0) certainly is not in the new support
    new_support = []
    for g in Points:
        i = r
        while g[i] == 0:
            i -= 1
        f = gcd(psi(t,1).numerator() - gi/g[i] * psi(t,1).denominator() for gi,psi in zip(g,Psi/Psi[i]))
        new_support.append(P1(f.roots()[0][0],1))

    # no need to recover the multiplier since we already have the Goppa polynomial

    return ([p[0] for p in new_support],Psi[-1](t,1).numerator())