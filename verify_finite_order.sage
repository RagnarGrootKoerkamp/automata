# Check the minimal polynomial for sigma.
# This is used to check whether automata/polynomials that pass the 
# Pass in the output from main.cpp, e.g. results/size_to_5_break_sequence_check.txt

import fileinput

class Automaton:
    pass

def read_automata():
    automatons = []
    for line in fileinput.input():
        if line != 'Candidate automaton\n': continue

        val = Automaton()

        val.number = int(input().split()[1])
        val.order = int(input().split()[1])
        line = input()
        prefix = 'Vertex labels: '
        assert line.startswith(prefix)
        line = line[len(prefix):]
        val.labels = list(map(int, line.split()))
        assert input() == 'Edges:'
        val.n = len(val.labels)
        val.edges = [list(map(int, input().split())) for _ in range(val.n)]
        assert input() == ''

        def read_sigma(exp):
            sigma = input()
            sigma_prefix = f'sigma^{exp}: '
            assert sigma.startswith(sigma_prefix)
            sigma = sigma[len(sigma_prefix):].replace('x', 't')
            idx = sigma.find('+ O(')
            assert idx >= 0
            return (sigma[:idx], int(sigma[idx + 6:sigma.find(')')]))
        val.sigma, val.sigma_prec = read_sigma(1)
        val.sigma2, val.sigma2_prec = read_sigma(2)
        val.sigma4, val.sigma4_prec = read_sigma(4)

        automatons.append(val)

    return automatons

def verify_finite_order(a):
    n = a.n

    # Handle order 1 separately.
    if a.order == 1:
        print(f'Automaton #{a.number} passes the checks!')
        return True
    assert a.order in [2,4]

    # 1. Determine the polynomial in t for the current automaton.

    # Make n variables: one A_i for each node. s=A0 is the start node.
    names = ['t', 's'] + [f'A{i}' for i in range(1, n)]
    R = PolynomialRing(GF(2), n+1, names)
    variables = R.gens()
    t = variables[0]
    A = variables[1:]

    # Convert the string representation to an actual polynomial in R.
    # New ring to cover the precision of sigma.
    a.sigma = R(a.sigma)
    a.sigma2 = R(a.sigma2)
    a.sigma4 = R(a.sigma4)
    deg = a.sigma.degree(t)
    R2 = R.quotient(t^(deg+1))
    a.sigma = R2(a.sigma)

    # Make n equations, one for each node.
    equations = [A[i] + A[a.edges[i][0]]^2 + t*A[a.edges[i][1]]^2 for i in range(n)]
    I = R.ideal(equations)

    # Reduce the system of equations to a single equation in A_0 = sigma.
    J = I.elimination_ideal(A[1:])
    assert len(J.gens()) == 1
    gen = J.gens()[0]
    # Factorize the polynomial.
    factors = list(gen.factor())
    #print(factors)

    # If the order of sigma is finite, it can only be a root of factors with deg(t) = deg(s).
    found = False
    for F, exp in factors:
        deg_t = F.degree(t)
        deg_s = F.degree(A[0])
        if deg_t != deg_s: continue

        # Check that sigma is a root of this factor.
        Fs = F.subs(s = a.sigma)
        #print('Substitute s in F:', F, Fs, sep='\n')
        if Fs != 0: continue

        # 2. We found the minimal polynomial F(t, sigma) = 0.
        #    We also have F(sigma, sigma(sigma)) = 0 after substituting t=sigma(t).
        #    From this we extract G(t, sigma^2) = 0.

        # Check that G(t,t) = 0 and sigma^order = t + O(t^valuation(Res_t(G(t, sigma), (dG/dt)(t, sigma))))
        if a.order==2:
            # ss = sigma^2
            R3 = PolynomialRing(GF(2), n+2, names+['ss'])
            s = R3.gens()[1]
            s2 = R3.gens()[-1]
            ss = s2
            s4 = None
            # Equation for sigma^2.
            F2 = F.subs(s = s2, t=A[0])
            K = R3.ideal(F, F2)
            # Eliminate sigma to keep t and sigma^2.
            L = K.elimination_ideal([s] + list(A[1:]))
        else:
            # ss = sigma^4
            R3 = PolynomialRing(GF(2), n+4, names+['s2', 's3', 'ss'])
            s = R3.gens()[1]
            s2 = R3.gens()[-3]
            s3 = R3.gens()[-2]
            s4 = R3.gens()[-1]
            ss = s4
            # Equation for sigma^2.
            F2 = F.subs(s = s2, t=A[0])
            # Equation for sigma^3.
            F3 = F.subs(s = s3, t=s2)
            # Equation for sigma^4.
            F4 = F.subs(s = s4, t=s3)
            K = R3.ideal(F, F2, F3, F4)
            # Eliminate sigma and sigma^2 to keep t and sigma^4.
            L = K.elimination_ideal([s, s2, s3] + list(A[1:]))
        assert len(L.gens()) == 1
        G = L.gens()[0]

        # 3. Check G(t,t) = 0
        Gtt = G.subs(ss= t)
        if Gtt != 0:
            # G(t, t) must be 0.
            #print('FAILED G(t,t)=0 check!')
            continue

        #print('order', a.order)
        #print('F', F)

        #print('G',G)
        #print('Gens', R3.gens())

        # 4. Check sigma^2 = t + O(t^(m+1)) where m = Res_X(G(t, X), (dG)/(dX)(t, X)).
        dG = derivative(G, ss)
        #print('dG', derivative(G, ss))
        res = G.resultant(dG, ss)
        #print('res', res)
        # t-Valuation of res.
        m = min(map(lambda p: p.degree(), res.monomials()))
        #print('m', m)
        assert res.is_univariate()
        if a.order==2:
            assert a.sigma2_prec > m 
            #print(a.sigma2)
            assert a.sigma2 == t
        else:
            assert a.sigma4_prec > m 
            #print(a.sigma4)
            assert a.sigma4 == t
        print(f'Automaton #{a.number} passes the checks!')
        return True
    return False

for automaton in read_automata():
    assert verify_finite_order(automaton)

print('All automata passed the checks.')
