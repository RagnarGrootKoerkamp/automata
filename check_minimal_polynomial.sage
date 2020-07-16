# Check the minimal polynomial for sigma.
# This is used to check whether automata/polynomials that pass the 
# Pass in the output from main.cpp, e.g. results/size_to_5_break_sequence_check.txt

import fileinput

def read_automata():
    for line in fileinput.input():
        #print(line)
        prefix = 'Vertex labels: '
        if line.startswith(prefix):
            line = line[len(prefix):]
            labels = list(map(int, line.split()))
            assert input() == 'Edges:'
            n = len(labels)
            edges = [list(map(int, input().split())) for _ in range(n)]
            assert input() == ''
            sigma = input()
            sigma_prefix = 'sigma^1: '
            assert sigma.startswith(sigma_prefix)
            sigma = sigma[len(sigma_prefix):].replace('x', 't')

            yield (n, labels, edges, sigma)

def check_automaton(a):
    n, l, edges, s = automaton
    #print(n,l,edges,s)

    # Make n variables A_i for each node. sigma=A_0 is the start node.
    names = ['t'] + [f'A{i}' for i in range(n)]
    R = PolynomialRing(GF(2), n+1, names)
    variables = R.gens()
    t = variables[0]
    A = variables[1:]
    #print(A)
    #print(t)

    # Convert the string representation to an actual polynomial in R.
    # New ring to cover the precision of sigma.
    s = R(s)
    #print(s)
    deg = s.degree(t)
    R2 = R.quotient(t^(deg+1))
    s = R2(s)

    # Make n equations, one for each node.
    equations = [A[i] + A[edges[i][0]]^2 + t*A[edges[i][1]]^2 for i in range(n)]
    #print(equations)
    I = R.ideal(equations)
    J = I.elimination_ideal(A[1:])
    #print('IDEAL: ', J)
    assert len(J.gens()) == 1
    gen = J.gens()[0]
    #print('GENERATOR: ', gen)
    try:
        factors = list(gen.factor())
    except:
        print('ERROR WHILE FACTORIZING: XXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
        print(gen)
        return
    #print(factors)

    # If the order of sigma is finite, it can only be a root of factors with deg(t) = deg(s).
    found = False
    for f, exp in factors:
        deg_t = f.degree(t)
        deg_s = f.degree(A[0])
        if deg_t != deg_s: continue

        # Substitute sigma for A0 and see if it is a factor modulo the degree of sigma.
        if f.subs(A0 = s) != 0: continue
        print()
        print('FOUND A POTENTIAL FINITE ORDER CANDIDATE')
        print('FACTOR', f)
        print(f.subs(A0=s))
        found = True
        print(n, l, edges)
        print(s)
        print(factors)
        print()
        break
    if not found:
        print('THIS AUTOMATON HAS INFINITE ORDER')


for automaton in read_automata():
    check_automaton(automaton)
