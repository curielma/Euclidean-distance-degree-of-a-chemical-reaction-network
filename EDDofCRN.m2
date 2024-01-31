--GOAL: compute the Euclidean distance degree (EDD) of the steady-state variety of a CRN
--INPUT: a chemical reaction network G
--OUTPUT: Euclidean distance degree (integer)
--CREDIT: "The Euclidean distance degree of an algebraic variety" by
--    	      Jan Draisma, Emil Horobet, Giorgio Ottaviani, Bernd Strumfels, Rekha R. Thomas
loadPackage "ReactionNetworks"

G = reactionNetwork{"A+B --> B+C","B+C --> C+D","C+D --> A+D","A+D --> A+B"}

createRing(G,QQ)

f1 = subRandomReactionRates G
f2 = subRandomInitVals G

R = QQ[G.ConcentrationRates]
f1 = apply(f1, p -> sub(p,R))
f2 = apply(f2, p -> sub(p,R))

W = transpose gens ker transpose stoichiometricMatrix G

(Q,L,U) = LUdecomposition( sub(W,RR) )

f = f1;
for i from 0 to numRows W - 1 do(f = replace( position(flatten entries U^{i}, j -> j != 0), f2_i, f))

I = radical ideal f
u = random(ZZ^1,ZZ^(#G.ConcentrationRates)) -- EDD for generic input u

--if codimension of I is 'infinity', then Euclidean distance degree is zero
codim I

--else continue...
Ising = I + minors(codim I, jacobian(I))
M = (vars(R) - u) || (transpose jacobian(I))
J = saturate(I + minors(codim(I)+1,M), Ising)

degree J
