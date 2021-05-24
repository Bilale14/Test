#=
	Anne-CAtherine Letournel
	09/03/2021
	TP1 Calcul scientifique: méthode de Newton
=#

# fonction Newton, exercice 1

# calcul d'une solution particulière par la méthode de Newton
function Newton(x0, f, df, tol, itmax)
	N = 1
	x1 = x0
	x2 = x1 - f(x1)/df(x1)
	while ( abs(x2-x1) > tol && N <= itmax )
		x1 = x2
		x2 = x1 - f(x1)/df(x1)
		N = N+1
	end
	return(x2,N)
end

# script de test
# initialisation des variables d'entrée
x0 = -2.5
f = x -> x^3 -4*x + 1
df = x -> 3*x^2 - 4
tol = 1e-6
itmax = 100

# appel de la fonction et affichage des résultats
y, N = Newton(x0, f, df, tol, itmax)
println("f(",y,") = ", f(y), " trouvé en ", N, " itérations.")

x0 = 0.0

# appel de la fonction et affichage des résultats
y, N = Newton(x0, f, df, tol, itmax)
println("f(",y,") = ", f(y), " trouvé en ", N, " itérations.")

x0 = 2.0

# appel de la fonction et affichage des résultats
y, N = Newton(x0, f, df, tol, itmax)
println("f(",y,") = ", f(y), " trouvé en ", N, " itérations.")



# cas particulier d'un polynome de degré 3
function NewtonPoly(x0 ,a, b, c, d, tol, itmax)
	N = 1
	x1 = x0
	x2 = x1 - (a*x1^3 + b*x1^2 + c*x1 + d) / (3*a*x1^2 + 2*b*x1 + c)
	while ( abs(x2-x1) > tol && N <= itmax )
		x1 = x2
		x2 = x1 - (a*x1^3 + b*x1^2 + c*x1 + d) / (3*a*x1^2 + 2*b*x1 + c)
		N = N+1
	end
	return(x2, N)
end

# script de test
# initialisation des variables d'entrée
x0 = -2.5
a = 1
b = 0
c = -4
d = 1
tol = 1e-6
itmax = 100

# appel de la fonction et affichage des résultats
x1, N = Newton(x0, f, df, tol, itmax)
println("f(",x1,") = ", f(x1), " trouvé en ", N, " itérations.")



