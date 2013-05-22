# -*- coding: utf-8 -*-
"""
This script evaluates integration methods in python for a posible usage
inside peye qm as evaluation of local matrices in isoparametric quad 
elements.

Integration of three test functions is tested inside a square of sides 2*2
centered at the origin. 
The forementioned functions are: 
f(x, y) = x**5*y**3 
f(x, y) = sin(x)*cos(2*y)
f(x, y) = exp(-x**2 - y**2)

They are solved first analytically by using sympy, and later two  numerical 
methods are compared. The first of them is scipy's method "dblquad" which 
calls  the general purpose integration method: "scipy.integrate.quad" 
based on a technique from the Fortran library QUADPACK.
The second method has been implemented by me by using calls to two functions 
that generate legendre quadrature points and weights, one belongs to numpy 
and the other to scypy:

"numpy.polynomial.legendre.leggauss"
"scipy.special.orthogonal.p_roots"

Author: Santiago Echeverri 2013
"""


from scipy.integrate import dblquad
from scipy.special.orthogonal import p_roots
from numpy import sin, cos, exp
from numpy.polynomial.legendre import leggauss
from sympy import sin as sym_sin
from sympy import cos as sym_cos
from sympy import integrate, exp
from sympy.abc import a,b 
import time

#-------- Symbolic solution of test integrals---------------------

#sym_integral = integrate(a**5*b**3,(a,-1,1),(b,-1,1))
#sym_integral = integrate(sym_sin(a)*sym_cos(2*b),(a,-1,1),(b,-1,1))
tic=time.clock()
sym_integral = integrate(exp(-a**2 - b**2),(a,-1,1),(b,-1,1))
toc=time.clock()

print sym_integral
print 'Solution by symbolic evaluation: ', sym_integral.evalf()
print toc-tic
#------ Numérical solution of test integrals -----------------------

# Definition of an integrand

#integrand = lambda x,y: x**5*y**3 
#integrand = lambda x,y: sin(x)*cos(2*y) 
integrand = lambda x,y: exp(-x**2 - y**2)

# scypy's general purpose integration method for 2D integrals:
tic1=time.clock()

integral = dblquad(integrand,-1.0,1.0,lambda y:-1.0,lambda y: 1.0)
toc1=time.clock()

# My own non general purpose integration scheme:

# Definition of integration degree for each of the coordinates
deg_r = 3
deg_s = 3
# Generation of Gauss-Legendre points using numpy's function:
r_i, weights_r  = leggauss(deg_r)
s_j, weights_s  = leggauss(deg_s)
# Generation of Gauss-Legendre points using scypys's function:
tic2=time.clock()
r_i2, weights_r2  = p_roots(deg_r)
s_j2, weights_s2  = p_roots(deg_s)

# Cycle through each integration node and multiply weights by 
# function evaluation.
integral2 = 0
integral3 = 0
for i in range(deg_r):
    for j in range(deg_s):
        integral2 = integral2 + weights_r[i]*weights_s[j]*integrand(r_i[i],s_j[j])
        integral3 = integral3 + weights_r2[i]*weights_s2[j]*integrand(r_i2[i],s_j2[j])        

toc2=time.clock()

print 'Solution using dblquad: ', integral
print toc1-tic1
print 'Solution using leggauss: ', integral2
print 'Solution using p_roots: ', integral3
print toc2-tic2
print 'Roots taken from leggauss: ',r_i
print 'Roots taken from p_roots: ',r_i2 #, r_i-r_i2