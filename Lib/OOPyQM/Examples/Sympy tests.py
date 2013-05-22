# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 09:20:41 2013

@author: santiago
"""
from sympy import *
from sympy.abc import r, s
from numpy import shape
from copy import copy
h = {}
for i in range(1,9):
    h[i] = Symbol('h'+str(i))
x = {}
for i in range(1,9):
    x[i] = Symbol('x'+str(i))
y = {}
for i in range(1,9):
    y[i] = Symbol('y'+str(i))

h[5] = Rational(1,2)*(1-r**2)*(1+s)
h[6] = Rational(1,2)*(1-s**2)*(1-r)
h[7] = Rational(1,2)*(1-r**2)*(1-s)                                
h[8] = Rational(1,2)*(1-s**2)*(1+r)
h[1] = Rational(1,4)*(1+r)*(1+s)-Rational(1,2)*h[5] \
                                -Rational(1,2)*h[8]
h[2] = Rational(1,4)*(1-r)*(1+s)-Rational(1,2)*h[5] \
                                -Rational(1,2)*h[6]
h[3] = Rational(1,4)*(1-r)*(1-s)-Rational(1,2)*h[6] \
                                -Rational(1,2)*h[7]
h[4] = Rational(1,4)*(1+r)*(1-s)-Rational(1,2)*h[7] \
                                -Rational(1,2)*h[8]
X = 0
for i in range(1, 9):
    X = X + h[i]*x[i]   
Y = 0
for i in range(1, 9):
    Y = Y + h[i]*y[i]   
J = zeros((2,2))
J[0,0] = diff(X,r)
J[0,1] = diff(X,s)
J[1,0] = diff(Y,r)
J[1,1] = diff(Y,s)
jac_det = J.det()
J_inv = zeros((2,2))
J_inv[0,0] = J[1,1]
J_inv[1,1] = J[0,0]
J_inv[0,1] = -J[1,0]
J_inv[1,0] = -J[0,1]
#J_inv = J_inv*(1/jac_det)

# Test to see if weight functions work
nodes=[{r:1,s:1},{r:-1,s:1},{r:-1,s:-1},{r:1,s:-1},{r:0,s:1},{r:-1,s:0},\
          {r:0,s:-1},{r:1,s:0}]
for j in range(1,9):
    print 'function h_',j
    for i in range(8):
        pprint(h[j].subs(nodes[i]))     

#--------------------------------------
pprint(J.subs({x[1]:1, y[1]:1, x[2]:-1, y[2]:1, x[3]:-1, y[3]:-1, \
                 x[4]:1,y[4]:-1, x[5]:0,y[5]:1, x[6]:-1,y[6]:0,\
                 x[7]:0,y[7]:-1, x[8]:1,y[8]:0}))
pprint(J[1,1])
pprint(J[1, 1].subs({x[1]:1, y[1]:1, x[2]:-1, y[2]:1, x[3]:-1, y[3]:-1, \
                 x[4]:1,y[4]:-1, x[5]:0,y[5]:1, x[6]:-1,y[6]:0,\
                 x[7]:0,y[7]:-1, x[8]:1,y[8]:0}))
print 'something i dont get'
pprint(diff(h[1],r))
print shape(J)

def numeric_J(r,s,x,y):

    J_mat = zeros((2,2))
    dhdr = [1.0/4.0*(s**2 + s +2*r*(s+ 1)), \
            1.0/4.0*(-s**2 - s +2*r*(s+ 1)), \
            1.0/4.0*(-s**2 + s +2*r*(-s+ 1)), \
            1.0/4.0*(s**2 - s +2*r*(-s+ 1)), \
            -r*(s+ 1), \
            -1.0/2.0*(-s**2+ 1), \
            -r*(-s+ 1),\
            1.0/2.0*(-s**2+ 1)]
    
    dhds = [1.0/4.0*(r**2 + r +2*s*(r+ 1)), \
            1.0/4.0*(r**2 - r +2*s*(-r+ 1)), \
            1.0/4.0*(-r**2 + r +2*s*(-r+ 1)), \
            1.0/4.0*(-r**2 - r +2*s*(r+ 1)), \
            1.0/2.0*(-r**2 + 1), \
            -s*(-r+ 1), \
            -1.0/2.0*(-r**2+ 1), \
            -s*(r+ 1)]

    for i in range(8):
        J_mat[0,0] = J_mat[0,0] + x[i]*dhdr[i]
        J_mat[1,0] = J_mat[1,0] + y[i]*dhdr[i]
        J_mat[0,1] = J_mat[0,1] + x[i]*dhds[i]
        J_mat[1,1] = J_mat[1,1] + y[i]*dhds[i]

    return J_mat

x = [1.0, -1.0,-1.0,  1.0, 0.0,-1.0, 0.0, 1.0]
y = [1.0,  1.0,-1.0, -1.0, 1.0, 0.0,-1.0, 0.0]

Jacobian = numeric_J(r,s,x,y)
print  Jacobian

from sympy.abc import x,y
Q = integrate(integrate(x**5*y**3,(x,-1,1)),(y,-1,1))
print Q
