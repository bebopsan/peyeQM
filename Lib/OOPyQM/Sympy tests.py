# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 09:20:41 2013

@author: santiago
"""
from sympy import *
from sympy.abc import r, s

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
pprint(X)
J = zeros((2,2))
J[0,0] = diff(X,r)
J[0,1] = diff(X,s)
J[1,0] = diff(Y,r)
J[1,1] = diff(Y,s)

jac_det = J.det()

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
