#! /usr/bin/python
# -*- coding: utf-8 -*-
import numpy
import math
from numpy import zeros,savetxt
from scipy import linalg
import matplotlib.pyplot as plt

# Tomado de "Computing quantum eigenvalues made easy" de
# Korsch y Gluck [2002]

N = 50
s = 1.0

n = numpy.linspace(1,N,N)
m = zeros( (N) )

for i in n:
    m[i-1] = math.sqrt(i)


x = zeros( (N,N) ,dtype = numpy.cfloat)
p = zeros( (N,N) ,dtype = numpy.cfloat)

for i in n[0:N-1]:
        x[i-1,i]= m[i-1]
        p[i-1,i]= -m[i-1]
        x[i,i-1]= m[i-1]
        p[i,i-1]= m[i-1]

x = s/math.sqrt(2) * x
p = 1.j/s/math.sqrt(2) * p

H = (numpy.dot(p,p)/2 + numpy.dot(x,x)/2 ).real

vals = linalg.eigvalsh(H,eigvals = (0,8))

print 'The first 8 eigenvalues are:\n', vals
        
