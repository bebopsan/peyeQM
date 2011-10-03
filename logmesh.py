#! /usr/bin/python

from numpy import *
import matplotlib.pyplot as plt


a=1.0
b=3.0
N=21
base = 5.0
option ="R"


if(option=="B"):
    if(mod(N,2)==0):
        x = zeros((N+1))
    else:
        x = zeros((N))
    xaux = logspace(0, 1, N/2+1, endpoint=True,base=base)
    xaux = (xaux-1.0)/(base-1.0)
    for i in range(0,N/2+1):
        x[i] = -xaux[N/2-i]
        x[i+N/2] = xaux[i]
    x = x*(b-a)/2.0 + (b+a)/2.0
elif(option=="L"):
    x = zeros((N))
    xaux = logspace(0, 1, N, endpoint=True,base=base)
    x = (xaux-1.0)*(b-a)/(base-1.0)+a  
elif(option=="R"):
    x = zeros((N))
    xaux = logspace(0, 1, N, endpoint=True,base=base)
    for i in range(0,N):
        x[i] = -xaux[-i-1]+base
    x = x*(b-a)/(base-1.0)+a
else:
    print "Not known value of parameter 'option' ",option

print x
y = 0.0*x

plt.plot(x, y, 'o')
plt.show()
