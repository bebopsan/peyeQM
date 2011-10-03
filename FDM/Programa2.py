#! /usr/bin/python
'''
    Calculate first and second derivative numerically
    showing how to write differential operator as a matrix
'''
import numpy as np
import matplotlib.pyplot as plt
from math import*

##Parameters for solving problem in the interval 0 < x < L

L = 2*pi                    # Interval Length
N = 500                     # No. of coordinate points

x = np.linspace(0,L,N)      # Coordinate vector
dx = x[1]- x[0]            # Coordinate step
# Two-point finite-difference representation of Derivative

D=(np.diagflat(np.ones(N-1),1)-np.diagflat(np.ones(N-1),-1))/(2*dx)

# Next modify D so that it is consistent with f(0) = f(L) = 0

D[0,0]=0;     D[1,0]=0;      D[0,1]=0              # So that f(0) = 0
D[N-1,N-1]=0; D[N-2,N-1]=0;  D[N-1,N-2]=0          # So that f(L) = 0

# Three-point finite-difference representation of Laplacian

Lap = (np.diagflat(np.ones(N-1),1)-2*np.diagflat(np.ones(N))\
       +np.diagflat(np.ones(N-1),-1))/(dx**2)

# Next modify Lap so that it is consistent with f(0) = f(L) = 0

Lap[0,0]=0;     Lap[1,0]=0;      Lap[0,1]=0              # So that f(0) = 0
Lap[N-1,N-1]=0; Lap[N-2,N-1]=0;  Lap[N-1,N-2]=0          # So that f(L) = 0

'''
     To verify that D*f corresponds to taking the derivative of f
     and Lap*f corresponds to taking a second derviative of f,
     define f = sin(x) or choose your own f
'''
f=np.sin(x)
Df=np.dot(D,f); Lapf=np.dot(Lap,f)

plt.hold(True)
plt.plot(x,f,label='sin')
plt.plot(x,Df,label='cos')
plt.plot(x,Lapf,label='-sin')
plt.axis('tight')
plt.legend()
plt.show()

