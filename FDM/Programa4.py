#! /usr/bin/python
'''
    Find several lowest eigenmodes V(x) and eigenenergies E
    of 1D Schrodinger equation:
    
        -1/2*hbar^2/m(d2/dx2)V(x) + U(x)V(x) = EV(x)

    For arbitrary potentials U(x).
'''
import numpy as np
import matplotlib.pyplot as plt
from math import*
from scipy.sparse import lil_matrix as lil
from scipy.sparse.linalg import eigsh

# Parameters for solving problem in the interval -L < x < L

L = 1.                       # Interval Length
N = 1000                   # No. of coordinate points

x = np.linspace(-L,L,N)     # Coordinate vector
dx = x[1]- x[0]             # Coordinate step

# POTENTIAL, choose one or make your own:

##V = 1./2.*x**2       # quadratic harmonic oscillator potential

##V = 1./2.*x**4          # quartic potential

# Finite square well of width 2w and depth given

V0 = 0.
Vizq = 1.
Vder = 1.
w = L/1.
V=np.zeros(N)
for i in range(0,N):
   if  x[i]<-w:
       V[i]= Vizq
   elif x[i]>w:
       V[i]= Vder
   else:
       V[i]= V0

##% Two finite square wells of width 2w and distance 2a apart
##V0 = 0
##Vizq = 1
##Vcent= 1
##Vder = 1
##w = 1.
##a=4
##V=np.zeros(N)
##for i in range(0,N-1):
##   if  x[i]<-a-w:
##       V[i]= Vizq
##   elif x[i]>(-a)+w and x[i]<a-w:
##        V[i]=Vcent
##   elif x[i]>a+w:
##       V[i]= Vder
##   
##   else:
##       V[i]= V0

##plt.plot(x,V)     

# Three-point finite-difference representation of Laplacian
# using sparse matrices, where you save memory by only
# storing non-zero matrix elements

Lap=lil((N,N))
Lap.setdiag(np.ones(N-1),k=-1)
Lap.setdiag(np.ones(N-1),k=1)
Lap.setdiag(-2*np.ones(N-1),k=0)
Lap=Lap/(dx**2)
Lap.tocsr()
# Total Hamiltonian

hbar = 1.; m = 1.; # constants for Hamiltonian

VMat=lil((N,N))
VMat.setdiag(V,k=0)     #Potential matrix
VMat.tocsr()
H = (-1./2.)*(hbar**2/m)*Lap + VMat 

## Find lowest nmodes eigenvectors and eigenvalues of sparse matrix

nmodes = 3
evals, evecs = eigsh(H, nmodes, which='SM',sigma=0)
print 'Eigenvalues: ', evals
plt.plot(x,evecs[:,0])
plt.show()
