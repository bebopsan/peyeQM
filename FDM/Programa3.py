#! /usr/bin/python
'''
    Matrix representation of differential operators,
    solving for Eigenvectors & Eigenvalues of Infinite Square Well
'''
import numpy as np
import matplotlib.pyplot as plt
from math import*
from scipy import linalg

L = 2                    # Interval Length
N = 500                    # No. of coordinate points

x = np.linspace(0,L,N)      # Coordinate vector
dx = x[1]- x[0]            # Coordinate step

# Three-point finite-difference representation of Laplacian

Lap = (np.diagflat(np.ones(N-1),1)-2*np.diagflat(np.ones(N))\
       +np.diagflat(np.ones(N-1),-1))/(dx**2)

# Total Hamiltonian where hbar=1 and m=1
hbar = 1; m = 1
H = -(1./2)*(hbar**(2./m))*Lap
# Solve for eigenvector matrix V and eigenvalue matrix E of H
E,V = linalg.eigh(H,eigvals=(0,3))
# Plot lowest 3 eigenfunctions

plt.hold(True)
legend=[]
print E
for i in range(0,3):
    plt.plot(x,V[:,i])
    legend.append('n = '+str(i+1))
plt.axis('tight')
plt.legend(legend)
plt.show()

