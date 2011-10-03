#! /usr/bin/python
# -*- coding: utf-8 -*-
from numpy import zeros,savetxt
##"import numpy.savetxt as save"
from scipy import linalg
import matplotlib.pyplot as plt
#Vector de función de onda Psi debe estar evaluado sobre el ancho del pozo a
#y con divisiones dadas por el ancho dividido el número de divisiones:

#a logitud del pozo

a=2.

#N: número de divisiones
N=1000


#L: tamaño de cada elemento
L=a/(N-1)
print L

#Inicializo el vector.
X=zeros((N))
for i in range(0,N):
    X[i]=i*L


#Inicializa la matriz equivalente de rigidez 
K=zeros((N,N))
#Inicializa la matriz equivalente de masa
M=zeros((N,N))

K[0,0]=1/L
K[N-1,N-1]=1/L  

for i in range(0,N-1):
    if i!=0:
        K[i,i]=2/L
    K[i,i+1]=-1/L
    K[i+1,i]=-1/L

M[0,0]=L/3
M[N-1,N-1]=L/3
for i in range(0,N-1):
    if i!=0:
        M[i,i]=2*L/3
    M[i,i+1]=L/6    
    M[i+1,i]=L/6  

Kd=K[1:N-1,1:N-1]
Md=M[1:N-1,1:N-1]
print 'K shape is:\n',K.shape
savetxt('matrixK.txt',K,fmt='%0.5f')
savetxt('matrixM.txt',M,fmt='%0.5f')

nVals=4

V,Dd=linalg.eigh(Kd,Md,eigvals=(0,nVals-1))

D=zeros((N,nVals))
D[1:N-1,:]=Dd

print 'D shape is:\n', D.shape,'\n','The Eigenvalues are:\n',V

plt.hold(True)
legend=[]
for i in range(0,nVals):
    plt.plot(X,D[:,i])
    legend.append('n = '+str(i+1))

plt.legend(legend)
plt.show()
