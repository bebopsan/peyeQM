#! /usr/bin/python
# -*- coding: utf-8 -*-
# Falta una descrición del script
from numpy import zeros
from scipy import linalg
import matplotlib.pyplot as plt
from math import sin,pi,exp, cosh
#Vector de función de onda Psi debe estar evaluado sobre el ancho del pozo "a"
#y con divisiones dadas por el ancho dividido el número de divisiones:

#a logitud del pozo

a=4.

# N: número de divisiones
N=500


# L: tamaño de cada elemento
L=a/(N-1)

# Inicializo el vector de coordenadas.
X=zeros((N))

# Asigna los valores correspondientes de coordenadas a cada nodo.
for i in range(0,N):
    X[i]=i*L

# Inicializa el vector de potencial

V=zeros(N-1)

##### Potencial #####

### Constante
##V0 = 0.
##for i in range(0,N-1):
##    V[i]= V0

### Pozo finito
Vizq = -5.
Vder = -10.
ancho = 1
for i in range(0,N-1):
    if  X[i]<(a/2.-ancho/2.):
        V[i]= Vizq
    elif X[i]>(a/2.+ancho/2.) :
        V[i]= Vder
    else:
        V[i]= 0.

### x
##alpha = 1.
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = alpha*(x-a/2)
##    V[i]=v

### 1/sin(x) -1
##for i in range(0,N-1):
##    V[i]= 1/sin( abs(X[i]+X[i+1])/2*pi/X[N-1] ) -1
##    
### Posch-Teller simétrico
##V0 = 1.0
##b = a/pi
##lam=2.
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2
##    v = 2.*V0*lam*(lam-1)/( sin(x/b) )**2
##    V[i]=v

### Potencial de Landau-Lifhshitz
##V0 = -20.0
##alpha = 1.0
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = V0/( cosh( alpha*(x-a/2) ) )**2
##    V[i]=v

### Barrera de potencial
##V0 = 5.0
##ancho = 2.0
##for i in range(0,N-1):
##    if  X[i]<(a/2.-ancho/2.):
##        V[i]= 0.
##    elif X[i]>(a/2.+ancho/2.) :
##        V[i]= 0.
##    else:
##        V[i]= V0 
    
### Oscilador armónico
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = (x-a/2)**2
##    V[i]=v

### Oscilador armónico - de orden 4
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = (x-a/2)**4
##    V[i]=v

### Doble pozo de orden 4
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = -2.*(x-a/2)**2 + (x-a/2)**4
##    V[i]=v

# Potencial de Morse
##xe = a/2.
##D =1.
##a =1.
##for i in range(0,N-1):
##    x = (X[i]+X[i+1])/2.
##    v = D*(1. - exp(-a*(x-xe)) )**2
##    V[i]=v
##
plt.figure(1)
plt.plot( (X[0:N-1]+X[1:N])/2. ,V)
plt.figure(1).suptitle('Potencial')

print 'V[0]:\n',V[0],'\nL:\n','Element size:\n',L

# Inicializa la matriz equivalente de rigidez 
K=zeros((N,N))
# Inicializa la matriz equivalente de masa
M=zeros((N,N))

K[0,0]= 1/L + L*V[0]/3
K[N-1,N-1]= 1/L + L*V[N-2]/3
for i in range(0,N-1):
    if i!=0:
        K[i,i]=2/L + L*(V[i-1]+V[i])/3 
    K[i,i+1]= -1/L + L*V[i-1]/6
    K[i+1,i]= -1/L + L*V[i-1]/6
    
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

nVals=4

V,Dd=linalg.eigh(Kd,Md,eigvals=(0,nVals-1))

D=zeros((N,nVals))
D[1:N-1,:]=Dd

print 'D shape is:\n', D.shape,'\n','The Eigenvalues are:\n',V

plt.hold(True)
plt.figure(2)
legend=[]
for i in range(0,nVals):
    plt.plot(X,D[:,i])
    legend.append('n = '+str(i+1))

plt.legend(legend)
plt.show()
