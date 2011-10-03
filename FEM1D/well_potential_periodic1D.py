#! /usr/bin/python
# -*- coding: utf-8 -*-
# Falta una descrición del script
from numpy import zeros
from scipy import linalg
import matplotlib.pyplot as plt
import math
import cmath
import numpy
#Vector de función de onda Psi debe estar evaluado sobre el ancho del pozo "a"
#y con divisiones dadas por el ancho dividido el número de divisiones:

#a logitud del pozo

a=math.pi

# N: número de divisiones
N=20


# L: tamaño de cada elemento
L=a/(N-1)

print '## ELECTRON IN A PERIODIC POTENTIAL IN 1D ##\n', 'Cell Length: ',a

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
V0 = 1.0
Vizq = 10
Vder = 3
ancho = a/2.0
for i in range(0,N-1):
    if  X[i]<(a/2.-ancho/2.):
        V[i]= Vizq
    elif X[i]>(a/2.+ancho/2.) :
        V[i]= Vder
    else:
        V[i]= V0


##plt.figure(1)
##plt.plot( (X[0:N-1]+X[1:N])/2. ,V)
##plt.figure(1).suptitle('Potencial')

print 'V[0]:\n',V[0],'\nL:\n','Element size:\n',L

# Inicializa la matriz equivalente de rigidez 
K=zeros((N,N),dtype=numpy.cfloat)
# Inicializa la matriz equivalente de masa
M=zeros((N,N),dtype=numpy.cfloat)

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

print 'K shape is:\n',K.shape

    
# Bloch-Periodicity imposition

xi=X[0]    # initial x
xf=X[N-1]  # final x

nVals =  8  # number of eigenvales to compute

nk = 101  # number of k to sweep
kmax = 4.*math.pi/a
kmin = -0
k_range = numpy.linspace(kmin, kmax, num=nk)
omega = numpy.zeros( (len(k_range),nVals) )
E = numpy.zeros( (len(k_range),nVals) )
                  
print 'Number of eigenvales to compute: ', nVals,'\nNumber of wave numbers to sweep: ', nk, ' in ',  [k_range[0],k_range[nk-1]]

ll = 0

Kaux = K.copy()
Maux = M.copy()

for k in k_range:
    fi=cmath.exp(1.0j*k*xi)
    ff=cmath.exp(1.0j*k*xf)
    K = Kaux.copy()
    M = Maux.copy()


    for i in range(0,N):
        K[0,i]=K[0,i]*fi.conjugate()
        K[i,0]=K[i,0]*fi
        K[N-1,i]=K[N-1,i]*ff.conjugate()
        K[i,N-1]=K[i,N-1]*ff
        
        M[0,i]=M[0,i]*fi.conjugate()
        M[i,0]=M[i,0]*fi
        M[N-1,i]=M[N-1,i]*ff.conjugate()
        M[i,N-1]=M[i,N-1]*ff
        
    K[N-1,:] = K[0,:] + K[N-1,:]
    K[:,N-1] = K[:,0] + K[:,N-1]

    M[N-1,:] = M[0,:] + M[N-1,:]
    M[:,N-1] = M[:,0] + M[:,N-1]


    Kd=K[1:N,1:N]
    Md=M[1:N,1:N]

    vals = linalg.eigvalsh(Kd,Md,eigvals=(0,nVals-1) )
    
    for i in range(0,nVals):
        omega[ll,i] = math.sqrt( abs(vals[i]) )

    for i in range(0,nVals):
        E[ll,i] = vals[i]
        
    ll = ll + 1



plt.figure(2)
plt.hold(True)
legend=[]
for i in range(0,nVals):
    plt.plot(k_range,E[:,i])
    legend.append('n = '+str(i+1))

plt.plot(k_range,(k_range)**2,'--k')
legend.append('Free Electron')
plt.title('Dispersion relation')
plt.legend(legend,loc=2)
plt.xlabel('Adimensional wave number - $a\kappa/\pi$')
plt.ylabel('Adimensional energy - $2ma^2E/\hslash^2$')
plt.grid()
plt.show()

numpy.savetxt('Energy.txt', omega,fmt='%1.4e')

