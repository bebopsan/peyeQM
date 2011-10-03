## module PrePro
# -*- coding: latin-1 -*-
''' This module contains the necessary functions for the preprocessing
    phase of finite element procedure for the solution of differentia
    equations.
'''

__all__=['Mesh1D','Potential1D','well','finwell','xwell','poschl','oscil']
__author__="Santiago Echeverri Chacón"
import math
from numpy import zeros,array
def Mesh1D(Type,a=math.pi,N=100):
    '''
        In the mesh function the user can decide wether to use an
        homogeneus discretization of the line of lenght a in N elements,
        or give the name of a .msh file that contains the discretized line.

        Parameters:
        ----------
        Type:  Is a string variable that will tell the kind
                of Mesh to be made.

		"simple":   Will execute the discretization using
                            N divisions over the line [xmin,xmax]
		"left":     Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values on
                            the left.
		"right":    Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values
                            on the right side.
		"center":   Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values on
                            the middle.
        
        a:  Lenght of the line
        N:  Number of divisions

        The default values  for a and N are a=4 and N=100.
        The user can define new values by assigning them as arguments of the
        function e.g. a=10,N=100000.

        Returns:
        --------
        X:  numpy array like vector of the discretized domain with lenght N

        Raises:
        -------
        Error if the user gives a wrong Type argument
        
        Last modification: 2/Oct/2011
    '''

    if 'simple' in Type:

        # L: size of each element
        L=a/(N-1)
        # Initilization of the coordinate vector
        X=zeros((N))
        # Coordinates values for each node in the domain
        for i in range(0,N):
            X[i]=i*L
        
    elif 'left' in Type:
        from numpy import sqrt,logspace
        import matplotlib.pyplot as plt
        X=logspace(0.001,math.pi,num=N)
        X=X*math.pi/10**math.pi
        plt.plot(X,zeros(X.size),'+')
        plt.show()
        return X
    elif 'right' in Type:
        from numpy import sqrt,logspace
        import matplotlib.pyplot as plt
        
        X=logspace(0.001,math.pi,num=N)
        X=X*math.pi/10**math.pi
        disp=X[X.size-1]
        X=X[::-1]
        X=array(X)
        X=X*-1+disp
        plt.plot(X,zeros(X.size),'+')
        plt.show()
        return X
        
    elif 'center' in Type:
        from numpy import sqrt,concatenate,copy,hstack,logspace
        import matplotlib.pyplot as plt
        X=logspace(0.001,math.pi,num=N)
        X=X*math.pi/10**math.pi
        disp=X[X.size-1]
        Y=copy(X[::-1])
        Y=X*-1+disp
        X=X+disp
        X=hstack((Y,X))
        plt.plot(X,2*zeros(X.size),'+')
        plt.show()
        return X
    else:
        print 'You entered a wrong Type parameter. Please try again, or read the\
              documentation'
        
    
def Potential1D(Type,X,V0=1,Vleft=0.0,Vright=0.0,width=1,lam=2.,D=1):
    '''
        This funcion is meant to evaluate different types of one dimensional
        potentials for the solution of Schrödinger equation.

        Each potential represents a particular aproach to one dimensional
        quantum mechanical phenomena.

        The potential to be evaluated is selected by the user by giving
        a string variable as the first argument of the function. The list
        of available potentials is given below

        * "well"
        * "finwell"
        * "xwell"
        * "poschl"
        * "oscil"
        * "morse"

        The next argument is the position vector X output from the "Mesh1D"
        function of this module. this vector contains the coordinate for each
        node of the domain to be evaluated.

        All the other arguments are specific for each of the available
        potentials. They are given with predefined values in order to ease
        the procedure to the user, however this values can be modified
        by assigning new ones as arguments of the function.
        For example:
                        Potential1D("well",X,V0=5)
    '''
   
    
    if 'well' in Type:
       
        V=well(X,V0)
        return V
    elif 'finwell' in Type:
        V=finwell(X,V0,Vleft,Vright,widht)
        return V
    elif 'xwell' in Type:
        V=xwell(X,lam)
        return V
    elif "poschl" in Type:
        V=poschl(X,V0,lam)
        return V
    elif "oscil" in Type:
        V=oscil(X)
        return V
    elif "morse" in Type:
        V=morse(X,D,lam)
        return V
    else:
         print 'You entered a wrong Type parameter. Please try again, or read the\
              documentation'
       

def well(X,V0):
    """
well":          Defines a well with constant potential over the 1D domain
                 given by vector X

                 inf ____  _____inf
                        |  |
                        |__|___V0
                        
                 Parameters:
                 -----------
                 V0:  Real value of the potential distributed over the well
                 X:   numpy 1D array output of the Mesh1D function.

                 Returns:
                 --------
                 V: numpy array like vector with the values for the potential
                    over each node on X

                 Note: more info in square wells can be found at
                 http://en.wikipedia.org/wiki/Particle_in_a_box
    """
    N=size(X)
    V=zeros(N-1)
    V0 = 0.
    for i in range(0,N-1):
        V[i]= V0
    return V

def finwell(X,V0,Vleft,Vright,widht):
    """
"finwell":	  Works like "well" but  with control over the height of
                  the left and right walls of the well:

                             _____Vright
                  Vleft___   |
                         |   |     
                         |___|__V0
                         width
                    
                 Parameters:
                 -----------
                 V0:      Value of the potential over the bottom of the well
                 Vleft:   Value of the left wall
                 Vright:  Value of the right wall
                 width:   lenght between the left and right walls of the well 
                 X:   numpy 1D array output of the Mesh1D function.

                 Returns:
                 --------
                 V: numpy array like vector with the values for the potential
                    over each node on X
    """
    N=size(X)
    V=zeros(N-1)
    a=X[N-1]
    Vleft = -5.
    V0=0.
    Vright = -10.
    widht = 1
    for i in range(0,N-1):
    	if  X[i]<(a/2.-widht/2.):
		V[i]= Vleft
	elif X[i]>(a/2.+widht/2.) :
		V[i]= Vright
	else:
		V[i]= V0
    return V

def xwell(X,lam):
    """
"xwell":        Potential well defined as a straight line with slope given
                   by "lam"
                     inf___     ___
                           |   /          
                           |  /           
                           | /
                       0___|/
                       
                 Parameters:
                 -----------
                 lam:  Slope of the straight line
                 X:   numpy 1D array output of the Mesh1D function.

                 Returns:
                 --------
                 V: numpy array like vector with the values for the potential
                    over each node on X
    """
    N=size(X)
    V=zeros(N-1)
    a=X[N-1]
    lam = 1.
    for i in range(0,N-1):
        x = (X[i]+X[i+1])/2.
        v = lam*(x-a/2)
        V[i]=v
    return V

def poschl(X,V0,lam):
    """
"poschl":       Potential Well defined by G.Pöschl and Edward Teller
                  described in:
                    http://en.wikipedia.org/wiki/Pöschl-Teller_potential

                 Parameters:
                 -----------
                 V0:   Upper limit of the potential
                 lam:  lambda parameter
                 X:   numpy 1D array output of the Mesh1D function.

                 Returns:
                 --------
                 V: numpy array like vector with the values for the potential
                    over each node on X
    """

    N=size(X)
    V=zeros(N-1)
    a=X[N-1]
    V0 = 1.0
    b = a/pi
    lam=2.
    for i in range(0,N-1):
            x = (X[i]+X[i+1])/2
            v = 2.*V0*lam*(lam-1)/( sin(x/b) )**2
            V[i]=v
    return V
def oscil(X):
    """
"oscil" :         Is the potential asociated to the 1D quantum oscilator
                  where the value of the potential is proportional to the
                  distance from the middle of the domain.

                  More details about harmonic oscilators can be found at:

                  http://en.wikipedia.org/wiki/Harmonic_oscillator

                  Parameters:
                 -----------
                 X:   numpy 1D array output of the Mesh1D function.

                 Returns:
                 --------
                 V: numpy array like vector with the values for the potential
                    over each node on X
                 
    """
def morse(X,D,lam):
    """
"morse":          Evaluates the Morse potential over the domain given by D
                  with the constants D_e= V0 and a=lam. A partial description
                  of this particular potential can be found at:

                  http://en.wikipedia.org/wiki/Morse_potential
                 
                 Parameters:
                 -----------
                 V0:   Constant with Well depht
                 lam:  constant proportional to the well widht
                 X:   numpy 1D array output of the Mesh1D function. 
    """





        
