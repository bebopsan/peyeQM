## module PrePro
# -*- coding: latin-1 -*-
''' This module contains the necessary functions for the preprocessing
    phase of finite element procedure for the solution of differentia
    equations.
'''

__all__=['Mesh1D','meshPlot','Potential1D','well','finwell','xwell'\
          ,'Mesh2D','poschl','oscil']
__author__="Santiago Echeverri Chacón"
from math import sin,pi,exp, cosh,sqrt
from numpy import zeros,array,size
from meshUtils import*

# ---------------- Mesh construction related functions  -----------------#
def Mesh1D(Type,xmin,xmax,base=3,N=100):
    '''
        In the mesh function the user can decide wether to use an
        homogeneus discretization of the line of lenght [xmin,xmax] in N elements,
        or a logarithmically spaced discretization with more density of points in
        a certain point.

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
        
        xmin:  coordinate of the beginning of the line segment
        xmax:  coordinate of the end of the line segment
        N:     Number of divisions
        base:  For the case of a nonhomogeneus mesh, the parameter 'base' tells the
               base for the logarithm used to distribute coordinates.
               for base --> 1 like 1.0001 the logmesh tends to a equally
                distributed mesh

        The default value  for  N is N=100.
        The user can define a new value by assigning them as arguments of the
        function e.g. N=100000.

        Returns:
        --------
        Nodes:  numpy array like vector of the discretized domain with lenght N
        Elems:  numpy array like vector of the relations between nodes.

        Raises:
        -------
        Error if the user gives a wrong Type argument
        
        Last modification: 6/Oct/2011
    '''

    if 'simple' in Type:
        Nodes,Elems=mesh1D(xmin,xmax,N)
        return Nodes,Elems
    elif 'left' in Type:
        Nodes,Elems=logmesh1D(xmin,xmax,N,base,'L')
        return Nodes,Elems   
    elif 'right' in Type:
        Nodes,Elems=logmesh1D(xmin,xmax,N,base,'R')
        return Nodes,Elems
        
    elif 'center' in Type:
        Nodes,Elems=logmesh1D(xmin,xmax,N,base,'C')
        return Nodes,Elems
    else:
        print 'You entered a wrong Type parameter. Please try again, or read the\
              documentation'
        
def Mesh2D(Type,xmin,xmax,ymin,ymax,Nx=100,Ny=100,basex=3,basey=3,optionx='',optiony=''):
    '''
        Mesh2D works like Mesh 1D but with a rectangular domain rather than a
        line segment. The mesh is built using triangular elements, and the nodes
        can be arranged in a homogeneus distribution or a logarithm-like non
        homogeneus distribution for each coordinate. 

        Parameters:
        ----------
        Type:  Is a string variable that will tell the kind
                of Mesh to be made.

		"simple":   Will execute the discretization using
                            Nx divisions over the x side of lenght [xmin,xmax]
                            and Ny divisions over y side of lenght [ymin,ymax] 
		"center":   Creates a 2D domain with values spaced evenly
                            on a log scale with more density of values on
                            the middle. The user can define the base of the log
                            function by modifiying the parameters basex and basey
                "custom":   lets the user decide the orientation of the 2D logarithmic
                            mesh over each axis. This parameter enables full
                            customization over the parameters
        
        xmin:  initial value of the rectangular domain over x axis
        xmax:  final value of the rectangular domain over x axis
        ymin:  initial value of the rectangular domain over y axis
        ymax:  final value of the rectangular domain over y axis
        Nx:     Number of divisions over x
        Ny:     Number of divisions over y
        basex: the same base parameter but for the x coordinates of the rectangle
        basey: the same base parameter but for the y coordinates of the rectangle
        optionx:  Lets the user decide the orientation of the log distribution
                  for the custom mesh over the x axis. thi parameter can be
                  L for more points on the left
                  R for more points on the right
                  C for more points on the center
         optiony: Lets the user decide the orientation of the log distribution
                  for the custom mesh over the y axis. thi parameter can be
                  L for more points on the top
                  R for more points on the right
                  C for more points on the bottom
                    

        Returns:
        --------
        Nodes:  numpy array like matrix of the discretized domain with shape Nx Ny
        Elems:  numpy array like matrix of the relations between nodes.

        Raises:
        -------
        Error if the user gives a wrong Type argument
        
        Last modification: 6/Oct/2011
    '''

    if 'simple' in Type:
        Nodes,Elems= meshtr2D(xmin,xmax,ymin,ymax,Nx,Ny)
        return Nodes,Elems
    elif 'center' in Type:
        Nodes,Elems=logmeshtr2D(xmin,xmax,ymin,ymax,Nx,Ny,'C','C',basex,basey)
        return Nodes,Elems
    elif 'custom' in Type:
        if optionx ==''or optiony =='':
            print 'Please enter the custom parameters optionx and optiony.'
        else:
            Nodes,Elems=logmeshtr2D(xmin,xmax,ymin,ymax,Nx,Ny,optionx,optiony,basex,basey)
            return Nodes,Elems
    else:
        print 'You entered a wrong Type parameter. Please try again, or read the\
              documentation'

# ----------------  Potential related functions ---------------------------#
def Potential1D(Type,X,V0=1,Vleft=0.0,Vright=0.0,width=1.,lam=.2,D=1):
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
   
    error_flag=0
    
    if Type== 'well':
       
        V=well(X,V0)
  
    elif 'finwell' in Type:
        V=finwell(X,V0,Vleft,Vright,width)
        
    elif 'xwell' in Type:
        V=xwell(X,lam)
        
    elif "poschl" in Type:
        V=poschl(X,V0,lam)
        
    elif "oscil" in Type:
        V=oscil(X)
        
    elif "morse" in Type:
        V=morse(X,D,lam)
        
    else:
         print 'You entered a wrong Type parameter. Please try again, or read the\
              documentation'
         error_flag=1
    if not error_flag:
        z=zeros((V.size,1))
        z[:,0]=V
        V=z
        
        return V
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
    print 'V.size', V.size
def finwell(X,V0,Vleft,Vright,width):
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
    print a
    for i in range(0,N-1):
    	if  X[i]<(a/2.-width/2.):
		V[i]= Vleft
	elif X[i]>(a/2.+width/2.) :
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
"poschl":       Potential Well defined by G.Pï¿½schl and Edward Teller
                  described in:
                    http://en.wikipedia.org/wiki/Pï¿½schl-Teller_potential

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
    N=size(X)
    V=zeros(N-1)
    a=X[N-1]
    for i in range(0,N-1):
        x = (X[i]+X[i+1])/2.
        v = (x-a/2)**2
        V[i]=v

   
    return V
    
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
    N=size(X)
    V=zeros(N-1)
    a=X[N-1]
    xe = a/2.
    for i in range(0,N-1):
        x = (X[i]+X[i+1])/2.
        v = D*(1. - exp(-lam*(x-xe)) )**2
        V[i]=v
    return V
    
# Other --------------------------------------------------------------------

def meshPlot(Nodes,Elems,facecolor,coordnum,elemnum):

    if(Elems.shape[1]==2):
        y = Nodes*0
        plt.figure()
        plt.plot(Nodes,y,'-k')
        plt.plot(Nodes,y,'ob')
        plt.axis('equal')
        plt.show()
    else:
        plt.figure()
        plt.hold(True)
        for i in range(0,Elems.shape[0]):
            xcoord = zeros( (size(Elems[i,:])))
            ycoord = zeros( (size(Elems[i,:])))
            for j in range(0,size(Elems[i,:])):
                xcoord[j] = Nodes[Elems[i,j],0]
                ycoord[j] = Nodes[Elems[i,j],1]
            plt.fill(xcoord,ycoord,facecolor=facecolor,alpha=1.0, edgecolor='k')
        plt.plot(Nodes[:,0],Nodes[:,1],'ob')
        plt.axis('equal')
        plt.axis('off')
        plt.show()



        
