## module Solver
# -*- coding: utf-8 -*-
"""
Module Solver recieves data from the preprocessing stage and performs
the procesor stage of the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.
"""

__all__=['Schroedinger','matAssembly1D']
__author__=['Santiago Echeverri Chacón','Nicolas Guarin']

from ReadMesh import*
from Write import*
from numpy import zeros,shape,linspace,cfloat,array
from scipy import linalg
from math import sin,pi,exp, cosh,sqrt
import matplotlib.pyplot as plt

def Schroedinger(File,Nodes=0,Elems=0,parameter=[],Dimension=1,BCType='Dir'\
                 ,Type='Stationary',Eq='Schro',AnalisisParam=['y','y',4,4,101]):
    """
        Function Schroedinger computes the solution for Schroedinger
        equation using the Finite Element Method.

        For now it is able to solve the Stationary form with Dirichlet and
        Bloch boundary conditions in 1D.

        The function is made to be used either with direct assignment of
        the input parameters, or by reading those parameters from a file.

        Parameters:
        -----------

        File:        String with the name of the file that contains the
                     information regarding the geometry, mesh, border
                     conditions, and other parameters necessary for the
                     solution of the problem.

                     This file is a modified Gmsh output file with extension
                     .msh

        Nodes:	    Numpy array matrix of nodes containing the coordinates of
                    the nodes from the discretized domain.
                    Nodes is an array like matrix of dimension (nNodes,3).

                    Where each column represents the value of the nodes on
                    one of the three coordinate axes x,y,z.

        Elems:      Numpy array matrix of elements containing the relations
                    between nodes from the discretized domain.
                    Elems is an array like matrix of dimension (nElems,2)
                    in the case of 1D, and (nElems,3) for 2D problems.
                    
        parameter:  Is an array that describes the potential actuating over the
                    the elements of the domain given by Elems. For each element in
                    Elems there is an associated potential value on the same
                    position in the array parameter.

                    The potential in Scroedinger equation defines the specific
                    nature of the problem to be solved. For more details on how
                    to define a potential and what does it mean pleas read the
                    documentation of the Potential1D function in the module PrePro.

        Dimension:  int parameter that tells the program wether to solve for a
                    1D problem or a 2D problem (not supported yet)

        BCType:     String parameter for the selection of a border condition
                    that can be either:

                        'Dir'   For the Dirichlet border condition
                                (Infinite potential well).

                        'Bloch' For the periodic formulation of the problem.
                                (Electron in a periodic material )

        Type:       String that tells wether to solve the stationary version of
                    the equation or another not yet suported.

                    'Stationary'   

        AnalisisParam:   Array that contains the information regarding the number
                         of solutions to be computed and wether to save the values
                         or not.

                        AnalisisParam[0]:  String  answer to the question
                                                   save  Eigen Values?
                        AnalisisParam[1]:  String  answer to the question
                                                   save  Eigen Vectors?
                        AnalisisParam[2]:  Integer  number of Eigen Values to save
                        AnalisisParam[3]:  Integer  number of Eigen Vectors to save


	Last modification: date 25/10/2011
    """
#------------------------ Load from file if given -----------------------------------
    if File!='':      # If type is not blank 
        Nodes,Elems=Readmsh(File)                  # Import nodes and elements
        SolverInput=ReadSolverInput(File)           # Import input parameters for
                                                    # the solution    
        if SolverInput !=[]:
            Dimension=SolverInput[0]
            BCType=SolverInput[1]
            parameter=SolverInput[2]                # Reasign
            Eq=SolverInput[3]
            Type=SolverInput[4]
            AnalisisParam=SolverInput[5]

        # If the file does not contain solution parameters, the ones given as
        # arguments of the function are written in the file.
        else:
            WriteSolverInput(File,Dimension=Dimension,BCType=BCType,\
                             Type=Type,Eq=Eq,parameter=parameter)

    #------------------------------  1D problem ------------------------------------       
    
    if Dimension==1:
        N=shape(Nodes)[0]
        # Reinterpretation of the parameter as the potential
        V=parameter
        #---------------------- With dirichlet boundary conditions ----------------
        if 'Dir' in BCType and 'Stationary' in Type: 

               
            K,M = matAssembly1D(Nodes,N)                               

            # These two matrices are called the Dirichlet matrices
            # of both the stiffness equivalent and mass equivalent matrices.
            
            Kd=K[1:N-1,1:N-1]
            Md=M[1:N-1,1:N-1]
            
            print 'K shape is:\n',K.shape

            if 'y'in AnalisisParam[0] and 'n' in AnalisisParam[1]:
                nVals=int(AnalisisParam[2])
                V=linalg.eigvalsh(Kd,Md,eigvals=(0,nVals-1))
                V=V/2
                print 'The Eigenvalues are:\n',V
                return V
                
            elif 'y'in AnalisisParam[0] and 'y'in AnalisisParam[1]:
                nVals=int(AnalisisParam[2])
                nVects=int(AnalisisParam[3])
                n=max(nVals,nVects)
                V,Dd=linalg.eigh(Kd,Md,eigvals=(0,n-1))
                V=V/2
                D=zeros((N,nVects))
                D[1:N-1,:]=Dd
                return V,D

            elif 'n'in AnalisisParam[0] and 'y'in AnalisisParam[1]:
                nVects=int(AnalisisParam[3])
                V,Dd=linalg.eigh(Kd,Md,eigvals=(0,nVects-1))
                D=zeros((N,nVects))
                D[1:N-1,:]=Dd
                return D
            else:
                print 'error: If you dont want anything why do you solve?'

        #------------------ End Dirichlet conditions----------------------------
                
        #----------------- With Bloch periodic boundary conditions -------------        
        elif 'Bloch' in BCType:
            
            import cmath

            K,M = matAssembly1D(Nodes,N)                   

            print 'K shape is:\n',K.shape

                
            # Bloch-Periodicity imposition

            xi=Nodes[0,0]   # initial x
            xf=Nodes[N-1,0]  # final x
                                    
            nVals =int(AnalisisParam[2])  # number of eigenvales to compute

            nk = int(AnalisisParam[4]) # number of k to sweep
            kmax = 4.*pi/Nodes[N-1,0]
            kmin = -0.0
            k_range = linspace(kmin, kmax, num=nk)
            omega = zeros( (len(k_range),nVals) )
            E = zeros( (len(k_range),nVals) )
                              
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
                    omega[ll,i] = sqrt( abs(vals[i]) )

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
        else:
            print 'Only Dirichlet for now sorry'
    #------------------------- End 1D problem ------------------------------------
    else:
        print 'only 1D for now. Sorry'
    
           
def matAssembly1D(Nodes,N):
    """
        Assembly the equivalent stiffness and equivalent mass matrices for a right to left
        numbered 1D mesh.

        Parameters:
        -----------


        Nodes:	    Numpy array matrix of nodes containing the coordinates of
                    the nodes from the discretized domain.
                    Nodes is an array like matrix of dimension (nNodes,3).

                    Where each column represents the value of the nodes on
                    one of the three coordinate axes x,y,z.

        N:          Number of degree of freedom

  	Last modification: date 27/10/2011
    """

    # Initialization of the equivalent stiffness matrix
    K=zeros((N,N))
    # Initialization of the equivalent mass matrix
    M=zeros((N,N))
    
    
    
    # Value for the distance between the first 2 nodes and the last 2 nodes
    
    Li=abs(Nodes[1,0]-Nodes[0,0])
    Lf=abs(Nodes[N-1,0]-Nodes[N-2,0])


    # Matrices assembly
    K[0,0]= 1/Li + Li*V[0]/3
    K[N-1,N-1]= 1/Lf + Lf*V[N-2]/3

    for i in range(0,N-1):
        L=abs(Nodes[i+1,0]-Nodes[i,0])       #Distance between nodes
        if i!=0:
            K[i,i]=2/L + L*(V[i-1]+V[i])/3   # Central diagonal
        K[i,i+1]= -1/L + L*V[i-1]/6          # Upper diagonal
        K[i+1,i]= -1/L + L*V[i-1]/6          # Lower diagonal
        
    M[0,0]=Li/3
    M[N-1,N-1]=Lf/3
    for i in range(0,N-1):
        L=abs(Nodes[i+1,0]-Nodes[i,0])      #Distance between nodes
        if i!=0:
            M[i,i]=2*L/3                    # Central diagonal
        M[i,i+1]=L/6                        # Upper diagonal  
        M[i+1,i]=L/6                        # Lower diagonal

    return K,M












            
