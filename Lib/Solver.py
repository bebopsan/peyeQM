## module Solver
# -*- coding: utf-8 -*-
"""
Module Solver recieves data from the preprocessing stage and performs
the procesor stage of the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.
"""

__all__=['Schroedinger']
__author__='Santiago Echeverri Chacón'

from ReadMesh import*
from Write import*
from numpy import zeros,shape
from scipy import linalg
from math import sin,pi,exp, cosh,sqrt
import matplotlib.pyplot as plt

def Schroedinger(File,Nodes=0,Elems=0,parameter=[],Dimension=1,BCType='Dir'\
                 ,Type='Stationary',Eq='Schro',AnalisisParam=['y','y',4,4]):
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
                    Nodes is an array like matrix of dimension (Nes,3).

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

        Dimension:  String parameter that tells the program wether to solve for a
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


	Last modification: date 18/10/2011
    """
#------------------------ Load from file if given -----------------------------------
    if File!='':      # If type is not blank 
        Nodes,Elems=ReadMesh(File)                  # Import nodes and elements
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
    if '1' in Dimension:
        
        #---------------------- With dirichlet boundary conditions ----------------
        if 'Dir' in BCType and 'Stationary' in Type: 

            N=shape(Nodes)[0]   
                        
            # Initialization of the equivalent stiffness matrix
            K=zeros((N,N))
            # Initialization of the equivalent mass matrix
            M=zeros((N,N))
            
            # Reinterpretation of the parameter as the potential
            V=parameter
            
            # Value for the distance between the first 2 nodes and the last 2 nodes
            
            Li=abs(Nodes[1,0]-Nodes[0,0])
            Lf=abs(Nodes[N-1,0]-Nodes[N-2,0])
            
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

            # These two matrices are called the dirichlet matrices
            # of both the stiffness equivalent and mass equivalent matrices.
            
            Kd=K[1:N-1,1:N-1]
            Md=M[1:N-1,1:N-1]
            print 'K shape is:\n',K.shape

            if 'y'in AnalisisParam[0] and 'n' in AnalisisParam[1]:
                nVals=AnalisisParam[2]
                V=linalg.eigvalsh(Kd,Md,eigvals=(0,nVals-1))
                V=V/2
                print 'The Eigenvalues are:\n',V
                return V
                
            elif 'y'in AnalisisParam[1] and 'y'in AnalisisParam[1]:
                nVals=AnalisisParam[2]
                nVects=AnalisisParam[3]
                n=max(nVals,nVects)
                V,Dd=linalg.eigh(Kd,Md,eigvals=(0,n-1))
                V=V/2
                D=zeros((N,nVects))
                D[1:N-1,:]=Dd
                return V,D

            elif 'n'in AnalisisParam[1] and 'y'in AnalisisParam[1]:
                nVects=AnalisisParam[3]
                V,Dd=linalg.eigh(Kd,Md,eigvals=(0,nVects-1))
                D=zeros((N,nVects))
                D[1:N-1,:]=Dd
                return D
            else:
                print 'error: If you dont want anything why do you solve?'

        #------------------ End Dirichlet conditions----------------------------
                
        else:
            print 'Only Dirichlet for now sorry'
    #------------------------- End 1D problem ------------------------------------
    else:
        print 'only 1D for now. Sorry'
    
           
            
            














            
