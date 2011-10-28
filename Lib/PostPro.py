## module Solver
# -*- coding: utf-8 -*-
"""
Module PostPro recieves data from the processing stage and performs
operations such as visualization of the values and figures representing
the output from the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.

"""

__all__=['Plot1D']
__author__='Santiago Echeverri Chacón'

from meshUtils import meshPlot
from ReadMesh import ReadVTK
import matplotlib.pyplot as plt

def Plot1D(File,Nodes=0,Elems=0,parameter=[],BCType='Dir'):
    """
        Function plot 1D is intended as a tool for visualizing the results
        from the procesing stage of a 1D Finite Difference Method solver.

        For Schroedinger equation the plot would represent either the
        probability amplitude distribution function or the energy.

        Parameters:
        -----------

        File:	    String with the name of the file containing the
                    mesh and values.

        Nodes:	    Numpy array matriNodes of nodes containing the coordinates of
                    the nodes from the discretized domain.
                    Nodes is an array like matriNodes of dimension (nNodes,3).

                    Where each column represents the value of the nodes on
                    one of the three coordinate aNodeses Nodes,y,z.

        Elems:      Numpy array matriNodes of elements containing the relations
                    between nodes from the discretized domain.
                    Elems is an array like matriNodes of dimension (nElems,2)
                    in the case of 1D, and (nElems,3) for 2D problems.
                    
        parameter:  Is an array that describes the potential actuating over the
                    the elements of the domain given by Elems. For each element in
                    Elems there is an associated potential value on the same
                    position in the array parameter.

                    The potential in Scroedinger equation defines the specific
                    nature of the problem to be solved. For more details on how
                    to define a potential and what does it mean pleas read the
                    documentation of the Potential1D function in the module PrePro.

         BCType:    String parameter for the selection of a border condition
                    that can be either:

                        'Dir'   For the Dirichlet border condition
                                (Infinite potential well).

                        'Bloch' For the periodic formulation of the problem.
                                (Electron in a periodic material )
        Last modification: date 25/10/2011
    """

    
#------------------------ Load from file if given -----------------------------------
    if File !='':      # If type is not blank
        print File
        Nodes,Elems,V,Sol=ReadVTK(File)
        
        nVals=Sol.shape[1]
        N=Nodes.shape[0]
        Nodes=Nodes[:,0]
        
        Elems=Elems[:,1:]
        if BCType=='Dir':
            # -----------------------  Potential plot -----------------------------
            plt.figure(1)
            plt.plot( (Nodes[0:N-1]+Nodes[1:N])/2. ,V)
            plt.figure(1).suptitle('Potential')
            
            # ------------------------ Eigenvectors plot -----------------------
            plt.hold(True)
            plt.figure(2)
            legend=[]
            for i in range(0,nVals):
                plt.plot(Nodes,Sol[:,i])
                legend.append('n = '+str(i+1))

            plt.legend(legend)
            plt.show()
        
        















    
