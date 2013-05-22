## module Solver
# -*- coding: utf-8 -*-
"""
Module PostPro recieves data from the processing stage and performs
operations such as visualization of the values and figures representing
the output from the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.

"""

__all__=['plot_1d']
__author__='Santiago Echeverri Chacón'

from meshUtils import meshPlot
from read_mesh import read_vtk
import matplotlib.pyplot as plt
from math import pi

def plot_1d(filename, nodes = 0,elements = 0, parameter = [], bc_type = 'Dir', sol = 0):
    """
        Function plot 1D is intended as a tool for visualizing the results
        from the procesing stage of a 1D Finite Difference Method solver.

        For Schroedinger equation the plot would represent either the
        probability amplitude distribution function or the energy.

        Parameters:
        -----------

        filename:	    String with the name of the file contaoteining the
                    mesh and values.

        nodes:      Numpy array matrix of nodes containing the coordinates of
                    the nodes from the discretized domain.
                    nodes is an array like matrix of dimension (n_nodes,3).

                    Where each column represents the value of the nodes on
                    one of the three coordinate aNodeses Nodes,y,z.

        elements:      Numpy array matrix of elements containing the relations
                    between nodes from the discretized domain.
                    elements is an array like matrix of dimension (n_elements,2)
                    in the case of 1D, and (n_elements,3) for 2D problems.
                    
        parameter:  Is an array that describes the potential actuating over the
                    the elements of the domain given by elements. For each element in
                    Elems there is an associated potential value on the same
                    position in the array parameter.

                    The potential in Scroedinger equation defines the specific
                    nature of the problem to be solved. For more details on how
                    to define a potential and what does it mean pleas read the
                    documentation of the Potential1D function in the module PrePro.

         bc_type:    String parameter for the selection of a border condition
                    that can be either:

                        'Dir'   For the Dirichlet border condition
                                (Infinite potential well).

                        'Bloch' For the periodic formulation of the problem.
                                (Electron in a periodic material )
        Last modification: date 25/10/2011
    """

    
#------------------------ Load from file if given -----------------------------------
    if filename != '':      # If type is not blank
        print filename
        nodes, elements, v, sol = read_vtk(filename)
        
        n_vals = sol.shape[1]
        n = nodes.shape[0]
        nodes = nodes[:, 0]
        
        elements = elements[:, 1:]
        if bc_type == 'Dir':
            # -----------------------  Potential plot -----------------------------
            
            plt.subplot(121)
            plt.plot( (nodes[0:n-1]+nodes[1:n])/2., v)
            plt.xlabel('x position inside the well')
            plt.ylabel('Magnitude of the potential')
            plt.title('Potential $\hat{V}$')
            
            # ------------------------ Eigenvectors plot -----------------------
            
            plt.subplot(122)
            legend = []
            for i in range(0, n_vals):
                plt.plot(nodes, sol[:, i])
                legend.append('n = '+str(i+1))

            plt.legend(legend)
            plt.title('Wave functions $\Psi_n$')
            plt.xlabel('x position inside the well')
            plt.ylabel('Amplitude of the probability function $\psi$')
            plt.show()
            
    n_vals = sol.shape[1]
    n = nodes.shape[0]
    nodes = nodes[:, 0]
    elements = elements[:, 1:]
    v = parameter
    if bc_type == 'Dir':
        # -----------------------  Potential plot -----------------------------
        plt.figure(figsize=(11, 4))
        plt.subplot(121)
        plt.plot( (nodes[0:n-1]+nodes[1:n])/2., v)
        plt.xlabel('x position inside the well')
        plt.ylabel('Magnitude of the potential')
       
        plt.xlim(xmax=2*pi)
        plt.title('Potential $\hat{V}$')
        
        # ------------------------ Eigenvectors plot -----------------------
        
        plt.subplot(122)
        
        legend = []
        for i in range(0, n_vals):
            plt.plot(nodes, sol[:, i])
            legend.append(str(i+1))
    
        plt.legend(legend,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title('Wave functions $\Psi_n$')
        plt.xlabel('x position inside the well')
        plt.xlim(xmax=2*pi)
        plt.ylabel('Amplitude of the probability function $\psi$')
        plt.subplots_adjust(wspace = 0.4)
        plt.show()    















    
