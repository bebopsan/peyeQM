# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 17:31:09 2011

@author: santiago
"""

##module utils

"""
    This module contains small functions to be called by the main file
    in order to save space and improove readability.
"""

__all__ = ['create_points', 'calculate_area', 'substract_1', \
            'bloch_multiplication', 'bloch_sum' ]

__author__ = "Santiago Echeverri"

from numpy import zeros

def create_points(nodes, triangles, el):
    """
        Extracts the coordinates of the nodes in each triangular element and
        retrieves each vertice as a vector of two dimensions.
        
        Parameters:
        -----------
        nodes:  numpy array of dimension (n_nodes, 3), where n_nodes is the 
                number of nodes forming the mesh, and the three columns 
                represent coordinates (x, y, z).
        
        triangles:  numpy array of dimension (n_triangles, 4), where n_nodes 
                    is the number of nodes forming the mesh, the first is the 
                    label that refers to the physical entity where each node 
                    belongs. The remaining columns  tell which nodes belong
                    to each of the vertices of the triangle.
                    
        el:   integer value of the current element
        
        Returns:
        --------
        pt_a, pt_b, pt_c:  vector like array of two dimensions where the first 
                           is the position in x, and the second is the position 
                           of the node in y.  
    """
    pt_a = zeros(2)
    pt_b = zeros(2)
    pt_c = zeros(2)
    
    pt_a[0] = nodes[triangles[el, 1] - 1, 0]
    pt_a[1] = nodes[triangles[el, 1] - 1, 1]
    pt_b[0] = nodes[triangles[el, 2] - 1, 0]
    pt_b[1] = nodes[triangles[el, 2] - 1, 1]
    pt_c[0] = nodes[triangles[el, 3] - 1, 0]
    pt_c[1] = nodes[triangles[el, 3] - 1, 1]
    
    return pt_a, pt_b, pt_c

def calculate_area(lines):
    """
        Calculates the area of a triangle by using Heron's theorem
        
        Parameters:
        -----------
        line:  list of 3 2D arrays. Defines line segments of the triangle

        Returns:
        --------
        area:  float number.          
    """
    from math import sqrt
    from numpy.linalg import norm

    mg_ab = norm(lines[0]) # Calculate the magnitude (mg) of each line
    mg_bc = norm(lines[1])
    mg_ca = norm(lines[2])
        
    s = (mg_ab + mg_bc + mg_ca)/2    # Define the semiperimeter (s)
                                     # for each triangle
   
    area = sqrt(s*(s-mg_ab)*(s-mg_bc)*(s-mg_ca))
    return area
    
def substract_1(matrix):
   
    n_rows = matrix.shape[0]
    n_cols = matrix.shape[1]
    for i in range(n_rows):
        for j in range(1, n_cols):
            matrix[i,j] = matrix[i,j] - 1
    return matrix
    
    
def bloch_multiplication(k_x, k_y, nodes, ref_im, *matrices):
     """
     This function multiplies a given matrix 
     (in a future a given set of matrices), multiplies each Bloch boundary 
     node, by the phase factor correspondent with it's position.
     
     Parameters:
     -----------
     k_x:       Current value of the x component from the wavenumber vector
     
     k_y:       Current value of the y component from the wavenumber vector
     
     nodes:     Numpy array like matrix of node coordinates (n_nodes,3) where 
                coorsd(x,:)= x,y,z components of the node.
     
     ref_im:    A list of 2-column numpy arrays. Each array in the list 
                'ref_im', has in it's first column the reference node and on it's 
                second column the image node for that particular reference node.
     
     matrices:  Matrices to be operated.
     
     Returns:
     --------
     
     matrices:  The input matrix with all the phase multiplication operations 
                performed.
     """
     from cmath import exp
     
     # For each bloch condition in ref_im 
     im = ref_im[:,1]
     ref = list(set(list(ref_im[:, 0])))
     for i in im:
        x_im = nodes[ i - 1, 0]
        y_im = nodes[ i - 1, 1]
        
        fi = exp(1.0j*k_x*x_im)*exp(1.0j*k_y*y_im)
        for matrix in matrices:
            # Multiply the column of the image node by the phase factor
            matrix[:, i - 1] = fi * matrix[:, i - 1] 
            # Multiply the column of the image node by the comlex conjugate
            #phase factor            
            matrix[i - 1, :] = fi.conjugate() * matrix[:, i - 1]
     for i in ref:        
        x_ref = nodes[ i - 1, 0]
        y_ref = nodes[ i - 1, 1]
        
        ff = exp(1.0j*k_x*x_ref)*exp(1.0j*k_y*y_ref)
        for matrix in matrices:
           # and the same for the reference node:  
            matrix[:, i - 1] = ff * matrix[:, i - 1] 
            matrix[i - 1, :] = ff.conjugate() * matrix[:, i - 1]
     return matrices
     
def bloch_sum(ref_im, *matrices ):
    """
    This function takes the value of the image nodes in bloch periodicity 
    boundaries, and sums it to the value of the reference node.
    
    Parameters:
    -----------
    
    ref_im:    A list of 2-column numpy arrays. Each array in the list 
                'ref_im', has in it's first column the reference node and on it's 
                second column the image node for that particular reference node.
    
        
    matrices:  Matrices to be operated.
    
    Returns:
    -------
    
    matrices:  The input matrices but with the sums performed 
               and the image nodes columns and rows removed. 
    
    """
    from numpy import delete
    remove = []
    
    for k in range(len(ref_im[:, 0])):
        i = ref_im[k, 0]
        j = ref_im[k, 1]
        for matrix in matrices:
            # Sum image node row to reference node row 
            matrix[i - 1, :] = matrix[i - 1, :]+ \
                               matrix[j - 1, :]
            # Sum image node column to reference node column
            matrix[:, i - 1] = matrix[:, i - 1]+ matrix[:, j -1]
            #== stack the values of nodes in vertices for further removal===
            remove.count(j-1)       
            if remove.count(j-1) == 0:
                remove.append(j-1)
    remove.sort()
    for matrix in matrices:
        matrix = delete(matrix, remove, 0)
        matrix = delete(matrix, remove, 1)
    return matrices