# -*- coding: utf-8 -*-
## module vectors
"""
Created on Fri Nov 11 10:36:29 2011

@author: santiago
"""
__all__ = ['dirichlet_vector','sources_vector', 'newman_vector', \
            'image_reference_bloch_vectors']

from numpy import zeros, dot, delete

#======================== Dirichlet conditions ===============================
def dirichlet_vector(bc_lines, dirichlet, glo_stif):
    """
        This function computes the d vector associated to Dirichlet
        boundary conditions on certain nodes of the domain. 
        
        Parameters:
        -----------
        bc_lines:  numpy array of dimension (n_lines, 3), where the first column
                   is the label that refers to the physical entity where each 
                   node belongs, and the other two columns represent node 1 and
                   2 of the line element
                
        dirichlet: list of lists. the first index refers to lists that describe
                   one particular Dirichlet condition, where the first elements 
                   of this list is the tag of the physical entity, and the 
                   second is an array that has the values over each degree of 
                   freedom.
                   
        glo_stiff: matrix like array of size (n_nodes,n_nodes) called
                   the stiffness matrix of the system.
                   
           
        Returns:
        --------
        glo_stif_d:  global stiffness matrix reduced by removing the rows and 
                     columns asociated with the dirichlet conditions.
        d:  
        remove:   List with the numbers of columns to be removed.
        
    """ 
    n_nodes = glo_stif.shape[0]         # Number (n) of nodes
    n_lines = bc_lines.shape[0]         # Number of lines
    n_bc_d = len(dirichlet)  # Number of Dirichlet boundary conditions
    remove = []              # This vector will tell which rows and columns
                             # should be removed 
    g = zeros(n_nodes)
    for bc in range(n_bc_d):
        tag = dirichlet[bc][0]
        value = dirichlet[bc][1][0][0]
        
        for ln in range(n_lines):
            
            if bc_lines[ln, 0] == tag:
                g[bc_lines[ln, 1]-1] = value
                g[bc_lines[ln, 2]-1] = value
                
                if bc_lines[ln, 1]-1 not in remove: # makes a list for removing
                    remove.append(bc_lines[ln, 1]-1)# lines and columns.
                if bc_lines[ln, 2]-1 not in remove:
                    remove.append(bc_lines[ln, 2]-1)
        
    d = dot(glo_stif, g)
    d = delete(d, remove)
    
    glo_stif = delete(glo_stif, remove, 0)
    glo_stif = delete(glo_stif, remove, 1)
    
    return glo_stif, d, remove, g

#======================== Newman conditions ===============================
def newman_vector(nodes, bc_lines, newman, remove):
    """
        This function computes the d vector associated to Dirichlet
        boundary conditions on certain nodes of the domain. 
        
        Parameters:
        -----------
        bc_lines:  numpy array of dimension (n_lines, 3), where the first column
                   is the label that refers to the physical entity where each 
                   node belongs, and the other two columns represent node 1 and
                   2 of the line element.
                
        newman:   list of lists. the first index refers to lists that describe
                   one particular Newman condition, where the first elements 
                   of this list is the tag of the physical entity, and the 
                   second is an array that has the values over each degree of 
                   freedom.
                   
        remove:   List with the numbers of columns to be removed.
        
        n_nodes:  integer. represents the number of rows on the array of nodes.               
           
        Returns:
        --------
        q:             
    """ 
    from matrices import local_newman_matrix
    from numpy.linalg import norm
    
    n_lines = bc_lines.shape[0]         # Number of lines
    n_bc_d = len(newman)  # Number of Dirichlet boundary conditions
    n_nodes = nodes.shape[0]         # Number (n) of nodes
    q = zeros(n_nodes)
    
    for bc in range(n_bc_d):
        tag = newman[bc][0]
        value = newman[bc][1][0][0]
        
        for ln in range(n_lines):
            
            if bc_lines[ln, 0] == tag:
                q[bc_lines[ln, 1]-1] = value
                q[bc_lines[ln, 2]-1] = value
#    print q            
    for ln in range(n_lines):
        if bc_lines[ln, 0] == tag:
            vec_q = zeros(2)                 # q on nodes
            vec_q[0] = q[bc_lines[ln, 1]-1]  # fill heat element vec_q 
            vec_q[1] = q[bc_lines[ln, 2]-1]
                                 # calculate distance between nodes in the line       
            ln_h = nodes[bc_lines[ln, 2]-1][0:2] - nodes[bc_lines[ln, 1]-1][0:2]
            h =  norm(ln_h)
            lo_new = local_newman_matrix() # Call the matrix for interpolation
            vec_q =  h * dot(lo_new, vec_q)
            q[bc_lines[ln, 1]-1] = vec_q[0] # Reasign the interpolated 
            q[bc_lines[ln, 2]-1] = vec_q[1] # newman conditions
 #   print q
    q = delete(q, remove) # Remove Dirichlet columns
 #   print q
    return q
       
#======================== Sources conditions =============================== 
   
def sources_vector(nodes, triangles, remove):
    """
        This function computes the l vector associated to sources in the domain 
        
        Parameters:
        -----------
        nodes:  numpy array of dimension (n_nodes, 3), where n_nodes is the 
                number of nodes forming the mesh, and the three columns 
                represent coordinates (x, y, z).
        
        triangles:  numpy array of dimension (n_triangles, 4), where n_nodes 
                    is the number of nodes forming the mesh, the first is the 
                    label that refers to the physical entity where each node 
                    belongs. The remaining columns  tell which nodes belong
                    to each of the vertices of the triangle
        Returns:
        --------
        l:  right hand side vector of the equation associated to sources.   
    """    
    
    from utils import create_points
    from numpy.linalg import det
    from numpy import array
    
    n_nodes = nodes.shape[0]         # Number (n) of nodes
    n_triangles = triangles.shape[0] # Number (n) of triangles

    glo_l = zeros(n_nodes)   # Initiate global sources vector l 

    tr_mat = zeros((2,2))      # Initiate the transformation (tr) matrix (mat)
    
    for el in range(n_triangles):
        pt_a, pt_b, pt_c = create_points(nodes, triangles, el)
        
        tr_mat[:, 0] = pt_b - pt_a
        tr_mat[:, 1] = pt_c - pt_a
        
        jac = det(tr_mat)  
        
        lo_l = array([1./6, 1./6, 1./6]) # This is the local (lo) l vector fo 
                                         # a standard triangular element of 
                                         # base 1 and height 1
        for i in range(1, 4):
            glo_l[triangles[el, i]-1] = glo_l[triangles[el, i]-1] + \
                                        jac * lo_l[i - 1]            
    
    glo_l = delete(glo_l, remove) # Remove Dirichlet columns
    
    return glo_l    
    
 #============== build solution using g and the dirichlet conditions ==========

def build_solution(dir_solution, g, remove):
    n_sol = g.shape[0]
    j = 0
    for i in range(n_sol):
        if i not in remove:
            g[i] = dir_solution[j]
            j = j+1
    return g       
   
# Construct the image and reference vectors needed for imposing bloch conditions

def image_reference_bloch_vectors(bc_lines, bloch):
    """
    This function loads the lines of the boundary and the list that contains 
    bloch periodicity conditions, and builds a list of arrays that relate
    image nodes with their corresponding reference nodes.
    
    Paprameters:
    -----------
    
    bc_lines:   Array of line elements.

    bloch:      List of bloch conditions as read by a '.bc' file.__class__

    Returns:
    --------
    
    im_ref:     Each array in the list 'im_ref', has in it's first column the 
                image node and on it's second collumn the reference node for 
                that particular image node.
    """
    from numpy import zeros, copy
    n_bloch = bloch.shape[0]
    it_bloch = range(n_bloch)    
    n_lines = bc_lines.shape[0]
    print 'n_lines', n_lines
    #=========== Create a list of arrays ======================================
    im_ref = []    
    for i in it_bloch:    
        im_ref.append(zeros((n_lines/(2*n_bloch)+1,2), dtype = int))  
    #=============== Loop and fill the lists ================================== 
    count = 0
    i = 0
    j = 0
    
    while count == 0:
        
        for bl in it_bloch:
            
            
            if bc_lines[j, 0] == bloch[bl, 0]:
                im_ref[bl][i, 0] = bc_lines[j, 1]
                if i == n_lines/(2*n_bloch)-1:
                    im_ref[bl][i+1, 0] = bc_lines[j, 2]
                
            elif bc_lines[j, 0] == bloch[bl, 1]:
                im_ref[bl][i, 1] = bc_lines[j, 1]
                if i == n_lines/(2*n_bloch)-1:
                    im_ref[bl][i+1, 1] = bc_lines[j, 2]
        i = i+1
        j = j+1

	
        if i == n_lines/(2*n_bloch):
            i = 0
        if j >= n_lines:
            count = 1
    
    for bl in it_bloch:
        if bloch[bl, 2]*bloch[bl, 3] == -1:
            im_ref[bl][:,1] = copy(im_ref[bl][::-1,1])
    
    
    return im_ref
                
    
    
    
    
    
    
    
    
    
    
    
