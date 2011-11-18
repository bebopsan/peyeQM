# -*- coding: utf-8 -*-
## module matrices

"""
Created on Fri Nov 11 09:18:42 201
@author: santiago
"""
from numpy import zeros, dot, array

#======================== Local Stiffness matrix 2D P1 ========================
def local_stiffness_matrix(lines):
    n = len(lines)
    lo_stiff =  zeros((n, n))  # local (lo) stiffness (stiff) matrix
    for i in range(n):
        for j in range(n):
            lo_stiff[i, j] = dot(lines[i], lines[j]) # / (4*area) 
            
    return lo_stiff
#======================== Local Mass matrix  2D P1 ===========================
def local_mass_matrix():
    lo_mass = 1./24*array([[2, 1, 1],\
                           [1, 2, 1],\
                           [1, 1, 2]])
    return lo_mass
    
#======================== Local Newman matrix  1D P1 ===========================
def local_newman_matrix():
    lo_new = 1./6* array([[2, 1],\
                          [1, 2]])
    return lo_new
#======================== Local potential matrix  2D P1 =======================
def local_potential_matrix():
    lo_v = 1./60*array([[3, 1, 1],\
                          [1, 3, 1],\
                          [1, 1, 3]]) # Local (lo) potential (v) matrix
    return lo_v
#========================= Global Stiffness matrix 2D P1 ======================

def global_stiffness_matrix(nodes, triangles):
    """
        Funtion for the ensamble of the global stiffnes matrix.
        
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
        remove:   List with the numbers of columns to be removed

        Returns:
        --------
        
        glo_stif:  Matrix like array of size (n_nodes, n_nodes). Represents
                   the global stiffness matrix of the system
    """
    from matrices import local_stiffness_matrix
    from utils import create_points, calculate_area
    
    
    n_nodes = nodes.shape[0]         # Number (n) of nodes
    n_triangles = triangles.shape[0] # Number (n) of triangles
    
    glo_stif = zeros((n_nodes, n_nodes))# Initiate Global (glo) stiffness matrix
    
    for el in range(n_triangles):
         
        pt_a, pt_b, pt_c = create_points(nodes, triangles, el)
               
        ln_ab = pt_b - pt_a    
        ln_bc = pt_c - pt_b
        ln_ca = pt_a - pt_c    
        lines = [ln_ab, ln_bc, ln_ca]
        
        area = calculate_area(lines)
        
        lo_stif = local_stiffness_matrix(lines)
        lo_stif = lo_stif/(4*area)
       
        for i in range(1, lo_stif.shape[0]+1):
            for j in range(1, lo_stif.shape[0]+1):
                
                glo_stif[triangles[el, i]-1, triangles[el, j]-1] = \
                glo_stif[triangles[el, i]-1, triangles[el, j]-1] + \
                lo_stif[i-1,j-1]
                #                        ^ due to python's 0 base numbering
    return glo_stif

#========================= Global Mass matrix 2D P1 =========================

def global_mass_matrix(nodes, triangles, *remove):
    """
        funtion for the ensamble of the global mass matrix
        
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
        remove:   List with the numbers of columns to be removed
        Returns:
        --------
    """
    from matrices import local_mass_matrix
    from utils import create_points
    from numpy.linalg import det
    from numpy import delete
    
    n_nodes = nodes.shape[0]         # Number (n) of nodes
    n_triangles = triangles.shape[0] # Number (n) of triangles
    
    glo_mass = zeros((n_nodes, n_nodes))# Initiate Global (glo) mass matrix
    tr_mat = zeros((2,2))      # Initiate the transformation (tr) matrix (mat)
    
    for el in range(n_triangles):
         
        pt_a, pt_b, pt_c = create_points(nodes, triangles, el)
        tr_mat[:, 0] = pt_b - pt_a
        tr_mat[:, 1] = pt_c - pt_a
        
        jac = det(tr_mat)  
        
        lo_mass = local_mass_matrix()
        lo_mass = jac * lo_mass    # Escalate acording with the jacobian of the
                                   # transformation 
       
        for i in range(1, lo_mass.shape[0]+1):
            for j in range(1, lo_mass.shape[0]+1):
                
                glo_mass[triangles[el, i]-1, triangles[el, j]-1] = \
                glo_mass[triangles[el, i]-1, triangles[el, j]-1] + \
                lo_mass[i-1,j-1]
                #                        ^ due to python's 0 base numbering
    glo_mass_d = delete(glo_mass, remove, 0)
    glo_mass_d = delete(glo_mass_d, remove, 1)
    
    return glo_mass_d
#========================= Global potential matrix 2D P1 ======================

def global_potential_matrix(nodes, triangles, v, *remove):
    """
        funtion for the ensamble of the global potential matrix
        
        Parameters:
        -----------
        nodes:      numpy array of dimension (n_nodes, 3), where n_nodes is the 
                    number of nodes forming the mesh, and the three columns 
                    represent coordinates (x, y, z).
                    
        triangles:  numpy array of dimension (n_triangles, 4), where n_nodes 
                    is the number of nodes forming the mesh, the first is the 
                    label that refers to the physical entity where each node 
                    belongs. The remaining columns  tell which nodes belong
                    to each of the vertices of the triangle
                    
        v:          Vector with the potential associated to each node. 
        
        remove:     List with the numbers of columns to be removed
        
        Returns:
        --------
    """
    from matrices import local_potential_matrix
    from utils import create_points
    from numpy.linalg import det
    from numpy import delete
    
    n_nodes = nodes.shape[0]         # Number (n) of nodes
    n_triangles = triangles.shape[0] # Number (n) of triangles
    
    glo_v = zeros((n_nodes, n_nodes))# Initiate Global (glo) stiffness matrix
    vec_v = zeros((3,1))    
    tr_mat = zeros((2,2))      # Initiate the transformation (tr) matrix (mat)
    for el in range(n_triangles):    
        
        vec_v[0, 0] = v[triangles[el, 1]-1]
        vec_v[1, 0] = v[triangles[el, 2]-1]
        vec_v[2, 0] = v[triangles[el, 3]-1]
              
        pt_a, pt_b, pt_c = create_points(nodes, triangles, el)
        tr_mat[:, 0] = pt_b - pt_a
        tr_mat[:, 1] = pt_c - pt_a
        
        jac = det(tr_mat)  
        
        lo_v = local_potential_matrix()
        lo_v = jac * (lo_v * vec_v)  # Escalate acording with the jacobian 
                                       # of the transformation 
        
        
        
        for i in range(1, lo_v.shape[0]+1):
            for j in range(1, lo_v.shape[0]+1):
                
                glo_v[triangles[el, i]-1, triangles[el, j]-1] = \
                glo_v[triangles[el, i]-1, triangles[el, j]-1] + \
                lo_v[i-1,j-1]
                #                        ^ due to python's 0 base numbering
    glo_v_d = delete(glo_v, remove, 0)
    glo_v_d = delete(glo_v_d, remove, 1)
    
    return glo_v_d