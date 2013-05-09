## module solver
# -*- coding: utf-8 -*-
"""
Module Solver recieves data from the preprocessing stage and performs
the procesor stage of the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.
"""

__all__ = ['schroedinger','matAssembly1D']
__author__ = ['Santiago Echeverri Chacón','Nicolas Guarin']

from read_mesh import read_mesh, read_solver_input, read_bc
from write import write_solver_input
from numpy import zeros, shape, linspace, cfloat, array
from scipy import linalg
from math import sin, pi, exp, sqrt
import matplotlib.pyplot as plt

def schroedinger(filename, nodes = 0, elements = 0, parameter = [], \
                 dimension = 1, bc_type = 'Dir', sol_type = 'Stationary', \
                 eq = 'Schro', analysis_param = ['y', 'y', 4, 4, 20, 20, 2], \
                 bc_filename = ''):
    """
    Function Schroedinger computes the solution for Schroedinger
    equation using the Finite Element Method.

    For now it is able to solve the Stationary form with Dirichlet and
    Bloch boundary conditions in 1D.

    The function is made to be used either with direct assignment of
    the input parameters, or by reading those parameters from a file.

    Parameters:
    -----------

    filename:        String with the name of the file that contains the
                 information regarding the geometry, mesh, border
                 conditions, and other parameters necessary for the
                 solution of the problem.

                 This file is a modified Gmsh output file with extension
                 .msh

    nodes:	Numpy array matrix of nodes containing the coordinates of
                the nodes from the discretized domain.
                nodes is an array like matrix of dimension (n,3).

                Where each column represents the value of the nodes on
                one of the three coordinate axes x,y,z.

    elements:      Numpy array matrix of elements containing the relations
                between nodes from the discretized domain.
                elements is an array like matrix of dimension (n_elements,2)
                in the case of 1D, and (n_elements,3) for 2D problems.

    parameter:  Is an array that describes the potential actuating over the
                the elements of the domain given by elements. For each element 
                in elements there is an associated potential value on the same
                position in the array parameter.

                The potential in Scroedinger equation defines the specific
                nature of the problem to be solved. For more details on how
                to define a potential and what does it mean pleas read the
                documentation of the Potential1D function in the module 
                PrePro.

    dimension:  int parameter that tells the program wether to solve for a
                1D problem or a 2D problem.

    bc_type:    String parameter for the selection of a border condition
                that can be either:

                    'Dir'   For the Dirichlet border condition
                            (Infinite potential well).

                    'Bloch' For the periodic formulation of the problem.
                            (Electron in a periodic material )

    sol_type:   String that tells wether to solve the stationary version of
                the equation or another not yet suported.

                'Stationary'

    analysis_param:   Array that contains the information regarding the 
                    number of solutions to be computed and wether to save 
                    the values or not.

                    analysis_param[0]:  String  answer to the question
                                               save  Eigen Values?
                    analysis_param[1]:  String  answer to the question
                                               save  Eigen Vectors?
                    analysis_param[2]:  Integer number of Eigen Values 
                                               to save
                    analysis_param[3]:  Integer number of Eigen Vectors 
                                               to save
                    analysis_param[4]:  Integer number of wave numbers 
                                        in x to sweep
                    analysis_param[5]:  Integer number of wave numbers 
                                        in y to sweep
                    analysis_param[6]:  biggest value of k. it may be the lenght
                                        of the dommain
                    
	Last modification: date 16/11/2011
    """
    #------------------------ Load from file if given -------------------------
    assert isinstance(filename, str)     
    
    if filename != '':      # If type is not blank
        solver_input = read_solver_input(filename)        # Import input parameters for
        if solver_input != []:
            dimension = solver_input[0]
            bc_type = solver_input[1]
            parameter = solver_input[2]                # Reasign
            eq = solver_input[3]
            sol_type = solver_input[4]
            analysis_param = solver_input[5]         # the solution 
            bc_filename = solver_input[6]
        # If the file does not contain solution parameters, the ones given as
        # arguments of the function are written in the file.
        else:
           write_solver_input(filename, dimension = dimension, 
                              bc_type = bc_type, sol_type = sol_type, \
                              eq = eq, parameter = parameter, \
                              bc_filename = bc_filename)
        
          
        nodes, elements = read_mesh(filename)        # Import nodes and elements
        n = shape(nodes)[0]                 #Number of nodes      
#        el_points = elements[0]
        bc_lines = elements[1]      # Boundary condition is (bc) lines
        if dimension == 2:
            triangles = elements[2]
            
   
    #------------------------------  1D problem -------------------------------
 
    # Reinterpretation of the parameter as the potential

    v = parameter

    if dimension == 1:
        
        K, M = mat_assembly_1d(nodes, n, v)

        #---------------------- With dirichlet boundary conditions ------------
        if 'Dir' in bc_type and 'Stationary' in sol_type:
            
            # These two matrices are called the Dirichlet matrices
            # of both the stiffness equivalent and mass equivalent matrices.

            Kd = K[1:n-1, 1:n-1]
            Md = M[1:n-1, 1:n-1]

            print 'K shape is:\n', K.shape

            if 'y'in analysis_param[0] and 'n' in analysis_param[1]:
                n_vals = int(analysis_param[2])
                v = linalg.eigvalsh(Kd, Md, eigvals = (0, n_vals-1))
#                v = v/2
                print 'The Eigenvalues are:\n', v
                return v

            elif 'y'in analysis_param[0] and 'y'in analysis_param[1]:
                n_vals = int(analysis_param[2])
                n_vects = int(analysis_param[3])
                n_solutions = max(n_vals,n_vects)
                v, Dd = linalg.eigh(Kd, Md, eigvals = (0, n_solutions-1))
#                v = v/2
                print 'The Eigenvalues are:\n', v
                E1=zeros(4)
                for i in range(0,4):
                    E1[i]=(i+1)**2*pi**2/((2*pi)**2)
                error=linalg.norm(E1-v)/linalg.norm(E1)
                print 'error', error
                D = zeros((n, n_vects))
                print 'D', D.shape, 'Dd', Dd.shape
                D[1:n-1, :] = Dd
                return v, D

            elif 'n' in analysis_param[0] and 'y' in analysis_param[1]:
                n_vects = int(analysis_param[3])
                v, Dd = linalg.eigh(Kd, Md, eigvals = (0, n_vects-1))
                D = zeros((n, n_vects))
                D[1:n-1, :] = Dd
                return D
            else:
                print 'error: If you dont want anything why do you solve?'

     
        #----------------- With Bloch periodic boundary conditions ------------
        elif 'Bloch' in bc_type:

            import cmath 
            from numpy import asarray
            K = asarray(K, dtype = complex)
            M = asarray(M, dtype = complex)
            print 'K shape is:\n', K.shape


            # Bloch analysis parameters

            xi = nodes[0, 0]   # initial x
            xf = nodes[n-1, 0]  # final x

            n_vals = int(analysis_param[2])  # number of eigenvales to compute

            nk = int(analysis_param[4]) # number of k to sweep
            kmax = 4.*pi/xf
            kmin = -0.0
            k_range = linspace(kmin, kmax, num = nk)
            omega = zeros( (len(k_range), n_vals) )
            E = zeros( (len(k_range), n_vals) )

            print 'Number of eigenvales to compute: ', n_vals, \
                  '\nNumber of wave numbers to sweep: ', nk, \
                  ' in ',  [k_range[0],k_range[nk-1]]

            # Bloch-Periodicity imposition

            ll = 0
            Kaux = K.copy()
            Maux = M.copy()

            for k in k_range:   # Loop over the different wave numbers k
                fi = cmath.exp(1.0j*k*xi)
                ff = cmath.exp(1.0j*k*xf)
                K = Kaux.copy()
                M = Maux.copy()


                for i in range(0, n):
                    K[0, i] = K[0, i]*fi.conjugate()
                    K[i, 0] = K[i, 0]*fi
                    K[n-1, i] = K[n-1, i]*ff.conjugate()
                    K[i, n-1] = K[i, n-1]*ff

                    M[0, i] = M[0, i]*fi.conjugate()
                    M[i, 0] = M[i, 0]*fi
                    M[n-1 , i] = M[n-1, i]*ff.conjugate()
                    M[i, n-1] = M[i, n-1]*ff

                K[n-1, :] = K[0, :] + K[n-1, :]
                K[:, n-1] = K[:, 0] + K[:, n-1]

                M[n-1, :] = M[0, :] + M[n-1, :]
                M[:, n-1] = M[:, 0] + M[:, n-1]

                Kd = K[1:n, 1:n]
                Md = M[1:n, 1:n]

                vals = linalg.eigvalsh(Kd, Md, eigvals = (0, n_vals-1) )

                for i in range(0, n_vals):
                    omega[ll, i] = sqrt( abs(vals[i]) )

                for i in range(0, n_vals):
                    E[ll, i] = vals[i]

                ll = ll + 1
            plt.figure(2)
            plt.hold(True)
            legend = []
            for i in range(0,n_vals):
                plt.plot(k_range, E[:, i])
                legend.append('n = '+str(i+1))

            plt.plot(k_range, (k_range)**2, '--k')
            legend.append('Free Electron')
            plt.title('Dispersion relation')
            plt.legend(legend,loc=2)
            plt.xlabel('Nondimensional wave number - $a\kappa/\pi$')
            plt.ylabel('Nondimensional energy - $2ma^2E/\hslash^2$')
            plt.xlim(xmax=2)
            plt.grid()
            plt.show()
        else:
            print 'Only Dirichlet and Bloch for now, sorry'
    
    #------------------------------  2D problem -------------------------------
    elif dimension == 2:
        
        from matrices import global_stiffness_matrix, global_mass_matrix, \
                             global_potential_matrix
        from vectors import dirichlet_vector, build_solution
        
        h = 1.
        m = 1.        
        print 'Building matrices...\n'
        boundary = read_bc(bc_filename)# Read boundary conditions from .bc file
        #===================== Build stiffness matrix ==================
        stif = global_stiffness_matrix(nodes, triangles)
        stif = h/(2.*m)*stif
            
        #---------------------- With dirichlet boundary conditions ------------
        if 'Dir' in bc_type and 'Stationary' in sol_type:
            
            dirichlet = boundary[0] # Dirichlet conditions are the first list  
            
            
            #===================== Build vectors ===========================
            stif_d, d, remove, g = dirichlet_vector(bc_lines, dirichlet, \
                                                                           stif)    
            # ==================== Build mass matrix =======================
            mass_d = global_mass_matrix(nodes, triangles, remove) 
            
            #===================== Build potential matrix ==================
            v_d = global_potential_matrix(nodes, triangles, v, remove)
#            print "nodes",nodes
#            print "stif",stif
#            print "mass",mass_d
#            print "v_d", v_d
        
            print 'Solving eigenvalue problem...\n'
            if 'y'in analysis_param[0] and 'n' in analysis_param[1]:
                n_vals = int(analysis_param[2])
                v = linalg.eigvalsh(stif_d + v_d, mass_d, \
                                        eigvals = (0, n_vals-1))
#                v = v/2
                print 'The Eigenvalues are:\n', v
                return v
    
            elif 'y'in analysis_param[0] and 'y'in analysis_param[1]:
                n_vals = int(analysis_param[2])
                n_vects = int(analysis_param[3])
                n_solutions = max(n_vals,n_vects)
                v, dir_solution = linalg.eigh(stif_d + v_d, mass_d, \
                                         eigvals = (0, n_solutions-1))
#                v = v/2
    
                solution = zeros((n, n_vects))
                           
                for i in range(n_vects):
                   solution[:,i] = build_solution(dir_solution[:, i], g, remove)
                
                return v, solution
    
            elif 'n'in analysis_param[0] and 'y'in analysis_param[1]:
                n_vects = int(analysis_param[3])
                v, dir_solution = linalg.eigh(stif_d + v_d, mass_d, \
                                        eigvals = (0, n_vects-1))
                
                solution = zeros((n, n_vects))
                           
                for i in range(n_vects):
                   solution[:,i] = build_solution(dir_solution[:, i], g, remove)
                
                return solution
            else:
                print 'error: If you dont want anything why do you solve?'
        
        #---------------------- With Bloch boundary conditions ------------
        elif 'Bloch' in bc_type and 'Stationary' in sol_type:
            from numpy import asarray,vstack
            from vectors import  reference_image_bloch_vectors
            from utils import bloch_multiplication, bloch_sum
            # ==================== Build mass matrix =======================
            mass = global_mass_matrix(nodes, triangles) 
            #===================== Build potential matrix ====================
            v = global_potential_matrix(nodes, triangles, v)
            #================ Convert to complex matrices ====================
            stif = asarray(stif, dtype = complex)
            mass = asarray(mass, dtype = complex) 
            v = asarray(v, dtype = complex)
            #===== Sum kinetic and potential matrices of the Hamiltonian ======
            stif = stif + v
            #= Retrieve the list of image and reference bloch boundary nodes ==
            bloch = boundary[2]
            ref_im = reference_image_bloch_vectors(bc_lines, bloch)
            #============= Discretize the wavenumber dommain ================== 
            nk_x = int(analysis_param[4]) # number of k to sweep in x
            nk_y = int(analysis_param[5]) # number of k to sweep in y
            k_max = pi/float(analysis_param[6])
            k_min = -k_max
            k_range_x = linspace(k_min, k_max, num = nk_x)
            k_range_y = linspace(k_min, k_max, num = nk_y)
            n_vals = int(analysis_param[2])             
            #================ Initiate global variable =======================
            energy = zeros( (nk_x * nk_y, n_vals) )
            
#            print 'Number of eigenvales to compute: ', n_vals, \
#                  '\nNumber of wave numbers to sweep: ', 'in x: ', nk_x, \
#                  ' and in y: ', nk_y, ' in ', \
#                  [k_range_x[0],k_range_x[nk_x-1]], \
#                  [k_range_y[0],k_range_y[nk_y-1]], '\n'
            #======== Defin<xe a new image reference ===========================      
            ref_im_aux = ref_im[0] # image(im) reference(ref) for complex
                                      # multimplication (mul) operations 
            
            for bl in range(1, len(ref_im)):
                ref_im_aux = vstack((ref_im_aux, ref_im[bl]))
                
            #=============== Find who are the corners =======================
            for corner in list(ref_im_aux[:, 0]):
                if list(ref_im_aux[:, 0]).count(corner) == 2:
                    break
            for corner2 in list(ref_im_aux[:, 1]):
                if list(ref_im_aux[:, 1]).count(corner2) == 2:
                    break
            #==== Relate the corners and delete their previous links===========            
            
            im_mul = list(set(list(ref_im_aux[:, 1])))
            ref_mul = []
            for i in im_mul:
                ref_mul.append(list(ref_im_aux[:, 0])[list(ref_im_aux[:, 1]).index(i)])
            ref_im_mul = array([ref_mul, im_mul]).T # Rearranged without 
                                                    # repeated reference 
                                                    # to corner        
#            ref_im_sum = ref_im_mul
#            from numpy import delete
#            for i in range(2):#Delete the links between conrner2 and other nodes
#                ref_im_sum=delete(ref_im_sum, list(ref_im_sum[:, 1]).index(corner2), 0)
#            # Make corner2 the image of corner
#            ref_im_sum[list(ref_im_mul[:, 0]).index(corner), 1] = corner2    
            
            ref_im_mul[list(ref_im_mul[:, 1]).index(corner2), 0] = corner
           
            print "ref_im_mul", ref_im_mul          
            #======================= Main cycles ===============================
            i = 0         
            print 'Calculating each of the ', nk_x * nk_y, \
                   ' points in k plane...\n'
            for k_x in k_range_x:     # Loop over the different wave numbers k in x
                    
                for k_y in k_range_y: # Loop over the different wave numbers k in y
                    # Multiply boundary nodes of each matrix by phase factor                   
                    new_stif, new_mass = bloch_multiplication(k_x, \
                                                    k_y, nodes, \
                                                    ref_im_mul, stif.copy(), \
                                                    mass.copy())
                    new_stif, new_mass = bloch_sum(ref_im_mul, new_stif, new_mass)
                    
                    vals = linalg.eigvalsh(new_stif, new_mass, \
                                           eigvals = (0, n_vals-1))
                    for val in range(0, n_vals):
                        energy[i, val] = vals[val]
                    
                    print i+1, ' points out of :', nk_x * nk_y, '\n'  
                    i = i+1
            from meshUtils import meshtr2D
            k_mesh = meshtr2D(k_min, k_max, k_min, k_max, nk_x, nk_y)
            return k_mesh, energy       
       #------------------ If they ask for something we don't have-----------
        else: 
            print 'Hang on, we are currently working on more options.'
    else:
        print 'Only 1D and 2D for now. Sorry'
    print 'done'

def mat_assembly_1d(nodes, n, v):
    """
        Assembly the equivalent stiffness and equivalent mass matrices for a
        right to left numbered 1D mesh.

        Parameters:
        -----------


        nodes:	Numpy array matrix of nodes containing the coordinates of
                    the nodes from the discretized domain.
                    nodes is an array like matrix of dimension (n,3).

                    Where each column represents the value of the nodes on
                    one of the three coordinate axes x,y,z.

        n:          Number of degree of freedom

        v:	    numpy arraylike vector of dimension (n_lines,1).
                    For Schrödinger case represents the potential acting
                    over the domain.


  	Last modification: date 27/10/2011
    """

    # Initialization of the equivalent stiffness matrix
    K = zeros((n, n))
    # Initialization of the equivalent mass matrix
    M = zeros((n, n))



    # Value for the distance between the first 2 nodes and the last 2 nodes

    Li = abs(nodes[1, 0]-nodes[0, 0])
    Lf = abs(nodes[n-1, 0]-nodes[n-2, 0])


    # Matrices assembly
    K[0, 0] = 1/Li + Li*v[0]/3
    K[n-1, n-1] = 1/Lf + Lf*v[n-2]/3

    for i in range(0, n-1):
        L = abs(nodes[i+1, 0]-nodes[i, 0])       #Distance between nodes
        if i!=0:
            K[i, i] = 2/L + L*(v[i-1]+v[i])/3   # Central diagonal
        K[i, i+1] = -1/L + L*v[i-1]/6          # Upper diagonal
        K[i+1, i] = -1/L + L*v[i-1]/6          # Lower diagonal

    M[0, 0] = Li/3
    M[n-1, n-1] = Lf/3
    for i in range(0, n-1):
        L = abs(nodes[i+1, 0]-nodes[i, 0])      #Distance between nodes
        if i != 0:
            M[i, i] = 2*L/3                    # Central diagonal
        M[i, i+1] = L/6                        # Upper diagonal
        M[i+1, i] = L/6                        # Lower diagonal

    return K, M


        
        











