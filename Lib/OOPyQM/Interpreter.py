# -*- coding: utf-8 -*-
"""
Module interpreter. Will be a place to temorarily hold the class that
receives a simulation file and then assembles the linear equations to solve
"""
__author__ = ['Santiago Echeverri Chacón']

class Interpreter():
    """
    Instances from this class are meant to receive a simulation object and
    construct from it the linear equations that then get to be solved by 
    an appropiate solver.
    """
    def build_QM_dirichlet_eq(self,simulation):

        h = 1.
        m = 1.        
        print 'Building matrices...\n'        
        stif = self.global_stiffness_matrix(simulation)
        stif = h/(2.*m)*stif 
        [stif_d, d, remove, g] = self.dirichlet_vector(simulation, stif)
        mass_d = self.global_mass_matrix(simulation, remove)
        v_d = self.global_potential_matrix(simulation, remove)
        equation = {'left_side': stif_d + v_d, 'right_side':mass_d, \
                    'sol_vec':g, 'dir_positions':remove}
        return equation
        
#    def build_QM_bloch_eq(self, simulation):
#        from numpy import asarray,vstack
#        h = 1.
#        m = 1.        
#        print 'Building matrices...\n'        
#        stif = self.global_stiffness_matrix(simulation)
#        stif = h/(2.*m)*stif 
#        mass = self.global_mass_matrix(simulation)
#        v = self.global_potential_matrix(simulation)
#        #================ Convert to complex matrices ====================
#        stif = asarray(stif, dtype = complex)
#        mass = asarray(mass, dtype = complex) 
#        v = asarray(v, dtype = complex)
#         #===== Sum kinetic and potential matrices of the Hamiltonian ======
#        stif = stif + v
#        bloch = simu.domain.boundaries.bloch      
#        #Will continue...
        
    def global_stiffness_matrix(self, simulation, vectorial = False):
        """
        Funtion for the assembly of the global stiffnes matrix of a 
        simulation file. 
        
        Parameters:
        -----------
        simulation: Instance of class Simulation( ) which contains all 
                    of the information needed for an analysis.
                    Read more about the structure of this class somewhere.
        
        Returns:
        --------
        
        glo_stif:  Matrix like array of size (n_nodes, n_nodes). Represents
                   the global stiffness matrix of the system
        """
        from numpy import zeros
        # Discuss this error with Edward      
        #print isinstance(simulation, Simulation)
        #from Classes import Simulation
        #print isinstance(simulation, Simulation)
        if True:
            nodes_coords = simulation.domain.nodes.coords            
            n_nodes = simulation.domain.nodes.n
            
            sim_elements = simulation.domain.elements
            if simulation.dimension == 1:
                elements = sim_elements.lines
                # Pending to add 1D support
            elif simulation.dimension == 2:
                if 'triangles' in sim_elements.__dict__:
                    elements = sim_elements.triangles
                elif 'quads' in sim_elements.__dict__:
                    elements = sim_elements.quads
                else:
                    print 'Wait untill other elements are supported'
            else:
                print 'No 3 dimensional simulations supperted'
            if vectorial:
                n_elements = elements.n_elements   
                el_set = elements.el_set
                glo_stif = zeros((2*n_nodes, 2*n_nodes))
                for el in range(n_elements):
                    lo_stif = elements.build_local_stiffness(nodes_coords, el)
                    
                    for i in range(1, lo_stif.shape[0]/2+1):
                        for j in range(1, lo_stif.shape[0]/2+1):                            
                            glo_stif[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] = \
                            glo_stif[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] + \
                            lo_stif[2*(i-1),2*(j-1)]
                            glo_stif[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] = \
                            glo_stif[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] + \
                            lo_stif[2*(i-1)+1,2*(j-1)+1]
            else:
                n_elements = elements.n_elements
                el_set = elements.el_set
            
                glo_stif = zeros((n_nodes, n_nodes))            
                
                for el in range(n_elements):
                    lo_stif = elements.build_local_stiffness(nodes_coords, el)
                    for i in range(1, lo_stif.shape[0]+1):
                        for j in range(1, lo_stif.shape[0]+1):
                            
                            glo_stif[el_set[el, i]-1, el_set[el, j]-1] = \
                            glo_stif[el_set[el, i]-1, el_set[el, j]-1] + \
                            lo_stif[i-1,j-1]
                            #                      ^ due to python's 0 base numbering
            return glo_stif
        else:
            print 'Wrong class. input argument should be an instance of simulation.'
            
    def dirichlet_vector(self, simulation, glo_stif, vectorial = False):
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
        d:        Vector associated with siffness matrix over dirichlet nodes.
        remove:   List with the numbers of columns to be removed.
            
        """ 
        from numpy import zeros, dot, delete 
        
        #assert isinstance(simulation, Simulation)        
        assert 'lines' in simulation.domain.elements.__dict__
        bc_lines = simulation.domain.elements.lines
        dirichlet = simulation.domain.boundaries.dirichlet
        n_nodes = glo_stif.shape[0]         # Number (n) of nodes
        n_lines = simulation.domain.elements.n_lines         # Number of lines
        nodes_line = bc_lines.shape[1]
#not needed        n_bc_d = len(dirichlet)  # Number of Dirichlet boundary conditions
        remove = []              # This vector will tell which rows and columns
                                 # should be removed 
       
        g = zeros(n_nodes)
        if vectorial:
            for tag in dirichlet:
                xvalue = dirichlet[tag][0][0]
                yvalue = dirichlet[tag][0][1]
                for ln in range(n_lines):
                    #print bc_lines[ln, 0]
                    if bc_lines[ln, 0] == int(tag):
                        #print bc_lines[ln, 0],'look here'
                        print bc_lines[ln,1:nodes_line]
                        for node in bc_lines[ln,1:nodes_line]:
                            print node-1
                            g[2*(node - 1)] = xvalue
                            g[2*(node - 1)+1] = yvalue
                            if 2*(node-1) not in remove: # makes a list for removing
                                remove.append(2*(node- 1))# lines and columns.
                                remove.append(2*(node- 1)+1)
                            print remove
#                    else:      
                    # so here I found a bug...  if there are no tags defined 
                    # it will inmediatly assume a zero value.
        else:            
            g = zeros(n_nodes)
            for tag in dirichlet:
                value = dirichlet[tag][0][0]
                
                for ln in range(n_lines):
                    
                    if bc_lines[ln, 0] == int(tag):
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
        
    def global_mass_matrix(self, simulation, vectorial = False, *remove):
        """
        funtion for the assembly of the global mass matrix.
        
        Parameters:
        -----------
        simulation: Instance of class Simulation( ) which contains all 
                    of the information needed for an analysis.
                    Read more about the structure of this class somewhere.
                    
        remove:     List with the numbers of columns to be removed
        Returns:
        --------
        glo_mass: Global mass matrix where the dirichlet positions have 
                  been removed.
        """
        from numpy import zeros, delete
        #If simulation is an instance of Classes.Simulation()
        if True:
            nodes_coords = simulation.domain.nodes.coords            
            n_nodes = simulation.domain.nodes.n
            sim_elements = simulation.domain.elements
            if simulation.dimension == 1:
                elements = sim_elements.lines
                # Pending to add 1D support
            elif simulation.dimension == 2:
                if 'triangles' in sim_elements.__dict__:
                    elements = sim_elements.triangles
                elif 'quads' in sim_elements.__dict__:
                    elements = sim_elements.quads
                else:
                    print 'Wait untill other elements are supported'
            else:
                print 'No 3 dimensional simulations supperted'
            n_elements = elements.n_elements   
            el_set = elements.el_set
            if vectorial:
                glo_mass = zeros((2*n_nodes, 2*n_nodes))
                for el in range(n_elements):
                    lo_mass = elements.local_mass_matrix(nodes_coords, el)
                    for i in range(1, lo_mass.shape[0]/2+1):
                        for j in range(1, lo_mass.shape[0]/2+1):                            
                            glo_mass[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] = \
                            glo_mass[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] + \
                            lo_mass[2*(i-1),2*(j-1)]
                            glo_mass[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] = \
                            glo_mass[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] + \
                            lo_mass[2*(i-1)+1,2*(j-1)+1]
            else:
                glo_mass = zeros((n_nodes, n_nodes))# Initiate Global (glo) mass matrix
                
                for el in range(n_elements):
                    #pending for iteration over elements }
                    # call to local_mass_matrix
                    lo_mass = elements.local_mass_matrix(nodes_coords, el)
                    for i in range(1, lo_mass.shape[0]+1):
                        for j in range(1, lo_mass.shape[0]+1):                       
                            glo_mass[el_set[el, i]-1, el_set[el, j]-1] = \
                            glo_mass[el_set[el, i]-1, el_set[el, j]-1] + \
                            lo_mass[i-1,j-1]
            glo_mass_d = delete(glo_mass, remove, 0)
            glo_mass_d = delete(glo_mass_d, remove, 1)
            return glo_mass_d
                 
        else:
            print 'Wrong class. input argument should be an instance of simulation.'
        
        
    def global_potential_matrix(self, simulation, *remove):
        """
            funtion for the ensamble of the global potential matrix
            
            Parameters:
            -----------
            simulation: Instance of class Simulation( ) which contains all 
                    of the information needed for an analysis.
                    Read more about the structure of this class somewhere.
                    
            remove:     List with the numbers of columns to be removed
                        
            Returns:
            --------
        """
        from numpy import zeros, delete
        #If simulation is an instance of Classes.Simulation()
        if True:
            nodes_coords = simulation.domain.nodes.coords            
            n_nodes = simulation.domain.nodes.n
            sim_elements = simulation.domain.elements
            v = simulation.body_parameter
            if simulation.dimension == 1:
                elements = sim_elements.lines
                # Pending to add 1D support
            elif simulation.dimension == 2:
                if 'triangles' in sim_elements.__dict__:
                    elements = sim_elements.triangles
                elif 'quads' in sim_elements.__dict__:
                    elements = sim_elements.quads
                else:
                    print 'Wait untill other elements are supported'
            else:
                print 'No 3 dimensional simulations supperted'
            n_elements = elements.n_elements
            el_set = elements.el_set
            glo_v = zeros((n_nodes, n_nodes))# Initiate Global (glo) mass matrix
            
            for el in range(n_elements):
                #pending for iteration over elements }
                # call to local_mass_matrix
                lo_v = elements.local_potential_matrix(nodes_coords, el_set, v, el)
                for i in range(1, lo_v.shape[0]+1):
                    for j in range(1, lo_v.shape[0]+1):                       
                        glo_v[el_set[el, i]-1, el_set[el, j]-1] = \
                        glo_v[el_set[el, i]-1, el_set[el, j]-1] + \
                        lo_v[i-1,j-1]
            glo_v_d = delete(glo_v, remove, 0)
            glo_v_d = delete(glo_v_d, remove, 1)
            return glo_v_d
             
        else:
            print 'Wrong class. input argument should be an instance of simulation.'
    def reference_image_bloch_vectors(simulation, bc_lines, bloch):
        """
        This function loads the lines of the boundary and the list that contains 
        bloch periodicity conditions, and builds a list of arrays that relate
        image nodes with their corresponding reference nodes.
        
        Paprameters:
        -----------
        
        bc_lines:   Array of line elements.
    
        bloch:      List of bloch conditions as read by a '.bc' file.
    
        Returns:
        --------
        
        ref_im:     Each array in the list 'ref_im', has in it's first column the 
                    reference node and on it's second collumn the image node for 
                    that particular reference node.
        """
        from numpy import zeros, copy
        bc_lines = simulation.domain
        n_bloch = bloch.shape[0]
        it_bloch = range(n_bloch)    
        n_lines = bc_lines.shape[0]
        print 'n_lines', n_lines
        #=========== Create a list of arrays ======================================
        ref_im = []    
        for i in it_bloch:    
            ref_im.append(zeros((n_lines/(2*n_bloch)+1,2), dtype = int))  
        #=============== Loop and fill the lists ================================== 
        count = 0
        i = 0
        j = 0
        
        while count == 0:
            
            for bl in it_bloch:
                
                
                if bc_lines[j, 0] == bloch[bl, 0]:
                    ref_im[bl][i, 0] = bc_lines[j, 1]
                    if i == n_lines/(2*n_bloch)-1:
                        ref_im[bl][i+1, 0] = bc_lines[j, 2]
                    
                elif bc_lines[j, 0] == bloch[bl, 1]:
                    ref_im[bl][i, 1] = bc_lines[j, 1]
                    if i == n_lines/(2*n_bloch)-1:
                        ref_im[bl][i+1, 1] = bc_lines[j, 2]
            i = i+1
            j = j+1
    
    	
            if i == n_lines/(2*n_bloch):
                i = 0
            if j >= n_lines:
                count = 1
        
        for bl in it_bloch:
            if bloch[bl, 2]*bloch[bl, 3] == -1:
                ref_im[bl][:,1] = copy(ref_im[bl][::-1,1])
        
        
        return ref_im