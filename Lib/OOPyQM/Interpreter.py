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
    def build_static_EM_eq(self, simulation):
        epsilon = 1
        mu = 1
        simulation.domain.read_bc_file(simulation.bc_filename)
        print 'Building matrices...\n'    
        stif = self.global_stiffness_matrix(simulation, True)
        stif = (1/mu) * stif
        [stif_d, d, remove, g] = self.dirichlet_vector(simulation, stif, True)
        n_nodes = stif_d.shape[0]
        q = self.newman_vector(simulation, remove, n_nodes, vectorial = True)
        equation = {'left_side': stif_d, 'right_side':-d-q, \
                    'sol_vec':g, 'dir_positions':remove, 'vectorial':True}
        return equation
    def build_harmonic_EM_eq(self, simulation):
        epsilon = 1
        mu = 1
        simulation.domain.read_bc_file(simulation.bc_filename)
        print 'Building matrices...\n'    
        stif = self.global_stiffness_matrix(simulation, True)
        stif = (1.0/mu) * stif
        [stif_d, d, remove, g] = self.dirichlet_vector(simulation, stif, True)
        mass_d = self.global_mass_matrix(simulation, True, remove)
        mass_d = epsilon * mass_d
        equation = {'left_side': stif_d, 'right_side': mass_d, \
                    'sol_vec':g, 'dir_positions':remove, 'vectorial':True}
        return equation

    def build_EM_bloch_eq(self, simulation):
        """
        Form the equation to solve for a Bloch boundary simulation.
        Bloch boundaries must be matched and have the same number of nodes.
        The user must state in the .bc file the image and references 
        physical entities and interpret with a - sign if they are mirrored.
        """
        from numpy import array, asarray, vstack
        mu = 1.
        epsilon = 1.        
        print 'Building matrices...\n'       
        simulation.domain.read_bc_file(simulation.bc_filename)
        
        
        stif = self.global_stiffness_matrix(simulation, True)
        stif = (1.0/mu) * stif
        mass = self.global_mass_matrix(simulation)
        #================ Convert to complex matrices ====================
        stif = asarray(stif, dtype = complex)
        mass = epsilon * asarray(mass, dtype = complex)   
        ref_im = self.reference_image_bloch_vectors(simulation)
        
        #======== Define a new image reference ===========================      
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
        ref_im_mul[list(ref_im_mul[:, 1]).index(corner2), 0] = corner
        print "ref_im_mul", ref_im_mul        
        
        equation = {'left_side': stif, 'right_side': mass, \
                    'ref_im_mul':ref_im_mul, 'vectorial':True}
        return equation
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
                elements = sim_elements.lines.el_set
                # Pending to add 1D support
            elif simulation.dimension == 2:
                if 'triangles' in sim_elements.__dict__:
                    elements = sim_elements.triangles
                elif 'quads' in sim_elements.__dict__:
                    elements = sim_elements.quads
                else:
                    raise NotImplementedError('Wait untill other elements\
                                               are supported')
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
            raise  TypeError('Wrong class. input argument should be\
                              an instance of simulation.')
            
    def dirichlet_vector(self, simulation, glo_stif, vectorial = False):
        """
        This function computes the d vector associated to Dirichlet
        boundary conditions on certain nodes of the domain. 
        
        Parameters:
        -----------
        simulation: Instance of Simulation() which contains many of the 
                    characteristics of the problem definition.
                   
        glo_stiff: matrix like array of size (n_nodes,n_nodes) called
                   the stiffness matrix of the system.
        vectorial: Indicates if the simulation has to be defined as scalar or
                   vectorial. 
           
        Returns:
        --------
        glo_stif_d:  global stiffness matrix reduced by removing the rows and 
                     columns asociated with the dirichlet conditions.
        d:        Vector associated with siffness matrix over dirichlet nodes.
        remove:   List with the numbers of columns to be removed.
        g:                  
        """ 
        from numpy import zeros, dot, delete 
        
        #assert isinstance(simulation, Simulation)        
        assert 'lines' in simulation.domain.elements.__dict__
        bc_lines = simulation.domain.elements.lines.el_set
        dirichlet = simulation.domain.boundaries.dirichlet
        n_nodes = glo_stif.shape[0] 
        n_lines = simulation.domain.elements.lines.n_elements        
        nodes_line = bc_lines.shape[1]
        remove = []              
        
                                 # should be removed 
       
        g = zeros(n_nodes)
        if vectorial:
            for tag in dirichlet:
                for ln in range(n_lines):
                    #print bc_lines[ln, 0]
                    if bc_lines[ln, 0] == int(tag):
                        #print bc_lines[ln, 0],'look here'
                        for node in bc_lines[ln,1:nodes_line]:
                            xvalue = dirichlet[tag][0][0]
                            yvalue = dirichlet[tag][0][1]
                            if type(xvalue) == str:
                                from math import sqrt 
                                x = simulation.domain.nodes.coords[node, 0]
                                y = simulation.domain.nodes.coords[node, 1]
                                xvalue = eval(xvalue)
                                yvalue = eval(yvalue)
                            g[2*(node - 1)] = xvalue
                            g[2*(node - 1)+1] = yvalue
                            if 2*(node-1) not in remove: # makes a list for removing
                                remove.append(2*(node- 1))# lines and columns.
                                remove.append(2*(node- 1)+1)
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
    def newman_vector(self, simulation, remove, n_nodes, vectorial = False):      
        """
        This function computes the q vector associated to Newman
        boundary conditions on certain nodes of the domain. 
        
        Parameters:
        -----------
        simulation:  Instance of class Simulation( ) which contains all 
                    of the information needed for an analysis.
                    Read more about the structure of this class somewhere.
                            
        remove:   List with the numbers of columns to be removed.
        
        Returns:
        --------
        q:      Newman vector defined by... either static values or a 
                function defined by an expresion on the bc file.
        """ 
        from numpy import zeros, dot, delete 
        #assert isinstance(simulation, Simulation)        
        assert 'lines' in simulation.domain.elements.__dict__
        bc_lines = simulation.domain.elements.lines.el_set
        lo_newman = simulation.domain.elements.lines.local_newman
        newman = simulation.domain.boundaries.newman
        n_lines = simulation.domain.elements.lines.n_elements         # Number of lines
        nodes_line = bc_lines.shape[1]
        nodes = simulation.domain.nodes.coords
#not needed        n_bc_d = len(dirichlet)  # Number of Dirichlet boundary conditions
        q = zeros(n_nodes)
        if vectorial:
            
            for tag in newman:
                for ln in range(n_lines):
                    if bc_lines[ln, 0] == int(tag):
                        lo_q = lo_newman(q, nodes, newman, tag, ln, vectorial)    
                        count = 0                    
                        for node in bc_lines[ln,1:nodes_line]:
                            q[2*(node - 1)] = lo_q[count]
                            q[2*(node - 1)+1] = lo_q[count]
                            count += 1
                  
        else:    
            for tag in newman:
                for ln in range(n_lines):
                    if bc_lines[ln, 0] == int(tag):
                        lo_q = lo_newman(q, nodes, newman, tag, ln)
                        for node in bc_lines[ln,1:nodes_line]:
                            q[bc_lines[ln, 1]-1] = lo_q[count]                        
        return q
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
                elements = sim_elements.lines.el_set
                # Pending to add 1D support
            elif simulation.dimension == 2:
                if 'triangles' in sim_elements.__dict__:
                    elements = sim_elements.triangles
                elif 'quads' in sim_elements.__dict__:
                    elements = sim_elements.quads
                else:
                    raise NotImplementedError('Wait untill other elements\
                                               are supported')
            else:
                raise NotImplementedError('No 3 dimensional simulations\
                                           supperted')
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
            raise TypeError('Wrong class. input argument should be an\
                             instance of simulation.')
        
        
#    def global_potential_matrix(self, simulation,vectorial = False, *remove):
#        """
#            funtion for the ensamble of the global potential matrix
#
#            This function is pending to develop because the electromagnetism 
#            related proyect I'm using, which doesn't necesarily calls 
#            potential matrices. 
#                        
#            Parameters:
#            -----------
#            simulation: Instance of class Simulation( ) which contains all 
#                    of the information needed for an analysis.
#                    Read more about the structure of this class somewhere.
#                    
#            remove:     List with the numbers of columns to be removed
#                        
#            Returns:
#            --------
#        """
#        from numpy import zeros, delete
#        #If simulation is an instance of Classes.Simulation()
#        if True:
#            v = simulation. 
#            nodes_coords = simulation.domain.nodes.coords            
#            n_nodes = simulation.domain.nodes.n
#            sim_elements = simulation.domain.elements
#            v = simulation.body_parameter
#            if simulation.dimension == 1:
#                elements = sim_elements.lines
#                # Pending to add 1D support
#            elif simulation.dimension == 2:
#                if 'triangles' in sim_elements.__dict__:
#                    elements = sim_elements.triangles
#                elif 'quads' in sim_elements.__dict__:
#                    elements = sim_elements.quads
#                else:
#                    print 'Wait untill other elements are supported'
#            else:
#                print 'No 3 dimensional simulations supperted'
#            n_elements = elements.n_elements
#            el_set = elements.el_setf
#            if vectorial:
#                glo_v = zeros((2*n_nodes, 2*n_nodes))
#                for el in range(n_elements):
#                    lo_v = elements.local_potential_matrix(nodes_coords, el_set, v el)
#                    for i in range(1, lo_v.shape[0]/2+1):
#                        for j in range(1, lo_v.shape[0]/2+1):                            
#                            glo_v[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] = \
#                            glo_v[2*(el_set[el, i]-1), 2*(el_set[el, j]-1)] + \
#                            lo_v[2*(i-1),2*(j-1)]
#                            glo_v[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] = \
#                            glo_v[2*(el_set[el, i]-1)+1, 2*(el_set[el, j]-1)+1] + \
#                            lo_v[2*(i-1)+1,2*(j-1)+1]
#            else:
#                glo_v = zeros((n_nodes, n_nodes))# Initiate Global (glo) mass matrix
#                
#                for el in range(n_elements):
#                    #pending for iteration over elements }
#                    # call to local_mass_matrix
#                    lo_v = elements.local_potential_matrix(nodes_coords, el_set, v, el)
#                    for i in range(1, lo_v.shape[0]+1):
#                        for j in range(1, lo_v.shape[0]+1):                       
#                            glo_v[el_set[el, i]-1, el_set[el, j]-1] = \
#                            glo_v[el_set[el, i]-1, el_set[el, j]-1] + \
#                            lo_v[i-1,j-1]
#                glo_v_d = delete(glo_v, remove, 0)
#                glo_v_d = delete(glo_v_d, remove, 1)
#                return glo_v_d
#                 
#        else:
#            print 'Wrong class. input argument should be an instance of simulation.'
    def reference_image_bloch_vectors(self, simulation):
        """
        This function loads the lines of the boundary and the list that contains 
        bloch periodicity conditions. It builds a list of arrays that relate
        image nodes with their corresponding reference nodes.
        
        Paprameters:
        -----------
        
        bc_lines:   Array of line elements.
    
        bloch:      List of bloch conditions as read by a '.bc' file.
    
        Returns:
        --------
        
        ref_im:     Each array in the list 'ref_im', has in it's first column the 
                    reference node and on it's second column the image node for 
                    that particular reference node.
        """
        from numpy import array, zeros, copy, vstack
        bloch = simulation.domain.boundaries.bloch
        bc_lines = simulation.domain.elements.lines.el_set
        n_lines = simulation.domain.elements.lines.n_elements
        order = simulation.domain.elements.lines.order
        vectorial = simulation.domain.elements.lines.vectorial
        
        n_bloch = bloch.shape[0]
        it_bloch = range(n_bloch)    
       
        lines_dict = {}
        for ln in bc_lines:
            if str(ln[0]) not in lines_dict:
                lines_dict[str(ln[0])] = array(ln[1:])
            else:
                lines_dict[str(ln[0])] = vstack((lines_dict[str(ln[0])],array(ln[1:])))
       
        ref_im = [] 
        for bl in bloch:
            try:
                lines_dict[str(bl[0])].shape == lines_dict[str(bl[1])].shape
            except:
                raise ValueError, 'Matching boundaries must have the same \
                                    number of elements'
            if order == 1:
                if vectorial:
                    ref_im.append(zeros((2*lines_dict[str(bl[0])].shape[0]+1,2), dtype = int))
                else:
                    ref_im.append(zeros((lines_dict[str(bl[0])].shape[0]+1,2), dtype = int))
            elif order ==2:
                if vectorial:
                    ref_im.append(zeros((2*2*lines_dict[str(bl[0])].shape[0]+2,2), dtype = int))
                else:
                    ref_im.append(zeros((2*lines_dict[str(bl[0])].shape[0]+1,2), dtype = int))
            else:
                raise NotImplementedError, 'No higher order support yet'

        tags = lines_dict.keys()
        ref_im_dicts = {}
        i=0
        for bl in bloch:
            ref_im_dict = {}            
            for tag in tags:
                if str(bl[0]) == tag:
                    for ln in range(lines_dict[tag].shape[0]):
                        for j in range(3):
                            tag1 = lines_dict[tag][ln, j]
                            tag2 = lines_dict[str(bl[1])][ln, j]
                            if tag1 and tag2 not in ref_im_dict.keys():
                                ref_im_dict[tag1] = tag2
                    tags.pop(tags.index(tag))   
            ref_im_dicts[i] = ref_im_dict
            i += 1
      
       
        if vectorial:
            for bl in range(len(ref_im)):
                i = 0
                keys = ref_im_dicts[bl].keys()
                for key in keys:
                    ref_im[bl][2*i,:] = [2*key, 2*ref_im_dicts[bl][key]]
                    ref_im[bl][2*i+1,:] = [2*key+1, 2*ref_im_dicts[bl][key]+1]
                    i += 1
        else:
            for bl in range(len(ref_im)):
                keys = ref_im_dicts[bl].keys()
                i = 0
                for key in keys:
                    ref_im[bl][i] = [key, ref_im_dicts[bl][key]]
                    i += 1
        print ref_im
        for bl in it_bloch:
            if bloch[bl, 2]*bloch[bl, 3] == -1:
                ref_im[bl][:,1] = copy(ref_im[bl][::-1,1])
        
        return ref_im