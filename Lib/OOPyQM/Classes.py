# -*- coding: utf-8 -*-
"""
Module Classes will temporalry hold the classes used by PeYeQM to
define instances.
"""
__author__ = ['Santiago Echeverri Chacón']

class Simulation():
    """
    Instances from simulation class will act as a container for the overal 
    problem definition of the simulation.
    It shall include information about:
        * The Domain and its discretization.
        * Options and configuration.
        *Type of simulations and solvers to be used.
    """
    def __init__(self):
        self.domain = Domain()
    sim_type = 'QM'
    time_dependency = 'Stationary'
    dimension = 2
    def __str__(self):
        return 'This is an instance of the Simulation Class'
        
    def read_solver_input(self, filename):
        """
        Reads the solver input from a file of gmsh ASCII format V2.2.
        This input is to be read by the Solver module.
        Parameters:
        -----------
        filename: String which contain filename and path of the output file.
        Returns:
        --------
        solver_input: List with the following parameters:
    
        [dimension,bc_type,parameter,eq,sol_type,analysis_param]
    
        dimension:  int parameter that tells the program wether to solve for a
                1D problem or a 2D problem (not supported yet)
    
        bc_type:     String parameter for the selection of a border condition
                    that can be either:
    
                        'Dir'   For the Dirichlet border condition
                                (Infinite potential well).
    
                        'Bloch' For the periodic formulation of the problem.
                                (Electron in a periodic material )
                
        body_parameter:  
                Is an array that describes the potential actuating over the
                the elements of the domain given by Elems. For each element in
                Elems there is an associated potential value on the same
                position in the array parameter.

                The potential in Scroedinger equation defines the specific
                nature of the problem to be solved. For more details on how
                to define a potential and what does it mean pleas read the
                documentation of the Potential1D function in the module PrePro.
    
       
        sim_type: Can be 'QM' from Quantum Mechanics, or 'EM' for 
                  electromagnetism.
        
    
        sol_type:       String that tells wether to solve the stationary version of
                    the equation or another not yet suported.
    
                    'Stationary'   
    
        solver_param:   Array that contains the information regarding the number
                         of solutions to be computed and wether to save the values
                         or not.
    
                        analysis_param[0]:  String  answer to the question
                                                   save  Eigen Values?
                        analysis_param[1]:  String  answer to the question
                                                   save  Eigen Vectors?
                        analysis_param[2]:  Integer  number of Eigen Values to save
                        analysis_param[3]:  Integer  number of Eigen Vectors to save 
                        analysis_param[4]:  Integer number of wave numbers 
                                            in x to sweep
                        analysis_param[5]:  Integer number of wave numbers 
                                            in y to sweep
                        analysis_param[6]:  biggest value of k. it may be the lenght
                                            of the dommain                                        
                        
        bc_filename:    string that tells where to look for the boundary 
                        conditions  
                            
         Last modification: date 12/02/2013
        """
        from read_mesh import read_solver_input
        solver_input = read_solver_input(filename) 
        if solver_input != []:
            self.dimension = solver_input[0]
            self.bc_type = solver_input[1]
            self.body_parameter = solver_input[2]                # Reasign
            self.sim_type = solver_input[3]
            self.time_dependency = solver_input[4]
            self.solver_param = solver_input[5]         # the solution 
            self.bc_filename = solver_input[6]
        else:
            raise NotImplementedError('Pending implementation')
    

class Domain():
    """ 
    Instances from Domain Class contain the attributes that define 
    characteristics of the simulation such as regions, meshing, and 
    boundary conditions.
    To be included here are instances of other clases such as:
        * Nodes.
        * Elements 
        * And boundaries.
    """
    def __init_(self):
        self.nodes = 0
        self.elements = 0
        self.regions = {}
    def read_mesh_file(self, filename, vectorial = False):
        from read_mesh import read_mesh
        
        self.nodes = Nodes()
        self.elements = Elements()
  
        _nodes, _elements = read_mesh(filename)        
        
        self.nodes.add(_nodes)
        self.elements.add(_elements, vectorial)
    def read_regions_file(self, filename):
        from numpy import delete
        from copy import deepcopy
        import pickle
        print filename
        f = open(filename +'.reg', 'r')
        regions = pickle.load(f)
        for region in regions:
            region.elements = deepcopy(self.elements.all)
            tag = int(region.tag)
            for class_iter in region.elements:
                row = 0
                remove =[]
                for i in region.elements[class_iter].el_set:
                    if i[0] != tag:
                        remove.append(row)
                    row += 1
                el_set = delete(region.elements[class_iter].el_set, remove, 0)
                region.elements[class_iter].el_set = el_set
                region.elements[class_iter].n_elements = el_set.shape[0]
        self.regions = regions
            
    def read_bc_file(self, filename):
        from read_mesh import read_bc
        self.boundaries = Boundaries()
        _bc = read_bc(filename)
        self.boundaries.add(_bc)                
                  
class Nodes():
    """
    Nodes, will be a instance fo handling nodes. Right now its only 
    two attributes are the total number of nodes and the array of 
    ccordinates.
    """
    # def __init__(self):
        
    def add(self, _nodes):        
        from numpy import shape
        self.n = shape(_nodes)[0]                 #Number of nodes
        self.coords = _nodes        
    
class Region():
    """
    Reperesents a region of the Domain. It is thought as a way to store 
    attributes that are specific to certain regions such as the kind of 
    elements that compose it, it's distinctive tag, and material properties
    """
    def __init__(self, tag = '', name = '',material_prop = {}, elements = {}):
        self.tag = tag
        self.name = name
        self.material_prop =  material_prop
        self.elements = elements
    def __str__(self):
        return "Region with tag %s is called %s and has the following"\
                "material properties: \n %s and sets of elements \n %s"\
                %(self.tag, self.name, self.material_prop,self.elements.keys())
    
class Elements():
    """
    Elements for now counts the number f line and triangle elements
    in a format imported from a file by using function read_mesh().
    
    The Idea is that the instance elements contains all the attributes 
    associated with elements
    """
    #Lines = 0    
    def __init__(self):
        self.all = {}
    def add(self, _elements, vectorial= False):
       # from numpy import shape
        if 'lin_Lines' in _elements:
            self.lines = Lines(_elements['lin_Lines'], vectorial)
            self.all['lin_Lines'] = self.lines
        elif 'cuad_Lines' in _elements:
            self.lines = Lines(_elements['cuad_Lines'], vectorial)
            self.all['cuad_Lines'] = self.lines
        else: 
            raise NotImplementedError("No higher order line elements of more"\
                                        "than second order. Maybe you did't"\
                                        "define physical lines?")
        if 'Triangles' in _elements:
            self.triangles = Triangles(_elements['Triangles'])
            self.all['Triangles'] = self.triangles
        if 'cuad_Quads' in _elements:
            self.quads = Quadrilaterals(_elements['cuad_Quads'], vectorial)
            self.all['cuad_Quads'] = self.quads
            
    #It doesn't look so hard to implement square elements

    # It would be cool to have a function and or attributes where 
    # The order of elements is taken into account. 
    #Things such as an element properties like:
    #   * Order
    #   * Area, or other shape attributes like jacobian
    #   * Local matrices...
class Lines():
    """ Line elements are used for solving 1D problems and handling 
        boundary conditions. The instance lines has been made in order
        to contain properties that are defined exlusively for line elements
        such as stiffness matrices for 1D problems and
        interpolation methods for newman and source vectors.
    """    
    def __init__(self, raw_lines, vectorial = False):
        from numpy import shape
        self.el_set = raw_lines
        self.n_elements = shape(self.el_set)[0]
        if shape(self.el_set)[1] == 3:
            self.order = 1            
            self.h = [lambda r: 1.0/2.0*(1-r),\
                      lambda r: 1.0/2.0*(1+r)]
        elif shape(self.el_set)[1] == 4:
            self.order = 2            
            self.h = [lambda r: 1.0/2.0*(1-r)-1.0/2.0*(1-r**2),\
                      lambda r: 1-r**2,\
                      lambda r: 1.0/2.0*(1+r)-1.0/2.0*(1-r**2)]
        else: 
            raise NotImplementedError('Up to second order only')
        self.vectorial = vectorial
    def extract_el_points(self, _nodes, el_id):
        """
        Extracts the coordinates of the nodes in each line element and
        retrieves an array of node coordinates.
        
        Parameters:
        -----------
        nodes:  numpy array of dimension (n_nodes, 3), where n_nodes is the 
                number of nodes forming the mesh, and the three columns 
                represent coordinates (x, y, z).
                    
        el_id:   integer value of the current element
        
        Returns:
        --------
        node_coor: Array with coordinates of coordinate node. 
                   If the element is 3 node then array is 2*3
                   If it is a 2 node element then 2*2.
        """
        from numpy import zeros
        if 'el_set' in self.__dict__:
            if self.order == 2:
                p2_lines = self.el_set   
                node_coords = zeros((2,3))
                for i in range(1,4):
                    node_coords[0,i-1] = _nodes[p2_lines[el_id, i] - 1,0]
                    node_coords[1,i-1] = _nodes[p2_lines[el_id, i] - 1,1]
                   
                return node_coords   
            elif self.order == 1:
                p1_lines = self.el_set                            
                node_coords = zeros((2,2))
                for i in range(1,3):
                    node_coords[0,i- 1] = _nodes[p1_lines[el_id, i] - 1,0]
                    node_coords[1,i- 1] = _nodes[p1_lines[el_id, i] - 1,1]
                return node_coords   
            else:
                raise NotImplementedError('No higher order elements yet, \
                                           there must be an error')
        else:
            raise IOError('No elements parsed, do something else')
    def numeric_J(self, node_coords):
        from numpy.linalg import norm     
        if self.order == 2:             
            L1 = norm(node_coords[:,2] - node_coords[:,0])
            L2 = norm(node_coords[:,1] - node_coords[:,2])
            J = L1 + L2
        elif self.order == 1:
            J = norm(node_coords[:,1] - node_coords[:,0])
        else:
            raise NotImplementedError('No higher order elements yet, \
                                           there must be an error')
        return J
    def local_newman(self, q, nodes, newman, tag, ln, vectorial = False):
        """ Returns the nodal values for a certain newman boundary line 
            element by using the specific definition from the.bc file.
        """
        bc_lines = self.el_set
        nodes_line = bc_lines.shape[1]
        if self.order == 2:
            from scipy.special.orthogonal import p_roots
            from numpy import zeros
            h3 = self.h
             # Definition of integration degree for each of the coordinates
            deg_r = 3
            r, weights_r  = p_roots(deg_r)
            if self.vectorial:
                lo_q = zeros(6)
                counter = 0
                for node in bc_lines[ln,1:nodes_line]:
                    xvalue = newman[tag][0][0]
                    yvalue = newman[tag][0][1]
                    if type(xvalue) == str:
                        from math import sqrt 
                        x = nodes[node, 0]
                        y = nodes[node, 1]
                        xvalue = eval(xvalue)
                        yvalue = eval(yvalue)
                    lo_q[2*counter] = xvalue
                    lo_q[2*counter+1] = yvalue
                    counter += 1
                node_coords = self.extract_el_points(nodes, ln)   
                det_J = self.numeric_J(node_coords)
                for i in range(deg_r):
                    H3 = zeros(6)
                    for k in range(3):
                        H3[2*k] = h3[k](r[i])
                        H3[2*k+1] = h3[k](r[i])
                    lo_q += weights_r[i]*det_J*H3*lo_q
                return lo_q
            else:
                if bc_lines[ln, 0] == int(tag):
                    lo_q = zeros(3)
                    counter = 0
                    for node in bc_lines[ln,1:nodes_line]:
                        value = newman[tag][0][0]
                        if type(xvalue) == str:
                            from math import sqrt 
                            x = nodes[node, 0]
                            value = eval(value)
                        lo_q[counter] = value
                        counter += 1
                    node_coords = self.extract_el_points(nodes, ln)   
                    det_J = self.numeric_J(node_coords)
                    for i in range(deg_r):
                        H3 = zeros(6)
                        for k in range(3):
                            H3[k] = h3[k](r[i])
                        lo_q += weights_r[i]*det_J*H3*lo_q
                return lo_q
                
        
class Quadrilaterals():
    """ Quadrilateral elements give better accuracy and resistance to 
        locking than triangular elements. They can also be more economic
        when meshing compared to their equivalent triangular peers.
        Option vectorial will initiate base functions and nodal 
        representation that accounts for components of the solution function.
        This means that there will be x and y components of the solution.
    """
    
    def __init__(self,raw_quads, vectorial = False):
        from numpy import shape
        self.el_set = raw_quads
        self.n_elements = shape(self.el_set)[0]
        #Definition of second shape functions according to order and type
        if shape(self.el_set)[1] == 5:
            self.order = 1
            self.h =[lambda r, s: 1.0/4.0*(1+r)*(1+s),\
                 lambda r, s: 1.0/4.0*(1-r)*(1+s),\
                 lambda r, s: 1.0/4.0*(1-r)*(1-s),\
                 lambda r, s: 1.0/4.0*(1+r)*(1-s)] 
        elif shape(self.el_set)[1] == 9:
            self.order = 2
            self.h =[lambda r, s: 1.0/4.0*(1+r)*(1+s)-1.0/4.0*(1-r**2)*(1+s)\
                                                 -1.0/4.0*(1-s**2)*(1+r),\
                 lambda r, s: 1.0/4.0*(1-r)*(1+s)-1.0/4.0*(1-r**2)*(1+s)\
                                                 -1.0/4.0*(1-s**2)*(1-r),\
                 lambda r, s: 1.0/4.0*(1-r)*(1-s)-1.0/4.0*(1-r**2)*(1-s)\
                                                 -1.0/4.0*(1-s**2)*(1-r),\
                 lambda r, s: 1.0/4.0*(1+r)*(1-s)-1.0/4.0*(1-r**2)*(1-s)\
                                                 -1.0/4.0*(1-s**2)*(1+r),\
                 lambda r, s: 1.0/2.0*(1-r**2)*(1+s),\
                 lambda r, s: 1.0/2.0*(1-s**2)*(1-r),\
                 lambda r, s: 1.0/2.0*(1-r**2)*(1-s),\
                 lambda r, s: 1.0/2.0*(1-s**2)*(1+r)] 
        else:
            raise NotImplementedError('Not a defined order...')
        if vectorial:
            self.vectorial = vectorial
            
    def extract_el_points(self, _nodes, el_id):
        """
        Extracts the coordinates of the nodes in each quad element and
        retrieves an array of node coordinates.
        
        Parameters:
        -----------
        nodes:  numpy array of dimension (n_nodes, 3), where n_nodes is the 
                number of nodes forming the mesh, and the three columns 
                represent coordinates (x, y, z).
                    
        el_id:   integer value of the current element
        
        Returns:
        --------
        node_coor: Array with coordinates of coordinate node. 
                   If the element is 8 node then array is 2*8
                   If it is a 4 node element then 2*4.
        """
        from numpy import zeros
        if 'el_set' in self.__dict__:
            if self.order == 2:
                p2_quads = self.el_set   
                node_coords = zeros((2,8))
                for i in range(1,9):
                    node_coords[0,i-1] = _nodes[p2_quads[el_id, i] - 1,0]
                    node_coords[1,i-1] = _nodes[p2_quads[el_id, i] - 1,1]
                   
                return node_coords   
            elif self.order == 1:
                p1_quads = self.el_set                            
                node_coords = zeros((2,4))
                for i in range(1,5):
                    node_coords[0,i- 1] = _nodes[p1_quads[el_id, i] - 1,0]
                    node_coords[1,i- 1] = _nodes[p1_quads[el_id, i] - 1,1]
                return node_coords   
            else:
                raise NotImplementedError('No higher order elements yet, \
                                           there must be an error')
        else:
            raise IOError('No elements parsed, do something else')
    
    def numeric_J(self, node_coords, r,s, dHdrs = False):
        """ I should document more often, now I don't get how this\
            function work.
            Parameters:
                node_coords:    Array that contains coordinates of nodes
                r, s:   
        """
        from numpy import dot, zeros
        from numpy.linalg import det
        if self.order == 2:
            J_mat = zeros((2,2))
            dhdr = [1.0/4.0*(s**2 + s +2*r*(s+ 1)), \
                    1.0/4.0*(-s**2 - s +2*r*(s+ 1)), \
                    1.0/4.0*(-s**2 + s +2*r*(-s+ 1)), \
                    1.0/4.0*(s**2 - s +2*r*(-s+ 1)), \
                    -r*(s+ 1), \
                    -1.0/2.0*(-s**2+ 1), \
                    -r*(-s+ 1),\
                    1.0/2.0*(-s**2+ 1)]
            
            dhds = [1.0/4.0*(r**2 + r +2*s*(r+ 1)), \
                    1.0/4.0*(r**2 - r +2*s*(-r+ 1)), \
                    1.0/4.0*(-r**2 + r +2*s*(-r+ 1)), \
                    1.0/4.0*(-r**2 - r +2*s*(r+ 1)), \
                    1.0/2.0*(-r**2 + 1), \
                    -s*(-r+ 1), \
                    -1.0/2.0*(-r**2+ 1), \
                    -s*(r+ 1)]
        else:
            raise NotImplementedError('linear quad functions pending')
    
        J_mat[0, 0] = dot(node_coords[0], dhdr)
        J_mat[1, 0] = dot(node_coords[1], dhdr)
        J_mat[0, 1] = dot(node_coords[0], dhds)
        J_mat[1, 1] = dot(node_coords[1], dhds)
        det_J = det(J_mat)
        assert det_J > 0, "You did something wrong with numbering, because this should be positive"
        inv_J = zeros((2,2))
        inv_J[0, 0] = J_mat[1,1]
        inv_J[1, 1] = J_mat[0, 0]
        inv_J[1, 0] = -J_mat[0, 1]
        inv_J[0, 1] = -J_mat[1, 0]
        inv_J = 1.0/det_J * inv_J
        if dHdrs:
            from numpy import array
            return J_mat, det_J, inv_J, array([dhdr,dhds])
        else:
            return J_mat, det_J, inv_J

    def local_potential_matrix(self, nodes, v, el_id):
        """ 
        This function calculates the local potential matrix for a quad
        element using Gauss Legendre   quadratures as means for 
        integration.
        It returns a matrix that gets added to the global mass matrix.
        
        Parameters:
            nodes:  Array of nodes, attribute node_coord from class Nodes()
            v:     Potential is a vector that contains values for field 
                   parameters such as gravitational forces in  
            el_id: Integer that points to a certain element in the el_set
        output: 
            lo_mass: Array defining the local mass matrix of the problem
                    given by:
                        int_{\Omega} H^T \bar{\bar{\epsilon}} H det(J) 
                        d \Omega_{el}
                    Where H = interpolation functions according to whether
                              the formulation of the problem is scalar or vectorial
        """
        from numpy import zeros, dot
        from scipy.special.orthogonal import p_roots
        node_coords = self.extract_el_points(nodes, el_id)
        if self.order == 2:
            h8 = self.h
            # My own non general purpose integration scheme:
            # Definition of integration degree for each of the coordinates
            deg_r = 3
            deg_s = 3 
            # Generation of Gauss-Legendre points using scypys's function:
            r, weights_r  = p_roots(deg_r)
            s, weights_s  = p_roots(deg_s)
            if self.vectorial:
                lo_v = zeros((16,16))     
            else:
                lo_v = zeros((8,8))
            for i in range(deg_r):
                for j in range(deg_s):
                    J,det_J,invJ = self.numeric_J(node_coords, r[i], s[j]) 
                    if self.vectorial:
                        H8 = zeros((2,16)) 
                        H8[0, 0] = h8[0](r[i], s[j])
                        H8[1, 1] = h8[0](r[i], s[j])
                        for k in range(1,8):
                            H8[0, 2*k] = h8[k](r[i], s[j])
                            H8[1,2*k+1]= h8[k](r[i], s[j])
                    else:
                        H8 = zeros(8)
                        for k in range(8):
                            H8[i] = h8[k](r[i], s[j])                      
                  
                    lo_v = lo_v + weights_r[i]*weights_s[j]*\
                                        dot(H8.transpose(),H8)*det_J 
            lo_v = dot(lo_v,v)
            return lo_v  
            
    def local_mass_matrix(self, nodes, el_id):
        """ 
        This function calculates the local matrix for a quad element using
        Gauss Legendre   quadratures as means for integration.
        It returns a matrix that gets added to the global mass matrix.
        
        Parameters:
            nodes:  Array of nodes, attribute node_coord from class Nodes()
            el_id: Integer that points to a certain element in the el_set
        output: 
            lo_mass: Array defining the local mass matrix of the problem
                    given by:
                        int_{\Omega} H^T \bar{\bar{\epsilon}} H det(J) d \Omega_{el}
                    Where H = interpolation functions according to whether
                              the formulation of the problem is scalar or vectorial
        """
        
        from numpy import zeros, dot
        from scipy.special.orthogonal import p_roots
        epsilon = 1
        node_coords = self.extract_el_points(nodes, el_id)
        if self.order == 2:               
            h8 = self.h
            # My own non general purpose integration scheme:
            # Definition of integration degree for each of the coordinates
            deg_r = 3
            deg_s = 3 
            # Generation of Gauss-Legendre points using scypys's function:
            r, weights_r  = p_roots(deg_r)
            s, weights_s  = p_roots(deg_s)
            if self.vectorial:
                lo_mass = zeros((16,16))     
            else:
                lo_mass = zeros((8,8))
            for i in range(deg_r):
                for j in range(deg_s):
                    J,det_J,invJ = self.numeric_J(node_coords, r[i], s[j]) 
                    if self.vectorial:
                        H8 = zeros((2,16)) 
                        H8[0, 0] = h8[0](r[i], s[j])
                        H8[1, 1] = h8[0](r[i], s[j])
                        for k in range(1,8):
                            H8[0, 2*k] = h8[k](r[i], s[j])
                            H8[1,2*k+1]= h8[k](r[i], s[j])
                    else:
                        H8 = zeros(8)
                        for k in range(8):
                            H8[k] = h8[k](r[i], s[j])                      
                  
                    lo_mass = lo_mass + weights_r[i]*weights_s[j]*epsilon*\
                                        dot(H8.transpose(),H8)*det_J 
            return lo_mass
    def build_local_stiffness(self, nodes, el_id):
        """
         This function calculates the local stiffness matrix for a 
         quad element using Gauss Legendre   quadratures as means for 
         integration.
         It returns a matrix that gets added to the global stiffness matrix.
        
        Parameters:
            nodes:  Array of nodes, attribute node_coord from class Nodes()
            el_id: Integer that points to a certain element in the el_set
        output: 
            lo_stiff: Array defining the local stiffness matrix of the problem
                    given by:
                        \[\integrate_\Omega\nabla^2\mathbf{E} \mathbf{E} d\Omega = \integrate_\Gamma \mathbf{W}\cdot (\grad\mathbf{E}\cdot \hat{n})d\Gamma
                        - \integrate_\Omega \nabla\mathbf{E} : \nabla\mathbf{W} d\Omega  \]
        """
        from numpy import zeros, dot
        from scipy.special.orthogonal import p_roots
        mu = 1
        node_coords = self.extract_el_points(nodes, el_id)
        if self.order == 2:               
            # My own non general purpose integration scheme:
            # Definition of integration degree for each of the coordinates
            deg_r = 3
            deg_s = 3 
            # Generation of Gauss-Legendre points using scypys's function:
            r, weights_r  = p_roots(deg_r)
            s, weights_s  = p_roots(deg_s)
            if self.vectorial:
                lo_stiff = zeros((16,16))     
            else:
                lo_stiff = zeros((8,8))
            for i in range(deg_r):
                for j in range(deg_s):
                    J,det_J,invJ, dHdr = self.numeric_J(node_coords, r[i], s[j], dHdrs = True) 
                    dHdx = dot(invJ, dHdr)                    
                    if self.vectorial:
                        B8 = zeros((4,16)) 
                        B8[0, 0] = dHdx[0, 0]
                        B8[1,1] = dHdx[1, 0]
                        B8[2,0] = dHdx[1, 0]
                        B8[3, 1] = dHdx[0, 0]
                        for k in range(1,8):
                            B8[0, 2*k] = dHdx[0, k]
                            B8[1, 2*k+1]= dHdx[1, k]
                            B8[2, 2*k] = dHdx[1, k]
                            B8[3, 2*k+1]= dHdx[0, k]
                    else:
                        B8 = zeros((2,8))
                        B8[0, 0] = dHdx[0, 0]
                        B8[1,1] = dHdx[1, 0]
                        for k in range(1,9):
                            B8[0, 2*k] = dHdx[0, k]
                            B8[1, 2*k+1]= dHdx[1, k]
                    lo_stiff = lo_stiff + weights_r[i]*weights_s[j]*mu*\
                                        dot(B8.transpose(), B8)*det_J 
            return lo_stiff                 
            
    
class Triangles():
    """
    Triangles are one kind of very important elements
    """    
    def __init__(self, raw_triangles):
        from numpy import shape
        self.el_set = raw_triangles
        self.n_elements = shape(self.el_set)[0]
        #== Assign a local stiffness for the elements given the order===
        if shape(self.el_set)[1] == 4:
            self.order = 1
            
    def extract_el_points(self, _nodes, el_id):
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
        from numpy import zeros
        if 'el_set' in self.__dict__:
            if self.order == 1:
                p1_triangles = self.el_set        
            else:
                raise NotImplementedError('No higher order elements yet,\
                                           there must be an error')
        else:
            raise IOError('No elements parsed, do something else')
            
        pt_a = zeros(2)
        pt_b = zeros(2)
        pt_c = zeros(2)
    
        pt_a[0] = _nodes[p1_triangles[el_id, 1] - 1, 0]
        pt_a[1] = _nodes[p1_triangles[el_id, 1] - 1, 1]
        pt_b[0] = _nodes[p1_triangles[el_id, 2] - 1, 0]
        pt_b[1] = _nodes[p1_triangles[el_id, 2] - 1, 1]
        pt_c[0] = _nodes[p1_triangles[el_id, 3] - 1, 0]
        pt_c[1] = _nodes[p1_triangles[el_id, 3] - 1, 1]
        return [pt_a, pt_b, pt_c]
        
    def extract_el_edges(self, points):
        """
        Receive the list of points given by "extract_el_points" and return 
        a list of line segments from which the area can be calculated
        """
        pt_a, pt_b, pt_c = points 
        #print pt_a, pt_b, pt_c
        ln_ab = pt_b - pt_a    
        ln_bc = pt_c - pt_b
        ln_ca = pt_a - pt_c    
        lines = [ln_bc, ln_ca, ln_ab]        
        return lines
        
    def calc_area(self, lines):
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
    
        
    def build_local_stiffness(self, nodes, el_id):
        if self.order == 1:
            from numpy import zeros, dot
            points = self.extract_el_points(nodes, el_id)
            lines = self.extract_el_edges(points)
            calc_area = self.calc_area
            area = calc_area(lines)
            n = len(lines)
            lo_stiff =  zeros((n, n))  # local (lo) stiffness (stiff) matrix
            for i in range(n):
                for j in range(n):
                    lo_stiff[i, j] = dot(lines[i], lines[j]) / (4*area) 
            return lo_stiff
        else:
            raise NotImplementedError('No higher order triangles for now')
        
    def local_mass_matrix(self, nodes, el_id):
        if self.order == 1:
            from numpy import array, zeros
            from numpy.linalg import det
            tr_mat = zeros((2,2))
            lo_mass = 1./24*array([[2, 1, 1],\
                                   [1, 2, 1],\
                                   [1, 1, 2]])
            points = self.extract_el_points(nodes, el_id)
            tr_mat[:, 0] = points[1] - points[0]
            tr_mat[:, 1] = points[2] - points[0]
            
            jac = det(tr_mat)  
            
            lo_mass = jac * lo_mass    # Escalate acording with the jacobian of the
                                       # transformation 
        else:
            raise NotImplementedError('No higher order triangles for now')
        return lo_mass
        
    #======================== Local Newman matrix  1D P1 ===========================
    def local_newman_matrix():
        from numpy import array
        lo_new = 1./6* array([[2, 1],\
                              [1, 2]])
        return lo_new
    #======================== Local potential matrix  2D P1 =======================
    def local_potential_matrix(self, nodes, el_set, v, el):
        if self.order == 1:
            from numpy import array, zeros
            from numpy.linalg import det
            vec_v = zeros((3,1))   
            vec_v[0, 0] = v[el_set[el, 1]-1]
            vec_v[1, 0] = v[el_set[el, 2]-1]
            vec_v[2, 0] = v[el_set[el, 3]-1]
            tr_mat = zeros((2,2))      # Initiate the transformation (tr) matrix (mat)
            pt_a, pt_b, pt_c = self.extract_el_points(nodes, el)
            tr_mat[:, 0] = pt_b - pt_a
            tr_mat[:, 1] = pt_c - pt_a
            jac = det(tr_mat)  
            lo_v = 1./60*array([[3, 1, 1],\
                                  [1, 3, 1],\
                                  [1, 1, 3]]) # Local (lo) potential (v) matrix
           
            lo_v = jac * (lo_v * vec_v)  # Escalate acording with the jacobian 
                                       # of the transformation 
        
        else:
            raise NotImplementedError('No higher order triangles for now')
        return lo_v
class Boundaries():
    """
    Define what tags from physical entities are associated with certain 
    boundary conditions.
    """
    def __init__(self, dirichlet = {}, newman = {}, bloch = {}):
        self.dirichlet = dirichlet
        self.newman = newman
        self.bloch = bloch           
    def add(self, _bc):
        self.dirichlet = _bc['Dir']
        self.newman = _bc['New']        
        if 'Bloch' in _bc:
            self.bloch = _bc['Bloch']
    def __str__(self):
        return '%s' % self.__dict__
    
    def bloch_multiplication(self, k_x, k_y, nodes, ref_im, *matrices):
        
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
        assert matrices != []
        # For each bloch condition in ref_im 
        im = ref_im[:,1]
        ref = list(set(list(ref_im[:, 0])))
        if matrices[0].shape[0] / 2 == nodes.shape[0]:    
            for i in im[::2]:
                x_im = nodes[ i/2, 0]
                y_im = nodes[ i/2, 1]
                fi_x = exp(1.0j*k_x*x_im)*exp(1.0j*k_y*y_im)
                fi_y = exp(1.0j*k_x*x_im)*exp(1.0j*k_y*y_im)
                for matrix in matrices:
                   # Multiply the column of the image node by the phase factor
                    matrix[:, i] = fi_x * matrix[:, i] 
                    matrix[:, i+1] = fi_y * matrix[:, i+1] 
                    # Multiply the row of the image node by the comlex conjugate
                    #phase factor            
                    matrix[i , :] = fi_x.conjugate()* matrix[i, :]
                    matrix[i + 1, :] = fi_y.conjugate()* matrix[i+ 1,:]
            for i in ref[::2]:      
                x_ref = nodes[ i/2, 0]
                y_ref = nodes[ i/2, 1]
                ff_x = exp(1.0j*k_x*x_ref)*exp(1.0j*k_y*y_ref)
                ff_y = exp(1.0j*k_x*x_ref)*exp(1.0j*k_y*y_ref)
                for matrix in matrices:
                   # and the same for the reference node:  
                    matrix[:, i] = ff_x * matrix[:, i] 
                    matrix[:, i+1] = ff_y * matrix[:, i+1]
                    matrix[i, :] = ff_x.conjugate() * matrix[i, :]
                    matrix[i + 1, :] = ff_y.conjugate() * matrix[i+ 1,:]
        else:
            for i in im:
                x_im = nodes[ i, 0]
                y_im = nodes[ i, 1]
                fi = exp(1.0j*k_x*x_im)*exp(1.0j*k_y*y_im)
                for matrix in matrices:
                   # Multiply the column of the image node by the phase factor
                    matrix[:, i] = fi * matrix[:, i] 
                    # Multiply the column of the image node by the comlex conjugate
                    #phase factor            
                    matrix[i, :] = fi.conjugate()* matrix[i, :]
            for i in ref:        
                x_ref = nodes[ i, 0]
                y_ref = nodes[ i, 1]
                ff = exp(1.0j*k_x*x_ref)*exp(1.0j*k_y*y_ref)
                for matrix in matrices:
                   # and the same for the reference node:  
                    matrix[:, i] = ff * matrix[:, i] 
                    matrix[i, :] = ff.conjugate() * matrix[i,:]
        return matrices
    
    def bloch_sum(self, ref_im, *matrices ):
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
        from numpy import delete, copy
        remove = []
        for k in range(len(ref_im[:, 0])):
            i = ref_im[k, 0]
            j = ref_im[k, 1]
            for matrix in matrices:
                # Sum image node row to reference node row 
                matrix[i, :] = matrix[i, :]+ matrix[j, :]
                # Sum image node column to reference node column
                matrix[:, i] = matrix[:, i]+ matrix[:, j]
                #== stack the values of nodes in vertices for further removal===
                remove.count(j)       
                if remove.count(j) == 0:
                    remove.append(j)
        remove.sort()
        new_matrices = []
        i = 0
        for matrix in matrices:
            new_matrices.append(copy(delete(copy(delete(matrix, remove, 0)), remove, 1)))
            i+=1
        return new_matrices 
class DOF():
    """
    	A degree of freedom to solve for in a explicit methodology.
    	Attributes:
    
    	node_id: 	Number that represents the row that corresponds to a certain node in 
    				the array of nodes.
    	boundary:  Boolean. df: false
    
    	region:  tells to which region it belongs
    
    	vectorial: Tells if the degree of freedom is part of a vectorial formulation.
    			   for 2D vectorial models we will have two dof's per node.
    	comp:  0,1   If vectorial, this will tell if the dof corresponds to x or y component of the field.	
    
    	s_elements:  List of surrounding elements 
    
    	k: 		Componento of the stiffness matrix built using the elements to which it belongs
    
    """
    value = 0.0
    def __init__(self, node_id, simulation, t = None, comp = 0, vectorial = False ):
        self.node_id = node_id
        self.F_i = 0
        self.region = ''
        self.comp = comp
        self.s_elements = []
        self.check_if_in_boundary(simulation, t = t)
        self.boundary = self.check_if_in_boundary(simulation, t = t)
    def __str__(self):
        return "Degree of freedom %s of node %s is shared by elements:\n %s\n"\
    		"and its corresondent stiffness constant is %s \n" %(self.comp, self.node_id,self.s_elements, self.k)		
    def find_surrounding_elements(self, dofs, simulation):
        """
        Loop over elements looking for reference to the current node_id
        """
        self.F_i = 0
        self.s_elements = []
        nodes_coords = simulation.domain.nodes.coords            
        for region in simulation.domain.regions:
            all_elements = region.elements
            for el_class in all_elements:
                elements = all_elements[el_class]
                if elements.vectorial:
                    n_elements = elements.n_elements   
                    el_set = elements.el_set
                    if el_set == []:
                        break
                    else:
                        for el in range(n_elements):
                            if self.node_id + 1 in el_set[el][1:]:
                                self.s_elements.append(region.name+' '+el_class+' '+str(el)+'\n')
                                lo_stif = elements.build_local_stiffness(nodes_coords, el)
                                if simulation.sim_type == 'EM':
                                    mu = region.material_prop['mu']
                                    lo_stif = 1.0/mu * lo_stif
                                elif simulation.sim_type == 'QM':
                                    h = region.material_prop['h']
                                    m = region.material_prop['m']
                                    lo_stif = h/(2.*m)*lo_stif
                                else:
                                    raise NotImplementedError('%s Not a simulation type'%simulation.sim_type) 
                                from numpy import array, dot, where, append
                                u = array([])
                                
                                for node in el_set[el][1:]:
                                    if self.comp == 0:
                                        u = append(u, dofs[2*(node-1)].value)
                                        u = append(u, 0)
                                    else: 
                                        u = append(u, 0)
                                        u = append(u, dofs[2*(node-1)+1].value)
#                                if self.node_id == 107:
#                                    print 'el_set[el][1:]',el_set[el][1:]
#                                    print u,
#                                    for dof in dofs:
#                                        print 'dof.node_id',dof.node_id, 'value', dof.value                                        
                                pivot = where(el_set[el][1:] == self.node_id+1)[0]    
                                k_row = lo_stif[pivot + self.comp,:][0]
                                
                                self.F_i += dot(k_row, u)
                                
                else:
                    n_elements = elements.n_elements
                    el_set = elements.el_set
                    if el_set == []:
                        break
                    else:
                        for el in range(n_elements):
                            if self.node_id + 1 in  el_set[el][1:]:
                                self.s_elements.append(region.name+' '+el_class+' '+str(el)+'\n')
                                lo_stif = elements.build_local_stiffness(nodes_coords, el)
                                if simulation.sim_type == 'EM':
                                    mu = region.material_prop['mu']
                                    lo_stif = 1.0/mu * lo_stif
                                elif simulation.sim_type == 'QM':
                                    h = region.material_prop['h']
                                    m = region.material_prop['m']
                                    lo_stif = h/(2.*m)*lo_stif
                                else:
                                    raise NotImplementedError('%s Not a simulation type'%simulation.sim_type) 
                                from numpy import sum, where       
                                pivot = where(el_set[el][1:] == self.node_id+1)[0]
                                row = lo_stif[pivot + self.comp,:]
                                column1 = lo_stif[:pivot + self.comp, pivot + self.comp]
                                column2 = lo_stif[pivot + self.comp+1:, pivot + self.comp]
                                self.k += sum(row) + sum(column1)+ sum(column2)
        return self.F_i
    def check_if_in_boundary(self, simulation, t = None):
        """
        Run throug line elements checking if this degree of freedom 
        belongs to a bc line defined by a bc condition.
        
        Parameters:
        -----------
        simulation:     Instance of class Simulation()
        
        Returns:
        -------
        
        ln:     the number of the line to which the dof belongs 
        """
        bc_lines = simulation.domain.elements.lines.el_set
        n_lines = simulation.domain.elements.lines.n_elements
        vectorial = simulation.domain.elements.lines.vectorial
        dirichlet = simulation.domain.boundaries.dirichlet
        if vectorial:
            for tag in dirichlet:
                for ln in range(n_lines):
                    if bc_lines[ln, 0] == int(tag):
                        if self.node_id + 1 in bc_lines[ln]:  
                            if self.comp == 0:                        
                                value = dirichlet[tag][0][0]
                            else:
                                value = dirichlet[tag][0][1]
                            if type(value) == str:
                                from math import sqrt, sin, pi 
                                from numpy import isnan
                                x = simulation.domain.nodes.coords[self.node_id, 0]
                                y = simulation.domain.nodes.coords[self.node_id, 1]
                                import re
                                value = re.sub("_"," ",value)
                                exec(value)
                                if isinstance(value, str):
                                    raise TypeError("value should be already" \
                                                    "evaluated something happened")
                                else:
                                    if isnan(value): 
                                        value = 0
                                        print 'value has been reassigned due to \
                                          division by zero'
                                
                                    self.value = value
                            else: 
                                self.value = value
                            return "This DOF belongs to line %s of bc with "\
                            "tag %s.\n Value %s has been assigned to the DOF"\
                            %(ln, tag, dirichlet[tag][0][self.comp])
        else:
        	raise NotImplementedError("Wait please")
        return False         