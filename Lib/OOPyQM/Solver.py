# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 08:41:20 2013

Module Solver shall contain all the solver algorithms needed to 
extract the solution out of a system of equations produced by the 
interpreter. 
"""
__author__ = ['Santiago Echeverri Chacón']

class Solver():
    """
    Instances from this class will take information out of the simulation
    and the system of equations and build a solution for the problem.
    """
    def build_solution(self, dir_solution, g, remove, vectorial = False):
        n_sol = g.shape[0]
        j = 0
        if vectorial:
            from numpy import zeros
            vec_field = zeros((n_sol/2,2))
            for i in range(n_sol):
                if i not in remove:
                    g[i] = dir_solution[j]
                    j = j+1
            for i in range(0,n_sol, 2):
               if i == 0:
                   vec_field[0] = [g[i], g[i+ 1]] 
               else:
                   vec_field[i/2] = [g[i], g[i+ 1]] 
               
            return vec_field
        else:
            for i in range(n_sol):
                if i not in remove:
                    g[i] = dir_solution[j]
                    j = j+1
            return g
    def build_solution_2(self, solution):
        """
        The other build solution method was good for cases where we had
        a dirichlet vector. This method is implemented for the
        explicit solution where I already have known and discovered 
        values in the same vector.
        The problem is assumed vectorial
        """
        n_sol = solution.shape[0]
        from numpy import zeros
        vec_field = zeros((n_sol/2,2))
        for i in range(0,n_sol, 2):
           if i == 0:
               vec_field[0] = [solution[i], solution[i+ 1]] 
           else:
               vec_field[i/2] = [solution[i], solution[i+ 1]] 
        return vec_field
    def solve_stationary(self, simulation, equation):
        from scipy.linalg import solve        
        g = equation['sol_vec']
        remove = equation['dir_positions']
        dir_solution = solve(equation['left_side'], equation['right_side'])
        print 'Solving...'
        if equation['vectorial']:
            solution = self.build_solution(dir_solution, g, remove, True)
        else:
            solution = self.build_solution(dir_solution, g, remove)
        return solution
        
    def solve_spectral(self, simulation, equation):
        """
        This method solves for the eigen-values (lam) and eigen_vectors 
        ({u}) of an equation in the frequency domain of the form:
            [K]{u}=lam[M]{u}.
        """
        from numpy import zeros
        from scipy.linalg import eigh, eigvalsh
        n = simulation.domain.nodes.n 
        solver_param = simulation.solver_param
        g = equation['sol_vec']
        remove = equation['dir_positions']
        print 'Solving eigenvalue problem...\n'
        if 'y'in solver_param[0] and 'n' in solver_param[1]:
            n_vals = int(solver_param[2])
            v = eigvalsh(equation['left_side'], equation['right_side'], \
                                    eigvals = (0, n_vals-1))
    #                v = v/2
            print 'The Eigenvalues are:\n', v
            return v
    
        elif 'y'in solver_param[0] and 'y'in solver_param[1]:
            n_vals = int(solver_param[2])
            n_vects = int(solver_param[3])
            n_solutions = max(n_vals,n_vects)
            v, dir_solution = eigh(equation['left_side'], equation['right_side'], \
                                     eigvals = (0, n_solutions-1))
    #                v = v/2
            if equation['vectorial']:
                solution = []
                for i in range(n_vects):
                    solution.append(self.build_solution(dir_solution[:, i], g, remove, True))
            else:
                solution = zeros((n, n_vects))
                for i in range(n_vects):
                    solution[:,i] = self.build_solution(dir_solution[:, i], g, remove)
            
            return v, solution
    
        elif 'n'in solver_param[0] and 'y'in solver_param[1]:
            n_vects = int(solver_param[3])
            v, dir_solution = eigh(equation['left_side'], equation['right_side'], \
                                    eigvals = (0, n_vects-1))
            
            if equation['vectorial']:
                solution = zeros((n/2, n_vects))
                for i in range(n_vects):
                    solution[:,i] = self.build_solution(dir_solution[:, i], g, remove, True)
            else:
                solution = zeros((n, n_vects))
                for i in range(n_vects):
                    solution[:,i] = self.build_solution(dir_solution[:, i], g, remove)
            
            return solution
        else:
            raise IOError('If you dont want anything why do you solve?')
    def solve_bloch(self, simulation, equation):
        """
            This method solves a 2D bloch periodic condition problem 
            given an equation and that was interpreted used an instance of
            the Interpreter Class.
            Right now I'm testing it for solution of EM probleems with 
            Bloch periodicity such as those encountered in perfect 
            photonic crystals. 
            
            Parameters:
            -----------    
            simulation:     An instance of class Simulation
                            from which all the information about the
                            domain and problem conditions will be extracted.
                            The method read_solver_input() must have been 
                            called with an appropiate input file in order 
                            to extract information such as solver parameters
                            from which the wave number is extracted.
            equation:   Output of one methods of an instance from 
                        class Interpreter(). Initialy only the output from 
                        method: build_EM_bloch_eq(self, simulation) works.
                        Future implementation of a simmilar method for QM 
                        is pending.
            Returns:
            --------
            k_mesh:     A simple triangular mesh made from x and y values 
                        of the wavenumber.
            energy: Set of eigenvalues assigned for each of the pairs 
                    (k_x, k_y) in the k_mesh. 
        """
        #============= Load methods and attributes ================== 
        from numpy import zeros, pi, linspace 
        from scipy.linalg import eigvalsh
        bloch_multiplication = simulation.domain.boundaries.bloch_multiplication
        bloch_sum = simulation.domain.boundaries.bloch_sum
        ref_im_mul = equation['ref_im_mul']
        stif = equation['left_side']
        mass = equation['right_side']
        analysis_param = simulation.solver_param
        nodes = simulation.domain.nodes.coords 
        analysis_param = simulation.solver_param
        #============= Discretize the wavenumber dommain ================== 
        analysis_param = simulation.solver_param
        nk_x = int(analysis_param[4]) # number of k to sweep in xpasare a windows entonces para usa
        nk_y = int(analysis_param[5]) # number of k to sweep in y
        k_max = pi/float(analysis_param[6])
        k_min = -k_max
        k_range_x = linspace(k_min, k_max, num = nk_x)
        k_range_y = linspace(k_min, k_max, num = nk_y)
        n_vals = int(analysis_param[2])        
        #================ Initiate global variable =======================
        energy = zeros( (nk_x * nk_y, n_vals) )
        
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
                vals = eigvalsh(new_stif, new_mass, eigvals = (0, n_vals-1))
                for val in range(0, n_vals):
                    energy[i, val] = vals[val]
                
                print i+1, ' points out of :', nk_x * nk_y, '\n'  
                i = i+1
        k_mesh = self.meshsq2D(k_min, k_max, k_min, k_max, nk_x, nk_y)
        return k_mesh, energy       
    def solve_bloch_Brillouin(self, simulation, equation):
        """
            This method solves a 2D bloch periodic condition problem 
            given an equation and that was interpreted used an instance of
            the Interpreter Class.
            Right now I'm testing it for solution of EM probleems with 
            Bloch periodicity such as those encountered in perfect 
            photonic crystals. 
            
            Parameters:
            -----------    
            simulation:     An instance of class Simulation
                            from which all the information about the
                            domain and problem conditions will be extracted.
                            The method read_solver_input() must have been 
                            called with an appropiate input file in order 
                            to extract information such as solver parameters
                            from which the wave number is extracted.
            equation:   Output of one methods of an instance from 
                        class Interpreter(). Initialy only the output from 
                        method: build_EM_bloch_eq(self, simulation) works.
                        Future implementation of a simmilar method for QM 
                        is pending.
            Returns:
            --------
            k_mesh:     A simple triangular mesh made from x and y values 
                        of the wavenumber.
            energy: Set of eigenvalues assigned for each of the pairs 
                    (k_x, k_y) in the k_mesh. 
        """
        #============= Load methods and attributes ================== 
        from numpy import zeros, pi
        from scipy.linalg import eigvalsh
        bloch_multiplication = simulation.domain.boundaries.bloch_multiplication
        bloch_sum = simulation.domain.boundaries.bloch_sum
        ref_im_mul = equation['ref_im_mul']
        stif = equation['left_side']
        mass = equation['right_side']
        analysis_param = simulation.solver_param
        nodes = simulation.domain.nodes.coords 
        analysis_param = simulation.solver_param
        #============= Discretize the wavenumber dommain ================== 
        analysis_param = simulation.solver_param
        nk_x = int(analysis_param[4]) # number of k to sweep in xpasare a windows entonces para usa
        nk_y = int(analysis_param[5]) # number of k to sweep in y
        assert nk_x == nk_y, "Must be the same number of divisions"
        k_max = pi/float(analysis_param[6])
        k_coords = self.mesh_BrillouinZone(k_max, k_max,nk_x, nk_y)
        print'k_coords.shape[0]', k_coords.shape[0]
        n_vals = int(analysis_param[2])        
        #================ Initiate global variable =======================
        energy = zeros( (k_coords.shape[0], n_vals) )
        
        #======================= Main cycles ===============================
        i = 0         
        print 'Calculating each of the ', k_coords.shape[0], \
               ' points in k plane...\n'
        for k in k_coords:
            k_x = k[0]
            k_y = k[1]
            new_stif, new_mass = bloch_multiplication(k_x, \
                                            k_y, nodes, \
                                            ref_im_mul, stif.copy(), \
                                            mass.copy())
            new_stif, new_mass = bloch_sum(ref_im_mul, new_stif, new_mass)
            vals = eigvalsh(new_stif, new_mass, eigvals = (0, n_vals-1))
            for val in range(0, n_vals):
                energy[i, val] = vals[val]
            
            print i+1, ' points out of :', k_coords.shape[0], '\n'  
            i = i+1
        
        return k_coords, energy           
        
                 
    def substract_1(self, matrix):
        """ 
            Substracts the number 1 from node and element positions 
            as to convert gmsh base 1 numbering to vtk zero base.
            
            Parameter: 
                matrix:     Arbitrary size array that may represent elements
                            or nodes from which the number one has to be 
                            extracted.
            Returns:
                matrix: Input matrix where the number 1 has been substracted 
                        from each position.
        """
        n_rows = matrix.shape[0]
        n_cols = matrix.shape[1]
        for i in range(n_rows):
            for j in range(1, n_cols):
                matrix[i,j] = matrix[i,j] - 1
        return matrix
    def meshsq2D(self, xmin, xmax, ymin, ymax, nxpoints, nypoints):
        """ Generate a 2D mesh where the points are equally spaced.

        :Parameters:
            xmin:  initial value of the rectangular domain over x axis
            xmax:  final value of the rectangular domain over x axis
            ymin:  initial value of the rectangular domain over y axis
            ymax:  final value of the rectangular domain over y axis
            nxpoints:     Number of divisions over x
            nypoints:     Number of divisions over y

        :Returns:
            coords:  numpy array like matrix of the discretized domain 
                     with shape Nx Ny
            elems:   numpy array like matrix of the relations between nodes.
    
        Last modification: April 30, 2013
    
        """    
        from numpy import zeros
        coordx, elem = self.mesh1D(xmin,xmax,nxpoints)
        coordy, elem = self.mesh1D(ymin,ymax,nypoints)
        npoints = nxpoints*nypoints
        coords = zeros ( (npoints,2),dtype=float)
        cont = 0
        for j in range(0,nypoints):
            for i in range(0,nxpoints):
                coords[cont,0] = coordx[i]
                coords[cont,1] = coordy[j]
                cont = cont+1
        nelems = (nxpoints - 1)*(nypoints - 1)
        elems = zeros ( (nelems,4) )
        cont = 0
        for j in range(0,nypoints - 1):
            for i in range(0,nxpoints - 1):
                elems[cont,0] = j*nxpoints + i
                elems[cont,1] = j*nxpoints + i + 1
                elems[cont,2] = nxpoints + j*nxpoints + i + 1
                elems[cont,3] = nxpoints + j*nxpoints + i
                cont = cont + 1
        return coords, elems             
    def meshtr2D(self, xmin,xmax,ymin,ymax,nxpoints,nypoints):
        """ Generate a 2D mesh where the points are equally spaced.

        Parameters:
        -----------
        xmin:  initial value of the rectangular domain over x axis
        xmax:  final value of the rectangular domain over x axis
        ymin:  initial value of the rectangular domain over y axis
        ymax:  final value of the rectangular domain over y axis
        nxpoints:     Number of divisions over x
        nypoints:     Number of divisions over y


        Returns:
        --------
        coords:  numpy array like matrix of the discretized domain with shape Nx Ny
        elems:   numpy array like matrix of the relations between nodes.
    
        Last modification: date 21/10/2011
    
        """    
        from numpy import zeros
        coordx, elem = self.mesh1D(xmin,xmax,nxpoints)
        coordy, elem = self.mesh1D(ymin,ymax,nypoints)
        npoints = nxpoints*nypoints
        coords = zeros ( (npoints,2),dtype=float)
        cont = 0
        for j in range(0,nypoints):
            for i in range(0,nxpoints):
                coords[cont,0] = coordx[i]
                coords[cont,1] = coordy[j]
                cont = cont+1
        nelems = 2*(nxpoints-1)*(nypoints-1)
        elems = zeros ( (nelems,3) )
        cont = 0
        for j in range(0,nypoints-1):
            for i in range(0,nxpoints-1):
                elems[cont,0] = j*nxpoints+i
                elems[cont,1] = j*nxpoints+i+1
                elems[cont,2] = nxpoints + j*nxpoints+i+1
                elems[cont+1,0] = j*nxpoints+i
                elems[cont+1,1] = nxpoints + j*nxpoints+i+1
                elems[cont+1,2] = nxpoints + j*nxpoints+i
                cont = cont+2
        return coords, elems
        
    def mesh1D(self, xmin, xmax, npoints):
        """
    
            Generate a 1D mesh where the points ares equally spaced.
    
            Parameters:
            -----------
            xmin:    coordinate of the beginning of the line segment
            xmax:    coordinate of the end of the line segment
            npoints: number of subdivisions
    
    
            Returns:
            --------
            coords:  numpy array like vector of the discretized domain with lenght N
            elems:  numpy array like vector of the relations between nodes.
        
        
            Raises:
    	-------
        
       
            Last modification: date 21/10/2011
        
        """
        from numpy import linspace, zeros
        coords = linspace(xmin,xmax,npoints,endpoint=True)
        elems = zeros( (npoints-1,2) ,dtype=int )
        for i in range(0,npoints-1):
            elems[i,0] = i
            elems[i,1] = i+1
        return coords, elems
    def mesh_BrillouinZone(self, xmax, ymax, nxpoints,nypoints):
        """
        This functions returns the contour of the Brillouin zone for a 
        square shaped Brillouin zone.
        """
        from numpy import zeros
        coordx, line_x = self.mesh1D(0, xmax,nxpoints)
        coordy, line_y = self.mesh1D(0, ymax,nypoints)
        assert nxpoints == nypoints 
        npoints = 3*nxpoints - 1
        coords = zeros ( (npoints,2),dtype=float)
        coords[:nxpoints, 0] = coordx
        coords[nxpoints:nypoints+nxpoints, 0] = coordx[-1]
        coords[nxpoints:nypoints+nxpoints, 1] = coordy
        coords[2*nxpoints:-1,0] = coordx[-2:0:-1]
        coords[2*nxpoints:-1,1] = coordy[-2:0:-1]
        coords[-1] = coords[0]
        return coords