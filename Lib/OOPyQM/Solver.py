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
    def build_solution(self, dir_solution, g, remove):
        n_sol = g.shape[0]
        j = 0
        for i in range(n_sol):
            if i not in remove:
                g[i] = dir_solution[j]
                j = j+1
        return g
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
    
            solution = zeros((n, n_vects))
                       
            for i in range(n_vects):
               solution[:,i] = self.build_solution(dir_solution[:, i], g, remove)
            
            return v, solution
    
        elif 'n'in solver_param[0] and 'y'in solver_param[1]:
            n_vects = int(solver_param[3])
            v, dir_solution = eigh(equation['left_side'], equation['right_side'], \
                                    eigvals = (0, n_vects-1))
            
            solution = zeros((n, n_vects))
                       
            for i in range(n_vects):
               solution[:,i] = self.build_solution(dir_solution[:, i], g, remove)
            
            return solution
        else:
            print 'error: If you dont want anything why do you solve?'
            
           
       