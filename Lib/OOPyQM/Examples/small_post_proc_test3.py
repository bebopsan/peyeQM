# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

@author: santiago
"""

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 

write_solver_input('square8_periodic.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'Periodic_Elec', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 4, 4, 20, 20, 2], bc_filename = 'square_bloch.bc')  

simu = Simulation()
simu.read_solver_input('square8_periodic.msh')
simu.domain.read_mesh_file('square8_periodic.msh',simu)

inter = Interpreter()
eq = inter.build_EM_bloch_eq(simu)
my_solver = Solver()
k_mesh, energy = my_solver.solve_bloch(simu, eq)


from numpy import zeros
nodes =  zeros((k_mesh[0].shape[0],3))
nodes[:,0:2] = k_mesh[0]
triangles = k_mesh[1]

write_vtk('Bloch_periodic'+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' , nodes,\
               triangles, ['VECTORS', ['sol'], [energy]])
print 'Finished.'