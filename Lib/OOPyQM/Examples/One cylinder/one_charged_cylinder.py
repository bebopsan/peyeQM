# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 18:58:29 2013

@author: santiago
"""


from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 
filename = 'eight_of_a_cylinder'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Dir', \
parameter = [], eq = 'Harm_Elec', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 4, 4, 20, 20, 2], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh', True)

inter = Interpreter()
eq = inter.build_static_EM_eq(simu)
g = eq['sol_vec']
my_solver = Solver()
fields = my_solver.solve_stationary(simu, eq)
quads = my_solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]
from numpy import zeros
field3 =  zeros((simu.domain.nodes.n,3))
field3[:,0:2] = fields
fields = field3


dir_sol = zeros(g.shape[0])
remove = eq['dir_positions']
g = my_solver.build_solution(dir_sol, g, remove, True)
field2 =  zeros((simu.domain.nodes.n,3))
field2[:,0:2] = g
g_sol = field2

write_vtk(filename +'_finer2.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields]])
