# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

@author: santiago
"""

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 

write_solver_input('two_cylinders.msh',dimension = 2, bc_type = 'Dir', \
parameter = [], eq = 'Harm_Elec', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 4, 4, 20, 20, 2], bc_filename = 'two_cylinders2.bc')  

simu = Simulation()
simu.read_solver_input('two_cylinders.msh')
simu.domain.read_mesh_file('two_cylinders.msh', True)

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
print g
field2 =  zeros((simu.domain.nodes.n,3))
field2[:,0:2] = g
g_sol = field2
from numpy import sqrt
expr_x = '12*x/sqrt(x**2+y**2)**3-12*(x-4)/sqrt((x-4)**2+y**2)**3'
expr_y = '12*y/sqrt(x**2+y**2)**3-12*y/sqrt((x-4)**2+y**2)**3'
for i in range(g_sol.shape[0]/2):
    if i not in remove:
        x = simu.domain.nodes.coords[i, 0]
        y = simu.domain.nodes.coords[i, 1]
        x_value = eval(expr_x)
        y_value = eval(expr_y)
        g_sol[i,0] = x_value
        g_sol[i,1] = y_value

write_vtk('cylinder2_finer'+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields]])
write_vtk('cylinder2_analy_finer'+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [g_sol]])
