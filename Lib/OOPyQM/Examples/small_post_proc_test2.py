# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

@author: santiago
"""

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 

write_solver_input('square8_cap.msh',dimension = 2, bc_type = 'Dir', \
parameter = [], eq = 'Harm_Elec', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 4, 4, 20, 20, 2], bc_filename = 'square.bc')  

simu = Simulation()
simu.read_solver_input('square8_cap.msh')
simu.domain.read_mesh_file('square8_cap.msh',simu)

inter = Interpreter()
eq = inter.build_static_EM_eq(simu)
my_solver = Solver()
fields = my_solver.solve_stationary(simu, eq)
quads = my_solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]
from numpy import zeros
field3 =  zeros((simu.domain.nodes.n,3))
field3[:,0:2] = fields
fields = field3

write_vtk('newman_cap'+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields]])
print 'Finished.'