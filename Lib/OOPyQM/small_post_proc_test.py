# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

@author: santiago
"""

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk
simu = Simulation()
simu.read_solver_input('square8_cap.msh')
simu.domain.read_mesh_file('square8_cap.msh')

inter = Interpreter()
eq = inter.build_harmonic_EM_eq(simu)

my_solver = Solver()
fields = my_solver.solve_stationary(simu, eq)
quads = my_solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]
from numpy import zeros
field3 =  zeros((simu.domain.nodes.n,3))
field3[:,0:2] = fields
fields = field3
#i = 0
#for field in fields:
#    field3 =  zeros((simu.domain.nodes.n,3))
#    field3[:,0:2] = field
#    fields[i] = field3
#    i=i+ 1
#print fields[0]
#for i in range(len(fields)):
#    write_vtk('cap'+str(i)+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
#               quads, ['VECTORS', ['sol'+str(i)], [fields[i]]])
write_vtk('capa'+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields]])
