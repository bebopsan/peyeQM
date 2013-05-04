# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 18:58:29 2013

@author: santiago
"""

import os, sys
lib_path = os.path.abspath('../..')
sys.path.append(lib_path)

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 

filename = 'unit_cell_1-02_two_reg'

write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Dir', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 20, 20, 2], bc_filename = filename +'_guide.bc')  


simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh', True)
inter = Interpreter()
eq = inter.build_harmonic_EM_eq(simu)
g = eq['sol_vec']
my_solver = Solver()


value, fields = my_solver.solve_spectral(simu, eq)
print len(fields), 'value', value
quads = my_solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]
from numpy import zeros
for i in range(len(fields)):
    field3 =  zeros((simu.domain.nodes.n,3))
    field3[:,0:2] = fields[i]
    fields[i] = field3
    

#dir_sol = zeros(g.shape[0])
#remove = eq['dir_positions']
#g = my_solver.build_solution(dir_sol, g, remove, True)
#field2 =  zeros((simu.domain.nodes.n,3))
#field2[:,0:2] = g
#g_sol = field2
    write_vtk(filename+'_' +str(i)+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields[i]]])

# Analitic solution
from numpy import pi
a = 2; b = 4
k=[]
for m in range(1,6):
    for n in range(1,6):
        k.append((m*pi/a)**2 +(n*pi/b)**2)
print 'k', k
