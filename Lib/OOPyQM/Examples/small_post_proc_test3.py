# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

@author: santiago
"""

import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 
filename = 'square_periodic'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 21, 21, 1], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh',simu)

inter = Interpreter(vectorial = True)
eq = inter.build_EM_bloch_eq(simu)
my_solver = Solver()
k_mesh, energy = my_solver.solve_bloch(simu, eq)


from numpy import zeros, sqrt, around
nodes =  zeros((k_mesh[0].shape[0],3))
nodes[:,0:2] = k_mesh[0]
triangles = k_mesh[1]

for i in range(energy.shape[1]):
    print i, energy[:,i]
    energy[:,i] = sqrt(around(energy[:,i],10))
    
    write_vtk('Bloch_periodic_'+str(i)+'.vtk','this shit','POLYDATA', nodes,\
               triangles, ['SCALARS', ['sol'], [energy[:,i]]])
from inform import inform
#inform('Iḿ done runing')
print 'Finished.'