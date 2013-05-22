# -*- coding: utf-8 -*-
""" This script solves the harmonic wave equation for electromagnetic fields in 
a circular waveguide.

It reads the files with name stated in variable filename.

The necessary files are files with following extensions:
    
    .msh    Mesh file has information about nodes and elements
    .reg    Information about regions and material properties
    .bc     Boundary conditions pertinent for the simulation
    .

@author: santiago
"""

import os, sys
lib_path = os.path.abspath('../..')
sys.path.append(lib_path)

from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 
from numpy import pi, sqrt
from numpy import zeros

filename = 'hexagonal_waveguide'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Dir', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 20, 20, 2], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh', True)
simu.domain.read_bc_file(simu.bc_filename)

reg_filename = simu.bc_filename.split('.bc')[0]
simu.domain.read_regions_file(reg_filename)

inter = Interpreter()
eq = inter.build_harmonic_EM_eq(simu)
g = eq['sol_vec']

my_solver = Solver()
value, fields = my_solver.solve_spectral(simu, eq)
print len(fields), 'value', sqrt(value)*2.0
quads = my_solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]

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
    write_vtk(filename +str(i)+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [fields[i]]])

# Analitic solution

a = 2; b = 4
k=[]
for m in range(1,6):
    for n in range(1,6):
        k.append((m*pi/a)**2 +(n*pi/b)**2)
print k