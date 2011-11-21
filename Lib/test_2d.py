#! /usr/bin/python
from read_mesh import read_mesh, read_solver_input
from PrePro import potential_2d
from write import write_solver_input, write_vtk
from numpy import shape, zeros
from solver import schroedinger
from utils import substract_1
from os import system
from inform import inform
glo_tag = 'test1' # Global tag for referencing all files

nodes, elements = read_mesh(glo_tag +'.msh')

n = shape(nodes)[0]                 #Number of nodes      
bc_lines = elements[1]      # Boundary condition is (bc) lines
triangles = elements[2]

potential = potential_2d('well', nodes, v0 = 2)



write_solver_input(glo_tag +'.msh', parameter = potential, dimension = 2, \
                   bc_type = 'Bloch', \
                   sol_type = 'Stationary', eq = 'Schro', \
                   analysis_param = ['y', 'y', 4, 1, 20, 20, 1], \
                   bc_filename = 'test2.bc')

#k = read_solver_input(glo_tag +'.msh')
#v, solution = schroedinger(glo_tag +'.msh')
#triangles = substract_1(triangles)
#triangles = triangles[:, 1:]

k, solution = schroedinger(glo_tag +'.msh')
z = zeros((k[0].shape[0], 3))
print k[0].shape, z[:, 0:1].shape
z[:, 0:2] = k[0]
nodes = z
triangles = k[1]
 

write_vtk(glo_tag + '_5.vtk', 'this shit', '', nodes, triangles, \
            ['SCALARS', ['solution'], [solution]])
            
#inform('Termine')
#system('shutdown -P now')
