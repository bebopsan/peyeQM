#! /usr/bin/python
"""
Created on Tue Feb 12 09:51:01 2013

@author: santiago
"""

from Classes import *
from Interpreter import *
from Solver import *
from utils import substract_1  
from write import write_vtk

simulation = Simulation()
simulation.read_solver_input('square.msh')
simulation.domain.read_mesh_file('square.msh')
simulation.domain.read_bc_file('square.bc')

interpreter = Interpreter()
equation = interpreter.build_QM_dirichlet_eq(simulation)

solver = Solver()
solution = solver.solve_spectral(simulation, equation)

nodes_coords = simulation.domain.nodes.coords
triangles = simulation.domain.elements.triangles.el_set
triangles = substract_1(triangles)
triangles = triangles[:, 1:]

write_vtk('square_from_OOP.vtk','this shit','',nodes_coords,\
            triangles,['SCALARS',['solution'],[solution[1]]])