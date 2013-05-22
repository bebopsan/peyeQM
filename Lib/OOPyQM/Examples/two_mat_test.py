# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 10:55:57 2013

@author: santiago
"""
import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
from Classes import Simulation
from Interpreter import Interpreter
from write import write_solver_input 

filename = 'periodic_2_mat'

write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 4, 4, 1], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh',simu)
simu.domain.read_regions_file(filename)
print simu.domain.regions[0].elements['cuad_Quads'].el_set
print simu.domain.regions[1].elements['cuad_Quads'].el_set
print simu.domain.regions[0].elements['cuad_Lines'].el_set

inter = Interpreter(vectorial = True)
#stif = inter.global_stiffness_matrix(simu)
#print stif

