# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 08:52:26 2013

@author: santiago
"""


import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 
import matplotlib.pyplot as plt
from matplotlib import rc

filename = 'unit_cell_1-02_two_reg'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 21, 21, 5], bc_filename = filename +'_guide.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh',simu)

inter = Interpreter(vectorial = True)
eq = inter.build_EM_bloch_eq(simu)
my_solver = Solver()
k_coords, energy = my_solver.solve_bloch_Brillouin(simu, eq)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.title("Dispersion relations for Square lattice")
k_div = k_coords.shape[0]/3

from numpy import sqrt, around
print max(energy[:,-1])
for i in range(energy.shape[1]):
    energy[:,i] = sqrt(around(energy[:,i],10))
    ax.plot(energy[:,i])
l = plt.axvline(x=k_div)
l = plt.axvline(x=k_div*2+1)
l = plt.axvline(x= k_coords.shape[0])
plt.ylim(0, max(energy[:,-1]))
plt.xlim(0, k_coords.shape[0])
ax.set_xticks((0,k_div,k_div*2+1, k_coords.shape[0]))
ax.set_xticklabels(('$\Gamma$','$\chi$','$M$','$\Gamma$'),fontsize = 20)
ax.set_ylabel(r'Frequency $ \frac{\omega}{c} $', fontsize = 20)
plt.show()
#from numpy import zeros, sqrt, around
#nodes =  zeros((k_mesh[0].shape[0],3))
#nodes[:,0:2] = k_mesh[0]
#triangles = k_mesh[1]
#
#for i in range(energy.shape[1]):
#    print i, energy[:,i]
#    energy[:,i] = sqrt(around(energy[:,i],10))
#    
#    write_vtk('Bloch_periodic_'+str(i)+'.vtk','this shit','POLYDATA', nodes,\
#               triangles, ['SCALARS', ['sol'], [energy[:,i]]])
#from inform import inform
#inform('Iḿ done runing')
print 'Finished.'