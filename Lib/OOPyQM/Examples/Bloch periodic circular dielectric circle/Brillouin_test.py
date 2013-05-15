# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 08:52:26 2013
Solves bloch periodic boundary value problem for 2D domains of periodic 
arrays of dielectric rods

@author: santiago
"""

import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)
from Classes import Simulation
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 
import matplotlib.pyplot as plt
from matplotlib import rc
from inform import inform
from numpy import sqrt, around, pi, savetxt
import pickle 

filename = 'unit_cell_1-04_two_reg'
#filename = 'empty_cell'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 8, 8, 31, 31, 1], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh',simu)

inter = Interpreter(vectorial = True)
eq = inter.build_EM_bloch_eq(simu)
print simu.domain.regions[1]

my_solver = Solver()
k_coords, energy = my_solver.solve_bloch_Brillouin(simu, eq)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.title("Dispersion relations for Square lattice")
k_div = k_coords.shape[0]/3.0

f = open('result.gz','w')
a = 1
print max(energy[:,-1])
for i in range(energy.shape[1]):
    energy[:,i] = sqrt(around(energy[:,i],10))*a/(2*pi)
    savetxt(f, around(energy[:,i],10),comments='energy profile '+str(i))
    ax.plot(energy[:,i])
l = plt.axvline(x=k_div)
l = plt.axvline(x=k_div*2+1)
l = plt.axvline(x= k_coords.shape[0])
plt.ylim(0, max(energy[:,-1]))
plt.xlim(0, k_coords.shape[0])
ax.set_xticks((0,k_div,k_div*2+1, k_coords.shape[0]))
ax.set_xticklabels(('$\Gamma$','$\chi$','$M$','$\Gamma$'),fontsize = 20)
ax.set_ylabel(r'Frequency $ \frac{\omega a}{2\pi c} $', fontsize = 20)
#plt.savefig("img/r_01a_inHomogeneus_material_true.svg")
plt.show()
f.close
pl = open(filename + '.pkl','wb')
pickle.dump([energy, k_div], pl)
pl.close()

inform('Iḿ done runing')
print 'Finished.'