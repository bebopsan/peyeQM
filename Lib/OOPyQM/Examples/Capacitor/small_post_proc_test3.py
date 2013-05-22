# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:31:57 2013

Pseudo:
-import stuff
- define filename root
- Modify mesh file in order to add solver parameters used as input for the 
  definition of class Simulation().
  Particularly define: filename, dimension, type of equation, kind of 
  boundary conditions, kind of solver, and a list of analysis parameters, 
  used to compute.
- Instantiate stuff:
    Simulation
        Read the msh file in order to save configuration parameters.
        Add a domain to the simulation by reading from the msh file.
    Interpreter
        Extract the equation to be solved from the simulation using 
        the interpreter.
    Solver
        Solve the equation passed from the interpreter, and return 
        what appeared.
- Make the wavenumber mesh preety for writing
- For each mode solution print a separate vtk file.
- Send a tweet when finished.
@author: santiago
"""

import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)
from Classes import Simulation, DOF
from Interpreter import Interpreter
from Solver import Solver
from write import write_vtk, write_solver_input 

filename = 'capacitor'
write_solver_input(filename +'.msh',dimension = 2, bc_type = 'Bloch', \
parameter = [], eq = 'EM', sol_type = 'Stationary',analysis_param \
= ['y', 'y', 15, 15, 11, 11, 1], bc_filename = filename +'.bc')  

simu = Simulation()
simu.read_solver_input(filename +'.msh')
simu.domain.read_mesh_file(filename +'.msh',simu)
simu.domain.read_bc_file(simu.bc_filename)
reg_filename = simu.bc_filename.split('.bc')[0]
simu.domain.read_regions_file(reg_filename)
node_coords = simu.domain.nodes.coords


dt = 0.1 
n_timesteps = 50
times = range(2, n_timesteps)

dofs_past = []
dofs_present = []
dofs_future = []

inter = Interpreter(vectorial =True)
mass = inter.lumped_mass_matrix(simu)
mass = mass[:,0]
print sum(mass)
for node in range(node_coords.shape[0]):
    dofs_past.append(DOF(node, simu, comp = 0, t = 0 * dt))
    dofs_past.append(DOF(node, simu, comp = 1, t = 0 * dt))
    dofs_present.append(DOF(node, simu, comp = 0, t = 1 * dt))
    dofs_present.append(DOF(node, simu, comp = 1, t = 1 * dt))
    dofs_future.append(DOF(node, simu, comp = 0, t = 2 * dt))
    dofs_future.append(DOF(node, simu, comp = 1, t = 2 * dt))

snapshot = [dofs_past, dofs_present, dofs_future]

from numpy import zeros, append, array
solver = Solver()
quads = solver.substract_1(simu.domain.elements.quads.el_set)
quads = quads[:,1:]
field3 =  zeros((simu.domain.nodes.n,3))
for n in times:
    print n
    field = array([])
#    for dof in snapshot[1]:
#        print 'dof.node_id',dof.node_id, 'value', dof.value    
    for dof in range(len(snapshot[0])):
        E_past =  snapshot[0][dof]      
        E_present = snapshot[1][dof]
        E_future = snapshot[2][dof]  
        if E_future.check_if_in_boundary(simu, t = n*dt) == False:
            Fi = E_present.find_surrounding_elements(snapshot[1], simu)
            F = -Fi + 1.0/(dt**2)*mass[dof]*(2*E_present.value - E_past.value)
            snapshot[2][dof].value = dt**2/mass[dof] * F

        snapshot[0][dof] = E_present
        snapshot[1][dof] =  E_future
        field = append(field, snapshot[2][dof].value)
        #print 'snapshot[2]['+str(dof)+'].value',snapshot[2][dof].value
        field_past = append(field, snapshot[0][dof].value)
        field_present =  append(field, snapshot[2][dof].value)
        
    field = solver.build_solution_2(field)
    field3[:,0:2] = field
    field = field3
    write_vtk(filename+'_' +str(n)+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
               quads, ['VECTORS', ['sol'], [field]])

#    print n, n*dt    
#    #print snapshot[2][100].check_if_in_boundary(simu)
#    field3[:,0:2] = fields[i]
#
#for i in range(len(fields)):
#    field3 =  zeros((simu.domain.nodes.n,3))
#    field3[:,0:2] = fields[i]
#    fields[i] = field3
#write_vtk(filename+'_' +str(i)+'.vtk', 'MyTitle', 'UNSTRUCTURED_GRID' ,simu.domain.nodes.coords,\
#               quads, ['VECTORS', ['sol'], [fields[i]]])
#eq = inter.build_EM_bloch_eq(simu)
#print simu.domain.regions[1]
#
#my_solver = Solver()
#k_mesh, energy = my_solver.solve_bloch(simu, eq)
#
#
#from numpy import zeros, sqrt, around
#nodes =  zeros((k_mesh[0].shape[0],3))
#nodes[:,0:2] = k_mesh[0]
#triangles = k_mesh[1]
#
#for i in range(energy.shape[1]):
#    energy[:,i] = sqrt(around(energy[:,i],10))
#    
#    write_vtk('Bloch_periodic_r-a_1-04_two_reg_true'+str(i)+'.vtk','this shit','POLYDATA', nodes,\
#               triangles, ['SCALARS', ['sol'], [energy[:,i]]])
#
from inform import inform
#inform('done runing time simulation ')
print 'Finished.'