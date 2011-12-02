#! /usr/bin/python
from read_mesh import read_mesh
from PrePro import Mesh1D,Potential1D,meshPlot,Mesh2D
from write import write_msh, write_vtk, write_solver_input
from solver import schroedinger
from math import pi
import numpy as np
from PostPro import plot_1d 
from time import time 
t0 = time()
nodes, elements = Mesh1D('simple',0,2*pi)

write_msh('lineSimple.msh', nodes, elements)

pot = Potential1D('well', nodes, V0=2)

write_solver_input('lineSimple.msh', parameter = pot, bc_type = 'Bloch',\
                      analysis_param = ['y', 'y', 4, 4, 50, 50, 1])

v, d = schroedinger('lineSimple.msh')
tf = time()
z = np.zeros((nodes.size, 3))
z[:, 0] = nodes
nodes = z
#write_vtk('test.vtk', 'thi shit', '', nodes, elements, \
#            ['SCALARS', ['Pot', 'Psi'], [pot, d]])
plot_1d('', nodes = nodes,elements = elements, parameter = pot, bc_type = 'Dir', \
         sol = d)
print 'time elapsed: \n', tf-t0, 'seconds'
##Nodes=Nodes[:,0]
##WriteMSH('lineSimple.msh',Nodes,Elems)
##Nodes,Elems=ReadVTK('test.vtk')


##Schroedinger('lineSimple.msh')

##X=3
##Potential1D('hk',X).__doc__
#Nodes,Elems=Mesh2D("center",0,2,0,2,Nx=30,Ny=30)
##facecolor = 'none'
##meshPlot(Nodes,Elems,facecolor,True,True)
##IntegralH("cone.msh")
##Integral("cone.msh")




