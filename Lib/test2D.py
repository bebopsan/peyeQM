#! /usr/bin/python
from ReadMesh import*
from PrePro import Mesh1D, potential_2d, meshPlot, Mesh2D
from Write import*
from Solver import*
from math import pi
import numpy as np
from PostPro import* 


nodes, lines, tria = Readmsh('square.msh')

pot = potential_2d('well', nodes, v0 = 2)

print pot, nodes, tria 

#WriteSolverInput('lineSimple.msh', parameter = pot, BCType='Dir')

v, d = Schroedinger('', Nodes = nodes, Elems = elems, parameter = pot, \
                    Dimension = 2, BCType = 'Dir', Type = 'Stationary', \
                    Eq = 'Schro', AnalisisParam = ['y', 'y', 4, 4, 101])


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




