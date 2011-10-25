#! /usr/bin/python
from ReadMesh import*
from PrePro import Mesh1D,Potential1D,meshPlot,Mesh2D
from Write import*
from Solver import*
from math import pi
import numpy as np


Nodes,Elems=Mesh1D('simple',0,2*pi)


Pot=Potential1D('morse',Nodes,V0=2)
z=np.zeros((Nodes.size,3))
z[:,0]=Nodes
Nodes=z
WriterVTK('test.vtk','thi shit','',Nodes,Elems,['SCALARS','Potential',Pot])

Nodes=Nodes[:,0]
WriteMSH('lineSimple.msh',Nodes,Elems)
Nodes,Elems=ReadVTK('test.vtk')

WriteSolverInput('lineSimple.msh',parameter=Pot,BCType='Bloch')
Schroedinger('lineSimple.msh')

##X=3
##Potential1D('hk',X).__doc__
#Nodes,Elems=Mesh2D("center",0,2,0,2,Nx=30,Ny=30)
##facecolor = 'none'
##meshPlot(Nodes,Elems,facecolor,True,True)
##IntegralH("cone.msh")
##Integral("cone.msh")




