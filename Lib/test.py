#! /usr/bin/python
from ReadMesh import*
from PrePro import Mesh1D,Potential1D,meshPlot,Mesh2D
from Write import*
from Solver import*
from math import pi
Nodes,Elems=Mesh1D('simple',0,2*pi)
Pot=Potential1D('xwell',Nodes,V0=0.02)
WriteMSH('lineSimple.msh',Nodes,Elems)
WriteSolverInput('lineSimple.msh',parameter=Pot,BCType='Bloch')
Schroedinger('lineSimple.msh')

##X=3
##Potential1D('hk',X).__doc__
#Nodes,Elems=Mesh2D("center",0,2,0,2,Nx=30,Ny=30)
##facecolor = 'none'
##meshPlot(Nodes,Elems,facecolor,True,True)
##IntegralH("cone.msh")
##Integral("cone.msh")




