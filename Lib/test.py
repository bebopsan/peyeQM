#! /usr/bin/python
from ReadMesh import*
from PrePro import Mesh1D,Potential1D,meshPlot,Mesh2D
from Write import*
Nodes,Elems=Mesh1D('simple',0,2)
Pot=Potential1D('poschl',Nodes)
WriteSolverInput('Parab.msh',params=Pot)

##X=3
##Potential1D('hk',X).__doc__
#Nodes,Elems=Mesh2D("center",0,2,0,2,Nx=30,Ny=30)
##facecolor = 'none'
##meshPlot(Nodes,Elems,facecolor,True,True)
##IntegralH("cone.msh")
##Integral("cone.msh")


