#! /usr/bin/python

from ReadMesh import*
from numpy import hstack,zeros,sqrt,sin
from math import pi
from Write import*  
Nodes,Triangles=ReadMesh("sencHigh.msh")

nNodes=Nodes.shape[0]
nTriangles=Triangles.shape[0]


for i in range(0,nNodes):
    if Nodes[i,0]==0 and Nodes[i,1]==0:
        Nodes[i,2]=1
    else:
        Nodes[i,2]= sqrt(abs(0.0002-(0.025-sqrt(Nodes[i,0]**2+Nodes[i,1]**2)**2)))
        
for i in range(0,nTriangles):
    Triangles[i,0]=3
    for j in range(1,4):
        Triangles[i,j]=Triangles[i,j]-1

WriteVTK("b.vtk",Nodes,Triangles)

WriteMSH("b.msh",Nodes,Triangles)
        
