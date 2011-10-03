#! /usr/bin/python

from ReadMesh import*
from numpy import hstack,zeros
from Write import*  
Nodes,Triangles=ReadMesh("parabol.msh")

nNodes=Nodes.shape[0]
nTriangles=Triangles.shape[0]


for i in range(0,nNodes):
    Nodes[i,2]=Nodes[i,0]**2+Nodes[i,1]**2 

for i in range(0,nTriangles):
    Triangles[i,0]=3
    for j in range(1,4):
        Triangles[i,j]=Triangles[i,j]-1

WriteVTK("Parab.vtk",Nodes,Triangles)

WriteMSH("Parab.msh",Nodes,Triangles)
        
