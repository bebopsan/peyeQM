#! /usr/bin/python

from ReadMesh import*
from numpy import hstack,zeros,sqrt,sin
from math import pi
from Write import*  
Nodes,Triangles=ReadMesh("sencMed.msh")

nNodes=Nodes.shape[0]
nTriangles=Triangles.shape[0]


for i in range(0,nNodes):

    if abs(-Nodes[i,0]**2-Nodes[i,1]**2)>0.9999:
            Nodes[i,2]=0

    else:
            Nodes[i,2]= sqrt(1-Nodes[i,0]**2-Nodes[i,1]**2)
    
for i in range(0,nTriangles):
    Triangles[i,0]=3
    for j in range(1,4):
        Triangles[i,j]=Triangles[i,j]-1

WriteVTK("sphere.vtk",Nodes,Triangles)

WriteMSH("sphere.msh",Nodes,Triangles)
        
