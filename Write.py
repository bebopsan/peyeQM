#! /usr/bin/python
## module Write
'''
    WriteVTK(points,triangles,gradient,Output).
   f
    points:
                Array of coordinates for each node.
    triangles:
                Array of conectivities between nodes. For this case
                sets of three nodes.
    gradient:
                Array of components from the gradient for each node.
            
    Output:
            String with the name of the VTK data fsile to be created.
            This file contains the information from the nodes and elements
            from the imput with the calculated gradient on each node.
            
'''
from numscan import numscan
import numpy as np
from numpy import linspace as lin
from numpy import hstack
def WriteVTK(Output,points,triangles,grad=0):

    nEle=str(np.shape(triangles)[0])
    nPoi=str(np.shape(points)[0])
    sizeEle=str(np.size(triangles))
    f = open(Output, 'w')
    f.write('# vtk DataFile Version 3.1 \n')
    f.write('Results\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    
    f.write(('POINTS '+nPoi+' DOUBLE\n'))
    np.savetxt(f,points,fmt='%6.6f')
    f.write(('TRIANGLE_STRIPS '+nEle+' '+sizeEle+' \n'))

    np.savetxt(f,triangles,fmt='%d')
    if grad!=0:
        f.write(('POINT_DATA '+nPoi+' \n'))
        f.write('VECTORS flechas double\n')
        np.savetxt(f,gradient,fmt='%+f %+f %d')
def WriteMSH(Output,points,triangles):

    nEle=str(np.shape(triangles)[0])
    nPoi=str(np.shape(points)[0])
    sizeEle=str(np.size(triangles))
    f = open(Output, 'w')
    f.write('$MeshFormat\n')
    f.write('2.2 0 8\n')
    f.write('$EndMeshFormat\n')
    f.write('$Nodes\n')
    vec=np.zeros((np.shape(points)[0],1))
    for i in range(0,np.shape(points)[0]):
        vec[i,0]=i+1
    
    points=hstack((vec,points))
    f.write(nPoi+'\n')
    np.savetxt(f,points,fmt='%d %f %f %f')
    f.write('$EndNodes\n')
    f.write('$Elements\n')
    f.write(nEle+'\n')
    vec=np.zeros((np.shape(triangles)[0],5))
    for i in range(0,np.shape(triangles)[0]):
        vec[i,0]=i+1
        vec[i,1:5]=[2,2,0,6]
    triangles =triangles[:,1:5]
    triangles=hstack((vec,triangles))
    print triangles.shape
    np.savetxt(f,triangles,fmt='%d %d %d %d %d %d %d %d')
    f.write('$EndElements\n')
