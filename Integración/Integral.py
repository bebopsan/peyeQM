#! /usr/bin/python
## module Integral
"""
    This modules contains functions for the evaluation of
    the integral of a domain defined by points stored in a file with
    GMSH format .msh.

    As input one has to give the name of the file as a string variable. like:

            Integral('parabol.msh')

    The result is the value of the sum of the integral over every element.
"""
from ReadMesh import*
from Read import*
from numpy import hstack,zeros,sqrt,array,dot
from numpy.linalg import det,inv
from types import*

def Integral(archivo):
    
    # Calling of the reading function and retrieve nodes and
    # triangles from the file "archivo".
    if "wrl"in archivo:
        Nodes,Triangles=ReadVRML(archivo)
    if ".msh" in archivo:
        Nodes,Triangles=ReadMesh(archivo)
    
    # Define number of nodes and number of triangles
    nNodes=Nodes.shape[0]
   
    if isinstance(Triangles,list):
        s=zeros((1,array(Triangles).size))
        print array(Triangles)
        s[0,:]=array(Triangles)
        Triangles =s
        nTriangles=1
    else:
        nTriangles=Triangles.shape[0]
        
    if ".msh" in archivo:
        Triangles=Triangles[:,1:4]
  
    # Initialize the matrix for scaling and rotation  needed for
    # the transformation.

    M  =zeros((2,2),dtype=float)
    iM =zeros((2,2),dtype=float)
    Det=0.
    intGlob=0.
    z=zeros(3)
##    print Nodes
    
    for i in range(0,nTriangles):
        M[0,0]=Nodes[Triangles[i,1],0]-Nodes[Triangles[i,0],0]
        M[0,1]=Nodes[Triangles[i,2],0]-Nodes[Triangles[i,0],0]
        M[1,0]=Nodes[Triangles[i,1],1]-Nodes[Triangles[i,0],1]
        M[1,1]=Nodes[Triangles[i,2],1]-Nodes[Triangles[i,0],1]
        
        Det=abs(M[0,0]*M[1,1]-M[1,0]*M[0,1])
       
        iM[0,0] = M[1,1]/Det
        iM[0,1] = -M[1,0]/Det
        iM[1,0] = -M[0,1]/Det
        iM[1,1] = M[0,0]/Det
        
        for j in range(0,3):
            z[j]=Nodes[Triangles[i,j],2]

        dzdr = z[1] - z[0]
        dzds = z[2] - z[0]
        deri=np.array([dzdr,dzds])
        zx=dot(iM[0,:],deri)
        zy=dot(iM[1,:],deri)
        area=sqrt(zx**2+zy**2+1.)/2.
        intGlob=intGlob+area*Det
   
    
    print "The value of the integral over the domain defined by the input is:  ",intGlob
    
