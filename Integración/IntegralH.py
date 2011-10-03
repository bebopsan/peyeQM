#! /usr/bin/python
## module IntegralH
"""
    This modules contains functions for the evaluation of
    the integral of a domain defined by points stored in a file with
    GMSH format .msh.

    As input one has to give the name of the file as a string variable. like:

            Integral('parabol.msh')

    The result is the value of the sum of the integral over every element.
"""
from ReadMesh import*
from numpy import hstack,zeros,sqrt
from numpy.linalg import det

def IntegralH(archivo):
    
    # Calling of the reading function and retrieve nodes and
    # triangles from the file "archivo".
    if "wrl"in archivo:
        [Nodes,Triangles]=ReadVRML(archivo)
    if ".msh" in archivo:
        Nodes,Triangles=ReadMesh(archivo)
    
    # Define number of nodes and number of triangles
    nNodes=Nodes.shape[0]
    nTriangles=Triangles.shape[0]
    Triangles=Triangles[:,1:4]

    # Initialize the matrix for scaling and rotation  needed for
    # the transformation.

    M=zeros((2,2))
    Det=0.
    intGlob=0.
##    print Nodes
    
    for i in range(0,nTriangles):
       
        magAx=(Nodes[Triangles[i,1],0]-Nodes[Triangles[i,0],0])**2
        magAy=(Nodes[Triangles[i,1],1]-Nodes[Triangles[i,0],1])**2
        magAz=(Nodes[Triangles[i,1],2]-Nodes[Triangles[i,0],2])**2
        magBx=(Nodes[Triangles[i,2],0]-Nodes[Triangles[i,0],0])**2
        magBy=(Nodes[Triangles[i,2],1]-Nodes[Triangles[i,0],1])**2
        magBz=(Nodes[Triangles[i,2],2]-Nodes[Triangles[i,0],2])**2
        magCx=(Nodes[Triangles[i,2],0]-Nodes[Triangles[i,1],0])**2
        magCy=(Nodes[Triangles[i,2],1]-Nodes[Triangles[i,1],1])**2
        magCz=(Nodes[Triangles[i,2],2]-Nodes[Triangles[i,1],2])**2

        a=sqrt(magAx+magAy+magAz)
        b=sqrt(magBx+magBy+magBz)
        c=sqrt(magCx+magCy+magCz)

        s=(a+b+c)/2
        area=sqrt(s*(s-a)*(s-b)*(s-c))
        
              
        intGlob=intGlob+area
   
    
    print "The value of the integral over the domain defined by the input is:  ",intGlob
    
