#! /usr/bin/python
## module Write
'''
    WriteVTK(Nodes,Elems,gradient,Output).
   f
    Nodes:
                Array of coordinates for each node.
    Elems:
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
def WriteVTK(Output,Nodes,Elems,grad=0):

    nEle=str(np.shape(Elems)[0])
    nPoi=str(np.shape(Nodes)[0])
    sizeEle=str(np.size(Elems))
    f = open(Output, 'w')
    f.write('# vtk DataFile Version 3.1 \n')
    f.write('Results\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    
    f.write(('Nodes '+nPoi+' DOUBLE\n'))
    np.savetxt(f,Nodes,fmt='%6.6f')
    f.write(('TRIANGLE_STRIPS '+nEle+' '+sizeEle+' \n'))

    np.savetxt(f,Elems,fmt='%d')
    if grad!=0:
        f.write(('POINT_DATA '+nPoi+' \n'))
        f.write('VECTORS flechas double\n')
        np.savetxt(f,gradient,fmt='%+f %+f %d')
def WriteMSH(Output,Nodes,Elems,physicalEnt='line'):
    """
        This funtion intends to create a file of gmsh format containing
        the characteristics of a given mesh and its inherent object.

        For now this program accepts nodes and elements from  a 1D problem
        and returns the files

        Parameters:
        -----------
                    Output      Name given for the output file
                    Nodes       Nodes
                    Elems   Elements
                    Physical entities   lines (for now)
        Returns:
        ----------
                       Output.msh           Nodes and elements in gmsh format 
                       Output.geo           Physical entities contained
    """

    nEle=str(np.shape(Elems)[0])
    sizeEle=str(np.size(Elems))
    nNodes=Nodes.shape[0]
    # Creates the output file for the mesh
    f = open(Output, 'w')
    # Creates the output file for the geometry by stripping the original name
    # and adding the .geo format
    g = open(Output.strip('.msh')+'.geo','w')
    g.write('// Gmsh proyect\n')
#-----------------Physical entities 1D-----------------------------------
    if Nodes.ndim==1:
        xmin=str(Nodes[0])
        xmax=str(Nodes[Nodes.size-1])
       
        g.write('Point(1)={'+xmin+',0,0,1.0};\n')
        g.write('Point(2)={'+xmax+',0,0,1.0};\n')
        g.write('Line(3)={1,2};\n')
    else:
        print 'dimension not supperoted yet'
    g.close()
#------------------------------------------------------------------------
#------------------------ Headers for .msh file------------------------
    f.write('$MeshFormat\n')
    f.write('2.2 0 8\n')
    f.write('$EndMeshFormat\n')
#------------------------ Nodes for the .msh file------------------------            
    
    f.write('$Nodes\n')
    vec=np.zeros((np.shape(Nodes)[0],1))
    for i in range(0,np.shape(Nodes)[0]):
        vec[i,0]=i+1      # Enumeration of nodes
    if Nodes.ndim==1:
        Nodes=np.array([Nodes]).T
        
        # Conversion for coordinates with y and z
        Nodes=hstack((Nodes,np.zeros((np.shape(Nodes)[0],2)))) 
    nPoi=str(np.shape(Nodes)[0])
    Nodes=hstack((vec,Nodes)) # Enumeration added to the array
    f.write(nPoi+'\n')
    np.savetxt(f,Nodes,fmt='%d %f %f %f')
    f.write('$EndNodes\n')
#------------------------ Elements for the .msh file------------------------  

    f.write('$Elements\n')
    f.write(nEle+'\n')
    nEle=int(nEle)
    vec=np.zeros((nEle,5))
    NodPerEl=Elems.ndim   # For line elements there should be 2 nodes per element
    # Point elements -------------------------------------------------------
    if NodPerEl==2:
        ePoints=np.zeros((2,6))
        ePoints[0,:]=[1,15,2,0,1,1]
        ePoints[1,:]=[2,15,2,0,2,nNodes]
        np.savetxt(f,ePoints,fmt='%d')
        
    #-----------------------------------------------------------------------

    
    # Line elements 1D -----------------------------------------------------
    
    
    if NodPerEl==2:
        for i in range(0,nEle):
            vec[i,0]=i+3
            
            vec[i,1:5]=[1,2,0,3]
##        Elems =Elems[:,1:5]
    #-----------------------------------------------------------------------

    # Line elements 2D -----------------------------------------------------        
    else:
        for i in range(0,nEle):
            vec[i,0]=i+1
            vec[i,1:5]=[2,2,0,6]
        Elems =Elems[:,1:5]
    #-----------------------------------------------------------------------
        
    Elems=hstack((vec,Elems))

    
   
    if NodPerEl==2:
        
        np.savetxt(f,Elems,fmt='%d %d %d %d %d %d %d')
    else:
        np.savetxt(f,Elems,fmt='%d %d %d %d %d %d %d %d')
    f.write('$EndElements\n')
def WriteSolverInput(Output,Dimension=1,BCType='Dir',parameter=[],Eq='Schro',\
                     Type='Stationary',AnalisisParam=['y','y',4,4,101]):

    f=open(Output,'r+')
    line=f.readline()
    while '$Solver' not in line:
        here=f.tell()
        line=f.readline()
        
        if '$Solver' in line:
            f.seek(here)
            f.truncate()
            break
        
        if line=='':
            break
    f.write('$Solver input\n')
    Dimension =str(Dimension)
    f.write(Dimension +'\n')
    f.write(BCType +'\n')
    if parameter !=[]:
        np.savetxt(f,parameter,fmt='%f')
    f.write(Eq +'\n')
    f.write(Type +'\n')
    b=str(AnalisisParam)
    b=b.replace("'","")
    f.write(b+'\n')
    f.write('$End Solver input\n')

    
    
        
