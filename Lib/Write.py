#! /usr/bin/python
## module Write
# -*- coding: utf8 -*-
"""
	This module contains functions to write mesh files (pre-defined file formats).
"""
__all__=['WriterVTK','WriteMSH','WriteSolverInput']
__author__="Edward Y. Villegas and Santiago Echeverri"

from numscan import numscan
import numpy as np
from numpy import linspace as lin
from numpy import hstack
def WriterVTK(filename,title,SET,points,cells,data):
    
    """
            WriterVTK write VTK ASCII files from matrixes of point coords and list of element nodes whit their data vector/matrix.
            Parameters:
            -----------
            filename: String which contain filename and path of the output file.
            title:    Title of data (256 characters).
            SET:      Geometry/topology. STRUCTURED_POINTS STRUCTURED_GRID UNSTRUCTURED_GRID POLYDATA RECTILINEAR_GRID FIELD.
                      POLYDATA by default (other still not implemented).
            points:   Matrix of point coords.
            cells:    Matrix of list of element nodes.
            data:     Matrixes of data with arguments. Each set of data has 3 arguments in a list (first add point_data and then cell_data)
                      * datatype:   SCALARS, VECTORS, NORMALS, TENSORS (from matrix future)
                      * dataname:   List of string names of data (default option future)
                      * datamatrix: list of matrices  with data
            Returns:
            --------

    """
    fid = open(filename,'w')
    fid.write('# vtk DataFile Version 2.0\n')
    if title == '':
	title = 'Data'
    if SET == '':
	SET = 'POLYDATA'
    fid.write(title+'\n')
    fid.write('ASCII\n')
    fid.write('DATASET '+SET+'\n')
    n = points.shape[0] # number-of-points
    datatype = 'double' # Future datatype will be extracted from ndarray.dtype
    fid.write('POINTS '+str(n)+' '+datatype+'\n')
    np.savetxt(fid,points,fmt='%6.6f')
    m = cells.shape # number-of-elms and number-of-nodes by elm
    fid.write('POLYGONS '+str(m[0])+' '+str(m[0]*(m[1]+1))+'\n') # elm, elm*(nodbyelm+1) +1 is because include number-of-nodes of polygon
    count = 0
    while count < m[0]:
    	new_elm = np.array(m[1],dtype=int)
	new_elm = np.hstack([new_elm, cells[count,:]])
	if count:
		cellsvtk = np.vstack([cellsvtk,new_elm.copy()])
	else:
		cellsvtk = new_elm.copy()
	count = count + 1
    np.savetxt(fid,cellsvtk,fmt='%d')
    ndata = len(data)/3 # sets of data to visualize

    
    
    if not ndata:
	print len(data),"arguments, but you need to put 3 arguments by set of data."
	return
    count = 0
    point_data = 0
    cell_data = 0

    # From this point it writes for each  type of data  SCALARS, VECTORS, NORMALS, TENSORS

   
    while count < ndata:
        
        
        
        if len(data[1])!=len(data[2]):
            print 'Wrong labeling, or... something else'
            break
        else:   # This else is made to avoid wrong labels

            DataPerType=len(data[2])  # Is the number of data for each type

            #------ Write over each element on a given category----------------------
            k=0
            while k<DataPerType:
            
                print k,DataPerType
                p = data[3*count+2][k].shape
                
                #---------------------------Headers---------------------------------
                if (p[0] == n) and (not point_data):
            
                    fid.write("POINT_DATA "+str(n)+"\n")
                    point_data = 1
                elif (p[0]==m[0]) and (not cell_data):
                    fid.write("CELL_DATA "+str(m[0])+"\n")
                elif (not point_data) and (not cell_data):
                    print "Data Matrix", count,"does not match size with nodes neither elements"
                    count = count + 1
                    continue     
         
                
                #----------------- When we have many sets: ---------------------------

                fid.write("\\\ The following number tells the amount of Scalar vectors to read")
                fid.write('\n'+'\\\  '+str(p[1])+'\n')
                
                count2=0
                while count2<p[1]:
                    
                    if data[3*count+1] == '':
                        data[3*count+1] = 'Data '+str(count+1)+str(count2)
                    TYPE = 'double'
                    fid.write(data[3*count]+' '+data[3*count+1][k]+'_'+str(count2)+' '+TYPE+'\n')
                    fid.write("LOOKUP_TABLE defaul\n")

                    
                    np.savetxt(fid,data[3*count+2][k][:,count2],fmt='%6.6f')
                    count2=count2+1
                k=k+1
        count = count + 1
    fid.close()
    return 0

def WriteMSH(Output,Nodes,Elems,physicalEnt='line'):
    """
        This funtion intends to create a file of gmsh format containing
        the characteristics of a given mesh and its inherent object.

        For now this program accepts nodes and elements from  a 1D problem
        and returns the files

        Parameters:
        -----------
        Output:      Name given for the output file
        Nodes:       Nodes
        Elems:	     Elements
        Physical:    entities   lines (for now)

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

    """
        Write the solver input into a file of gmsh ASCII format V2.2.
        This input is to be read by the Solver module.

        Parameters:
        -----------
        Output:      String with the name of the file that contains the
                     information regarding the geometry, mesh, border
                     conditions, and other parameters necessary for the
                     solution of the problem.

                     This file is a modified Gmsh output file with extension
                     .msh

        Dimension:  int parameter that tells the program wether to solve for a
                    1D problem or a 2D problem (not supported yet)
                    
        parameter:  Is an array that describes the potential actuating over the
                    the elements of the domain given by Elems. For each element in
                    Elems there is an associated potential value on the same
                    position in the array parameter.

                    The potential in Scroedinger equation defines the specific
                    nature of the problem to be solved. For more details on how
                    to define a potential and what does it mean pleas read the
                    documentation of the Potential1D function in the module PrePro.

       

        BCType:     String parameter for the selection of a border condition
                    that can be either:

                        'Dir'   For the Dirichlet border condition
                                (Infinite potential well).

                        'Bloch' For the periodic formulation of the problem.
                                (Electron in a periodic material )

        Type:       String that tells wether to solve the stationary version of
                    the equation or another not yet suported.

                    'Stationary'   

        AnalisisParam:   Array that contains the information regarding the number
                         of solutions to be computed and wether to save the values
                         or not.

                        AnalisisParam[0]:  String  answer to the question
                                                   save  Eigen Values?
                        AnalisisParam[1]:  String  answer to the question
                                                   save  Eigen Vectors?
                        AnalisisParam[2]:  Integer  number of Eigen Values to save
                        AnalisisParam[3]:  Integer  number of Eigen Vectors to save

        Returns:
	--------
        
	
	Last modification: date 25/10/2011
    """
    
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
    np.savetxt(f,parameter,fmt='%f')
    f.write(Eq +'\n')
    f.write(Type +'\n')
    b=str(AnalisisParam)
    b=b.replace("'","")
    f.write(b+'\n')
    f.write('$End Solver input\n')

    
    
        
