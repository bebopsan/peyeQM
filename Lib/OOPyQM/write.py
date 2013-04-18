#! /usr/bin/python
## module Write
# -*- coding: utf8 -*-
"""
	This module contains functions to write mesh files (pre-defined file formats).
"""
__all__=['write_vtk','write_msh','write_solver_input']
__author__="Edward Y. Villegas and Santiago Echeverri"

import numpy as np
from numpy import linspace as lin
from numpy import hstack

def write_vtk(filename, title, SET, points, cells, data):
    
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
    fid.write(title+'\n')
    fid.write('ASCII\n')
    if SET == '':
	SET = 'POLYDATA'

    fid.write('DATASET '+SET+'\n')
    n = points.shape[0] # number-of-points
    datatype = 'double' # Future datatype will be extracted from ndarray.dtype
    fid.write('POINTS '+str(n)+' '+datatype+'\n')
    np.savetxt(fid, points, fmt = '%6.6f')
    m = cells.shape # number-of-elms and number-of-nodes by elm
    if SET == 'POLYDATA':
        fid.write('POLYGONS '+str(m[0])+' '+str(m[0]*(m[1]+1))+'\n') 
    # elm, elm*(nodbyelm+1) +1 is because include number-of-nodes of polygon
    elif SET == 'UNSTRUCTURED_GRID':
        fid.write('CELLS '+str(m[0])+' '+str(m[0]*(m[1]+1))+'\n') 
    else:
        print 'Other types of DATASETS have not been implemented'
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
    if SET == 'UNSTRUCTURED_GRID':
        fid.write('CELL_TYPES '+str(m[0])+'\n')
        cell_ids = np.ones(m[0])
        if m[1] == 8:
            element_tag = 23
        else:
            print 'just cuad Quads for now, have some patience'
        for i in range(m[0]):
            cell_ids[i] = element_tag*cell_ids[i]
        np.savetxt(fid, cell_ids, fmt = '%d')
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

            DataPerType = len(data[2])  # Is the number of data for each type

            #------ Write over each element on a given category----------------------
            k=0
            while k < DataPerType:
            
                p = data[3*count+2][k].shape
                #---------------------------Headers---------------------------------
                if (p[0] == n) and (not point_data):
            
                    fid.write("POINT_DATA "+str(n)+"\n")
                    point_data = 1
                elif (p[0] == m[0]) and (not cell_data):
                    fid.write("CELL_DATA "+str(m[0])+"\n")
                elif (not point_data) and (not cell_data):
                    print "Data Matrix", count, "does not match size with nodes neither elements"
                    count = count + 1
                    continue     
         
                
                #----------------- When we have many sets: ---------------------------

#                fid.write("\\\ The following number tells the amount of Scalar vectors to read")
#                fid.write('\n'+'\\\  '+str(p[1])+'\n')                
                count2 = 0
                if data[0] == 'SCALARS':
                    if data[3*count+2][k].ndim == 1:
                        if data[3*count+1] == '':
                            data[3*count+1] = 'Data '+str(count+1)+str(count2)
                        TYPE = 'double'
                        fid.write(data[3*count]+' '+data[3*count+1][k]+ \
                                  '_'+str(count2)+' '+TYPE+'\n')
                        fid.write("LOOKUP_TABLE defaul\n")
                        np.savetxt(fid, data[3*count+2][k], \
                                    fmt = '%6.6f')
                        count2 = 1
                    else:
                        while count2 < p[1]:
                            
                            if data[3*count+1] == '':
                                data[3*count+1] = 'Data '+str(count+1)+str(count2)
                            TYPE = 'double'
                            fid.write(data[3*count]+' '+data[3*count+1][k]+ \
                                      '_'+str(count2)+' '+TYPE+'\n')
                            fid.write("LOOKUP_TABLE defaul\n")
                            np.savetxt(fid, data[3*count+2][k], \
                                        fmt = '%6.6f')
                            count2 = count2+1
                elif data[0] == 'VECTORS':
                    if data[3*count+2][k].ndim == 1:
                        if data[3*count+1] == '':
                            data[3*count+1] = 'Data '+str(count+1)+str(count2)
                        TYPE = 'double'
                        fid.write(data[3*count]+' '+data[3*count+1][k]+ \
                                  '_'+str(count2)+' '+TYPE+'\n')
                        fid.write("LOOKUP_TABLE defaul\n")
                        np.savetxt(fid, data[3*count+2][k][:, count2], \
                                    fmt = '%6.6f')
                        count2 = 1
                    else:
                        if data[3*count+1] == '':
                            data[3*count+1] = 'Data '+str(count+1)
                        TYPE = 'double'
                        fid.write(data[3*count]+' '+data[3*count+1][k]+ \
                                  ' '+TYPE+'\n')
                        #fid.write("LOOKUP_TABLE defaul\n")
    
                        
                        np.savetxt(fid, data[3*count+2][k], \
                                    fmt = '%+6.6f %+6.6f %d')
                        
                        
                else:
                    print 'Not a supported type  yet'
                k = k+1
        count = count+1  
    fid.close()
    return 0

def write_msh(output, nodes, elements, physicalEnt = 'line'):
    """
        This funtion intends to create a file of gmsh format containing
        the characteristics of a given mesh and its inherent object.

        For now this program accepts nodes and elements from  a 1D problem
        and returns the files

        Parameters:
        -----------
        output:      Name given for the output file
        nodes:       nodes
        elements:	     Elements
        Physical:    entities   lines (for now)

        Returns:
        ----------
        output.msh           Nodes and elements in gmsh format 
        output.geo           Physical entities contained
    """

    n_elements = str(np.shape(elements)[0])
    size_elements = str(np.size(elements))
    n_nodes = nodes.shape[0]
    # Creates the output file for the mesh
    f = open(output, 'w')
    # Creates the output file for the geometry by stripping the original name
    # and adding the .geo format
    g = open(output.strip('.msh')+'.geo','w')
    g.write('// Gmsh proyect\n')
#-----------------Physical entities 1D-----------------------------------
    if nodes.ndim == 1:
        xmin = str(nodes[0])
        xmax = str(nodes[nodes.size-1])
       
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
#------------------------ nodes for the .msh file------------------------            
    
    f.write('$Nodes\n')
    vec = np.zeros((np.shape(nodes)[0], 1))
    for i in range(0,np.shape(nodes)[0]):
        vec[i,0] = i+1      # Enumeration of nodes
    if nodes.ndim==1:
        nodes = np.array([nodes]).T
        
        # Conversion for coordinates with y and z
        nodes = hstack((nodes,np.zeros((np.shape(nodes)[0], 2)))) 
    nPoi = str(np.shape(nodes)[0])
    nodes = hstack((vec, nodes)) # Enumeration added to the array
    f.write(nPoi+'\n')
    np.savetxt(f, nodes, fmt = '%d %f %f %f')
    f.write('$EndNodes\n')
#------------------------ Elements for the .msh file------------------------  

    f.write('$Elements\n')
    f.write(n_elements+'\n')
    n_elements=int(n_elements)
    vec=np.zeros((n_elements,5))
    NodPerEl=elements.ndim   # For line elements there should be 2 nodes per element
    # Point elements -------------------------------------------------------
    if NodPerEl==2:
        ePoints=np.zeros((2,6))
        ePoints[0,:]=[1,15,2,0,1,1]
        ePoints[1,:]=[2,15,2,0,2,n_nodes]
        np.savetxt(f,ePoints,fmt='%d')
        
    #-----------------------------------------------------------------------

    
    # Line elements 1D -----------------------------------------------------
    
    
    if NodPerEl==2:
        for i in range(0,n_elements):
            vec[i,0]=i+3
            
            vec[i,1:5]=[1,2,0,3]
##        elements =elements[:,1:5]
    #-----------------------------------------------------------------------

    # Line elements 2D -----------------------------------------------------        
    else:
        for i in range(0,n_elements):
            vec[i, 0] = i+1
            vec[i, 1:5] = [2,2,0,6]
        elements = elements[:, 1:5]
    #-----------------------------------------------------------------------
        
    elements = hstack((vec, elements))


   
    if NodPerEl==2:
        
        np.savetxt(f, elements, fmt = '%d %d %d %d %d %d %d')
    else:
        np.savetxt(f, elements, fmt = '%d %d %d %d %d %d %d %d')
    f.write('$EndElements\n')
    
def write_solver_input(output, dimension = 1, bc_type = 'Dir', parameter = [],\
                       eq = 'Schro', sol_type = 'Stationary', \
                       analysis_param = ['y', 'y', 4, 4, 20, 20, 2], \
                       bc_filename = ''):
    """
    Write the solver input into a file of gmsh ASCII format V2.2.
    This input is to be read by the Solver module.

    Parameters:
    -----------
    output:      String with the name of the file that contains the
                 information regarding the geometry, mesh, border
                 conditions, and other parameters necessary for the
                 solution of the problem.

                 This file is a modified Gmsh output file with extension
                 .msh

    dimension:  int parameter that tells the program wether to solve for a
                1D problem or a 2D problem (not supported yet)
                
    parameter:  Is an array that describes the potential actuating over the
                the elements of the domain given by Elems. For each element in
                Elems there is an associated potential value on the same
                position in the array parameter.

                The potential in Scroedinger equation defines the specific
                nature of the problem to be solved. For more details on how
                to define a potential and what does it mean pleas read the
                documentation of the Potential1D function in the module PrePro.

   

    bc_type:     String parameter for the selection of a border condition
                that can be either:

                    'Dir'   For the Dirichlet border condition
                            (Infinite potential well).

                    'Bloch' For the periodic formulation of the problem.
                            (Electron in a periodic material )
                                
    sol_type:       String that tells wether to solve the stationary version of
                the equation or another not yet suported.

                'Stationary'   

    analysis_Param:  Array that contains the information regarding the number
                     of solutions to be computed and wether to save the values
                     or not.

                    analysis_param[0]:  String  answer to the question
                                               save  Eigen Values?
                    analysis_param[1]:  String  answer to the question
                                               save  Eigen Vectors?
                    analysis_param[2]:  Integer  number of Eigen Values to save
                    analysis_param[3]:  Integer  number of Eigen Vectors to save
                    analysis_param[4]:  Integer number of wave numbers 
                                        in x to sweep
                    analysis_param[5]:  Integer number of wave numbers 
                                        in y to sweep
                    analysis_param[6]:  biggest value of k. it may be the lenght
                                        of the dommain
                                        
    bc_filename:    string that tells where to look for the boundary 
                    conditions            
    Returns:
    --------
    
	
	Last modification: date 25/10/2011
    """
    
    f = open(output, 'r+')
    line = f.readline()
    while '$Solver' not in line:
        here = f.tell()
        line = f.readline()
        
        if '$Solver' in line:
            f.seek(here)
            f.truncate()
            break
        
        if line == '':
            break
    f.write('$Solver input\n')
    dimension  = str(dimension)
    f.write(dimension +'\n')
    f.write(bc_type +'\n')
    if parameter == []:
        f.write('\n')
    else:
        np.savetxt(f, parameter, fmt = '%f')
    f.write(eq +'\n')
    f.write(sol_type +'\n')
    b = str(analysis_param)
    b = b.replace("'", "")
    f.write(b+'\n')
    f.write(bc_filename+'\n')
    f.write('$End Solver input\n')

    
    
        
