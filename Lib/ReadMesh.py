## module ReadMesh
"""
	This module contains functions to read mesh files (pre-defined file formats).
"""

__all__=['Readmsh','ReadSolverInput','ReadVTK']
__author__="Edward Y. Villegas and Santiago Echeverri"

from numscan import numscan
import numpy as np
def Readmsh(filename):
    """
            Read a mesh file of gmsh ASCII format V2.2

            Parameters:
            -----------
            filename:  String with filename and path of the file

            Returns:
            --------
            coords:     numpy array like matrix of node coordinates
                        (nNodes,3) where coorsd(x,:)= x,y,z components of the node 

            elm_lines:  numpy array like matrix of physical entities and list of nodes of element
                        elm_lines is returned as (nElm,3) where
                        elm_lines[x,:]=[2,node1,node2 ]

            elm_triangles:  numpy array like matrix of physical entities and list of nodes of element (if triangle elements exist)
                            elm_triangles is returned as (nElm,4) where
                            elm_lines[x,:]=[3,node1,node2,node4]
            Last modification: date 25/10/2011
    """
    # _________________________ Open file_________________________________
    fid = open(filename, 'r')
    stop = 0 # flag to continue (0) or stop (1) reading
    number = 0 # initialize variable (recognize as global)
    
    #_____________________ Search header file________________________________

    while stop == 0:
        line = fid.readline()
        if '$Mesh' in line:
                line = fid.readline() # Gmsh information format
                stop = 1
    #_______________________ Search nodes __________________________________
    stop = 0
    while stop == 0:
        line = fid.readline()
        if '$Nodes' in line:
                line = fid.readline()
                number = int(line) # Number-of-nodes
                stop = 1
    count = 0
    coords = np.zeros((number,3),dtype=float) #coord-x-y-z
    while count < number:
        line = fid.readline()
        values = line.split()
        coords[count][0:] = values[1:] # Omit Node-number
        count = count + 1
    #______________________ Search elements_________________________________
    stop = 0
    while stop == 0:
        line = fid.readline()
        if '$Elements' in line:
                line = fid.readline()
                number = int(line) # Number-of-elm
                stop = 1
    count = 0
    flag_points = 0 # zero if points elements not added
    flag_lines = 0 # zero if line elements not added
    flag_triangles = 0 # zero if triangle elements not added 
    elm_points = np.zeros((1,2),dtype=int)
    elm_lines = np.zeros((1,3),dtype=int) #phys_ent node1-2
    elm_triangles = np.zeros((1,4),dtype=int) #phys_ent node1-2-3
    while count < number:
        line = fid.readline()
        values = line.split()
        number_tags = int(values[2]) # Number-of-tags (tags still not used)
        new_elm = [values[3]]
        new_elm.extend(values[3+number_tags:])
        new_elm = np.array(new_elm,dtype=int)
        
        if values[1] == '15':
                if flag_lines == 1:
                        elm_points = np.vstack((elm_lines,new_elm))
                        
                else:
                        elm_points = new_elm
                        flag_points = 1
                number=number+1
        elif values[1] == '1':
                if flag_lines == 1:
                        elm_lines = np.vstack((elm_lines,new_elm))
                else:
                        elm_lines = new_elm
                        flag_lines = 1
        elif values[1] == '2':
                if flag_triangles == 1:
                        elm_triangles = np.vstack((elm_triangles,new_elm))
                else:
                        elm_triangles = new_elm
                        flag_triangles = 1
        else:
                print "Type", values[1], "in element", values[0], "is not supported element type."
        count = count + 1
    fid.close()
    if not flag_triangles:
        
        return coords, elm_lines
            
    else:
        return coords, elm_lines, elm_triangles



def ReadVTK(filename):
    """
            Read a mesh file of VTK ASCII format 

            Parameters:
            -----------
            filename:  String which contain filename and path of the file
            Returns:
            --------
            coords:        numpy array like matrix of node coordinates and size
                           (nNodes,3) where coorsd(x,:)= x,y,z components of the node 
            
            elm_lines:     numpy array like matrix of physical entities and list of nodes of element
                           elm_lines is returned as (nElm,3) where
                           elm_lines[x,:]=[2,node1,node2 ]
                                
            elm_triangles: numpy array like matrix of physical entities and list of nodes of
                            element (if triangle elements exist)
                            elm_triangles is returned as (nElm,4) where
                            elm_lines[x,:]=[3,node1,node2,node4]
            elm_scalars:    numpy array like matrix of scalar values asociated with the elements
                                        
            Last modification: date 25/10/
    """
    
    # Open file
    fid = open(filename, 'r')
    stop = 0 # flag to continue (0) or stop (1) reading
    number = 0 # initialize variable (recognize as global)
    # _______________________Search header file   ____________________________
    while stop == 0:
        line = fid.readline()
        if '# vtk' in line:
                line = fid.readline() # Gmsh information format
                Title =line
                stop = 1    
    #___________________________ Search nodes __________________________________
    stop = 0
    while stop == 0:
        line = fid.readline()
        if 'POINTS' in line:
                number =int(line.split()[1]) #Number-of-nodes
                stop = 1
    count = 0
    coords = np.zeros((number,3),dtype=float) #coord-x-y-z

    while count < number:
        line = fid.readline()
        values = line.split()
        coords[count,0:] = values[:] # Omit Node-number
        count = count + 1
    # ________________________ Search elements ______________________________
    stop = 0
    while stop == 0:
            line = fid.readline()
            if 'POLYGONS' in line:
                    
                    number =int(line.split()[1]) # Number-of-elm
                    number2 =int(line.split()[2])# elm*(nodbyelm+1)
                    stop = 1
    count = 0
    flag_lines = 0 # zero if line elements not added
    flag_triangles = 0 # zero if triangle elements not added 
    elm_lines = np.zeros((1,3),dtype=int) #phys_ent node1-2
    elm_triangles = np.zeros((1,4),dtype=int) #phys_ent node1-2-3
    
    while count < number:
            
            line = fid.readline()
            
            values = line.split()
##		number_tags = int(values[2]) # Number-of-tags (tags still not used)
            
            new_elm = [values[0]]
            new_elm.extend(values[1:])
            new_elm = np.array(new_elm,dtype=int)
            
            if values[0] == '2':
                    if flag_lines == 1:
                            elm_lines = np.vstack((elm_lines,new_elm))
                    else:
                            elm_lines = new_elm
                            flag_lines = 1
            elif values[0] == '3':
                    if flag_triangles == 1:
                            elm_triangles = np.vstack((elm_triangles,new_elm))
                    else:
                            elm_triangles = new_elm
                            flag_triangles = 1
                    print 'bip'
            else:
                    print "Type", values[0], "in element", count, "is not a supported element type."
            count = count + 1
            
    # _______________________ Search Scalars _________________________________
    stop = 0
    flag_scalars = 0
    while stop == 0:
            line = fid.readline()
            if 'CELL_DATA' in line:
                number =int(line.split()[1]) # Number of scalars
                stop = 1
            else:
                flag_scalars = 1
                
    if flag_scalars ==0:      
        elm_scalars=np.zeros((number,1),dtype=float)

        stop = 0
        while stop == 0:
                line = fid.readline()
                if 'LOOKUP' in line:
                    stop = 1
        
        count = 0
        while count < number:
            line = fid.readline()
            values = line.split()
            elm_scalars[count,:] = values[:]
            count = count + 1
            
     
    fid.close()
    # ----------------------  Return  arrays -_________________________________
    if not flag_triangles and not flag_scalars:
        
        return coords, elm_lines[1:],elm_scalars
    elif not flag_triangles and  flag_scalars:
        
        return coords,elm_lines[1:]
    elif flag_triangles and not flag_scalars:
        
        return coords, elm_lines[1:], elm_triangles[1:],elm_scalars
    else:
        return coords, elm_lines[1:], elm_triangles[1:]
        
    
def ReadSolverInput(filename):
    """
        Reads the solver input from a file of gmsh ASCII format V2.2.
        This input is to be read by the Solver module.
        Parameters:
        -----------
            filename: String which contain filename and path of the output file.
        Returns:
	--------
            SolverInput: List with the following parameters:

            [Dimension,BCType,parameter,Eq,Type,AnalisisParam]

            Dimension:  int parameter that tells the program wether to solve for a
                    1D problem or a 2D problem (not supported yet)

            BCType:     String parameter for the selection of a border condition
                        that can be either:

                            'Dir'   For the Dirichlet border condition
                                    (Infinite potential well).

                            'Bloch' For the periodic formulation of the problem.
                                    (Electron in a periodic material )
                    
            parameter:  Is an array that describes the potential actuating over the
                        the elements of the domain given by Elems. For each element in
                        Elems there is an associated potential value on the same
                        position in the array parameter.

                        The potential in Scroedinger equation defines the specific
                        nature of the problem to be solved. For more details on how
                        to define a potential and what does it mean pleas read the
                        documentation of the Potential1D function in the module PrePro.

           
            Eq:         Can be 'Schroedinger' or another not implemented yet
            

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
            
         Last modification: date 25/10/2011
    """
    f=open(filename,'r')
    line=f.readline()
    while '$Solver' not in line:
            line=f.readline()
            if '$Solver' in line:
                
                Dimension=int(f.readline())
                BCType=f.readline()
                here=f.tell()
                if f.readline ==[]:
                    params=[]
                else:
                    f.seek(here)
                    parameter=np.array(numscan(f))
                    while True:
                        here=f.tell()
                        b=numscan(f)
                        if b==[]:
                            break
                        parameter=np.vstack((parameter,b))
                f.seek(here)
                Eq=f.readline()
                Type=f.readline()
                AnalisisParam=f.readline()
                AnalisisParam=AnalisisParam.strip('[ ]\n')
                AnalisisParam=AnalisisParam.split(',')
                SolverInput=[Dimension,BCType,parameter,Eq,Type,AnalisisParam]
                return SolverInput
                break
            if line=='':
                SolverInput=[]
                print 'Found no input for solver module' 
                return SolverInput
                break
            
        
        
        
    
    
    


























    
        
