## module read_mesh
"""
	This module contains functions to read mesh files (pre-defined file formats).
"""

__all__ = ['read_mesh','read_solver_input', 'read_vtk', 'read_bc']
__author__ = "Edward Y. Villegas and Santiago Echeverri"

from numscan import numscan
import numpy as np
def read_mesh(filename):
    """
        Read a mesh file of gmsh ASCII format V2.2

        Parameters:
        -----------
        filename:  String with filename and path of the file

        Returns:
        --------
        coords:     numpy array like matrix of node coordinates
                    (n_nodes,3) where coorsd(x,:)= x,y,z components 
                    of the node 

        elements:   list with the following arrays inside:  
        
            elm_lines:  numpy array like matrix of physical entities and list 
                        of nodes of element
                        elm_lines is returned as (nElm,3) where
                        elm_lines[x,:]=[2,node1,node2 ]
    
            elm_triangles:  numpy array like matrix of physical entities and 
                            list of nodes of element (if triangle elements exist)
                            elm_triangles is returned as (nElm,4) where
                            elm_lines[x,:]=[3,node1,node2,node4]
                            
        Last modification: date 09/11/2011
    """
    assert isinstance(filename, str) 
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
    #------------------------ Flags definition-----------------------------
    flag_points = 0         # zero if points elements not added
    flag_lines = 0          # zero if line elements not added
    flag_triangles = 0      # zero if triangle elements not added 
    
    #------------------------- Array initialization ---------------------    
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
                #number=number+1
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
                print "sol_type", values[1], "in element", values[0], \
                      "is not supported element type."
        count = count + 1
    fid.close()
    if not flag_triangles:
        
        elements = [coords, elm_lines]
            
    else:
        elements = [coords, elm_lines, elm_triangles]
       
    return coords, elements

def read_vtk(filename):
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
                    print "sol_type", values[0], "in element", count, "is not a supported element type."
            count = count + 1
            
    # _______________________ Search Scalars _________________________________
    
    
    
           
    
    flag_scalars = 0
    flag_point_data=0
    flag_cell_data=0
    while flag_scalars==0:
        stop=0
        while stop == 0:
            line = fid.readline()
            
            if 'POINT_DATA' in line and flag_point_data==0 :
                number =int(line.split()[1]) # Number of point data scalars
                flag_point_data=1  
                stop==1
                break
            elif 'CELL_DATA' in line and flag_cell_data==0:
                
                number2 =int(line.split()[1]) # Number of cell data scalars
                print number2
                flag_cell_data=1   
                stop = 1
            elif len(line) == 0:
                print 'If you want to plot scalars, you missed something...'
                stop=1
                flag_scalars = 1
            else:
                flag_scalars = 1
                
                
        line = fid.readline()
               
        if '\\ The' in line:
            line = fid.readline()
            nSol=int(line.split()[1])
        else:
            nSol=1
        
        if flag_scalars ==0:

            # ---------------- From here on are the point data scalars  -----------
            if flag_point_data==1: 
                elm_point_scalars=np.zeros((number,nSol),dtype=float)

                count=0
                while count<nSol:
                    
                    stop = 0
                    while stop == 0:
                            line = fid.readline()
                            if 'LOOKUP' in line:
                                stop = 1
                            elif len(line)==0:
                                stop=1
                    
                    count2 = 0
                    while count2 < number:
                        line = fid.readline()
                        values = line.split()[0]
                        elm_point_scalars[count2,count] = values
                        count2 = count2 + 1
                    count=count+1
                flag_point_data=0
            # ---------------- From here on are the cell data scalars  -----------
            if flag_cell_data==1:
                elm_cell_scalars=np.zeros((number2,nSol),dtype=float)
                count=0
                while count<nSol:
            
                    stop = 0
                   
                    while stop == 0:
                            line = fid.readline()
                            
                            if 'LOOKUP' in line:
                                stop = 1
                            elif len(line)==0:
                                stop=1
                                       
                    count2 = 0
                    
                    
                    while count2 < number2:
                        line = fid.readline()
                        values = line.split()[0]
                        elm_cell_scalars[count2,count] = values
                        count2 = count2 + 1
                    count=count+1
                flag_cell_data=0
          
    elm_point_scalars=np.array(elm_point_scalars)
   
    fid.close()
    # ----------------------  Return  arrays -_________________________________
    if not flag_triangles and  flag_scalars:
        
        return coords, elm_lines[1:],elm_cell_scalars,elm_point_scalars
    elif not flag_triangles and not  flag_scalars:
        
        return coords,elm_lines[1:]
    elif flag_triangles and not flag_scalars:
        
        return coords, elm_lines[1:], elm_triangles[1:],elm_cell_scalars,elm_point_scalars
    else:
        return coords, elm_lines[1:], elm_triangles[1:]
        
    
def read_solver_input(filename):
    """
    Reads the solver input from a file of gmsh ASCII format V2.2.
    This input is to be read by the Solver module.
    Parameters:
    -----------
    filename: String which contain filename and path of the output file.
    Returns:
    --------
    solver_input: List with the following parameters:

    [dimension,bc_type,parameter,eq,sol_type,analysis_param]

    dimension:  int parameter that tells the program wether to solve for a
            1D problem or a 2D problem (not supported yet)

    bc_type:     String parameter for the selection of a border condition
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

   
    eq:         Can be 'Schroedinger' or another not implemented yet
    

    sol_type:       String that tells wether to solve the stationary version of
                the equation or another not yet suported.

                'Stationary'   

    analysis_param:   Array that contains the information regarding the number
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
                        
     Last modification: date 14/11/2011
    """
    f = open(filename,'r')
    line = f.readline()
    while '$Solver' not in line:
            line = f.readline()
            if '$Solver' in line:
                
                dimension = int(f.readline())
                bc_type = f.readline().strip('\n')
                here = f.tell()
                if f.readline == []:
                    parameter = []
                else:
                    f.seek(here)
                    parameter = np.array(numscan(f))
                    while True:
                        here = f.tell()
                        b = numscan(f)
                        if b == []:
                            break
                        parameter = np.vstack((parameter, b))
                f.seek(here)
                eq = f.readline().strip('\n')
                sol_type = f.readline().strip('\n')
                analysis_param = f.readline()
                analysis_param = analysis_param.strip('[ ]\n')
                analysis_param = analysis_param.split(',')
                bc_filename = f.readline().strip('\n')
                solver_input = [dimension, bc_type, parameter, eq, sol_type,  \
                                analysis_param, bc_filename]
                return solver_input
                break
            if line == '':
                solver_input = []
                print 'Found no input for solver module' 
                return solver_input
                break

#==================  Boundary conditions reader ==============
            
def read_bc(filename):
    """
    Reads the boundary conditions given by a .bc file.    
    Where D stands for Dirichlet boundaries, N to Newman, and B for Bloch
    periodic. Physical lines have to be created in order to asociate boundary
    nodes to their corresponding lines.    
    Each boundary condistion consist in three lines like the following example:
    
    8
    B 1
    0 10 1          
    
    For Dirichlet and Newman conditions:    
    The initial number refers to the tag of the physical line, the letter of the 
    second line tells the type of condition, and the third line has the
    number of degrees of freedom and the value for each.
    
    However for Bloch periodic conditions:
    The initial number refers to the tag of the physical line, the first number 
    of the third line holds no meaning, the second number of the third line 
    tells which is the "partner" of the curent physical line, and the third 
    term of the third line tells if it is an "image" or "reference" entity. 
    
    Parameters:
    -----------
     filename: String which contain filename and path of the output file.
     
    Returns:
    --------
    
    bloch_2: Is a list of lists where each of the inner lists represents a pair 
             of image and reference physical line entities. 
    
    Note: This file should not have empty lines.
    
    Last modification: date 27/04/2012
    """
    from numpy import zeros
    assert isinstance(filename, str) 
    fid = open(filename, 'r')
    stop = 0  # flag to continue (0) or stop (1) reading
    #======Create empty lists =========
    dirichlet = []   
    newman = []
    bloch = []
    bloch_flag = 0
    bc_type = ''
    ndf = 0
    while stop == 0:      
        line = fid.readline().split()
        if len(line) == 0:  # Is this line empty? Is this the end of the file?
            stop = 1         
        elif len(line) == 1:    # This should tell that you are over the tag.     
           tag =  int(line[0]) 
           count = 0
           deg_of_fre = []
           
        elif len(line) == 2:  # This is the second line of a boundary condition.
            bc_type = line[0]    # Name of the border condition type 
            ndf = int(line[1])        # Number of degrees of freedom
            
        elif len(line) == 3:
            deg_of_fre.append(np.array(line[1:], dtype = 'float'))
            count = count + 1
        else:
           print 'error'
        #=============== Dirichlet conditions ===============================   
        if bc_type == 'D' and count == ndf and stop == 0:
            
            dirichlet.append([tag, deg_of_fre])
        #================ Newman conditions ================================    
        elif bc_type == 'N' and count == ndf and stop == 0:
            
            newman.append([tag, deg_of_fre])
        #================= Bloch periodicity conditions=====================
        elif bc_type == 'B' and count == ndf and stop == 0:
            bloch.append([tag, int(deg_of_fre[0][0]), int(deg_of_fre[0][1])])
            bloch_flag = 1
            
    # === This block of code rearranges the bloch  conditions into an array ===
    if bloch_flag == 1: 
        
        n_bloch = len(bloch)
        bloch_2 = zeros((n_bloch/2, 4), dtype = int) 
        print "bloch: ", bloch
        assert n_bloch % 2 == 0   # Make sure that bloch conditions are paired
        it_bloch = range(n_bloch)    # Iteration (it) list for bloch conditions   
        for i in it_bloch:
            bloch_2[i, 0] = int(bloch[i][0])
            bloch_2[i, 2] = int(bloch[i][2])
            for j in it_bloch:
                if bloch[i][1] == bloch[j][0]:
                    bloch_2[i, 1] = int(bloch[j][0])
                    bloch_2[i, 3] = int(bloch[j][2])
                    it_bloch.pop(j)
        
        print "bloch_2: ", bloch_2       
        return dirichlet, newman, bloch_2
    else:
        return dirichlet, newman
        