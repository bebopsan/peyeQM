#! /usr/bin/python
## module FIOTools
# -*- coding: utf8 -*-
"""
	This module contains functions to manipulate I/O strings by
	console and read and write files (pre-defined file formats).
"""

__all__=['Readmsh','WriterVTK']
__author__="Edward Y. Villegas"

import numpy as np

def Readmsh(filename):
	"""
		Read a mesh file of gmsh ASCII format V2.2

		Parameters:
        	-----------
		filename:  String which contain filename and path of the file
		Returns:
		--------
		coords:        numpy array like matrix of node coordinates
		elm_lines:     numpy array like matrix of physical entities and list of nodes of element
		elm_triangles: numpy array like matrix of physical entities and list of nodes of element (if triangle elements exist)

		Last modification: date 23/10/
	"""
	# Open file
	fid = open(filename, 'r')
	stop = 0 # flag to continue (0) or stop (1) reading
	number = 0 # initialize variable (recognize as global)
	# Search header file
	while stop == 0:
		line = fid.readline()
		if '$Mesh' in line:
			line = fid.readline() # Gmsh information format
			stop = 1
	# Search nodes
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
	# Search elements
	stop = 0
	while stop == 0:
		line = fid.readline()
		if '$Elements' in line:
			line = fid.readline()
			number = int(line) # Number-of-elm
			stop = 1
	count = 0
	flag_lines = 0 # zero if line elements not added
	flag_triangles = 0 # zero if triangle elements not added 
	elm_lines = np.zeros((1,3),dtype=int) #phys_ent node1-2
	elm_triangles = np.zeros((1,4),dtype=int) #phys_ent node1-2-3
	while count < number:
		line = fid.readline()
		values = line.split()
		number_tags = int(values[2]) # Number-of-tags (tags still not used)
		new_elm = [values[3]]
		new_elm.extend(values[3+number_tags:])
		new_elm = np.array(new_elm,dtype=int)
		if values[1] == '1':
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

def Readgeo(filename):
	"""
		in develop
	"""
	return 0
	
def Writergeo(filename, points, options, addphys):
	"""
		Writergeo write geo files (gmsh input file to meshing) from ordered points matrix.
		Parameters:
		-----------
		filename: String which contain filename and path of the output file.
		points  : numpy array like matrix of ordered points.
		options : Options of physical entity to create
			   * 'p'   : Create physical points
			   * 'l'   : Create physical lines
			   * 's'   : Create physical surface
			   * 'none': Not create phisical entities
		addphys : 
	"""
	fid = open(filename,'w')
	count = 1
	numpoints = points.shape[0]
	while count <= numpoints:
		fid.write('Point('+str(count)+') = {'+str(points[count-1,0])+', '+str(points[count-1,1])+', '+str(points[count-1,2])+'};\n')
		count = count + 1
	count =1
	while count <= numpoints:
		fid.write('Line('+str(count)+') = {')
		if count < points.shape[0]:
			fid.write(str(count)+', '+str(count+1)+'};\n')
		else:
			fid.write(str(count)+', 1};\n')
		count = count + 1
	lines = 1
	fid.write('Line Loop('+str(count)+') = {')
	while lines <= numpoints:
		if lines < numpoints:
			fid.write(str(lines)+', ')
		else:
			fid.write(str(lines)+'};\n')
		lines = lines + 1
	notphys = count + 1
	countphys = notphys + 1
	fid.write('Plane Surface('+str(notphys)+') = {'+str(count)+'};\n')
	count = 1
	while count <= numpoints:
		fid.write('Physical Line('+str(countphys)+') = {'+str(countphys-notphys)+'};\n')
		count = count + 1
		countphys = countphys + 1
	fid.write('Physical Surface('+str(countphys)+') = {'+str(countphys-notphys+1)+'};\n')
	fid.close()
	return 0

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
			  * dataname:   String name of data (default option future)
			  * datamatrix: Matrix of data
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
	shiftnodes = np.ones((1,m[1]),dtype=int)
	while count < m[0]:
		new_elm = np.array(m[1],dtype=int)
		new_elm = np.hstack([new_elm, cells[count,:]-shiftnodes[0,:]])
		if count:
			cellsvtk = np.vstack([cellsvtk,new_elm.copy()])
		else:
			cellsvtk = new_elm.copy()
		count = count + 1
	np.savetxt(fid,cellsvtk,fmt='%d')
	tdata = str(type(data)).split('\'')
	if tdata[1]== 'list':
		ndata = len(data)/3 # sets of data to visualize
	else:
		print "DATA argument should be a list"
		return
	if not ndata:
		print ndata,"sets. You need 3 arguments by set of data."
		return
	count = 0
	point_data = 0
	cell_data = 0
	while count < ndata:
		p = data[3*count+2].shape
		if (p[0] == n) and (not point_data):
			fid.write("POINT_DATA "+str(n)+"\n")
			point_data = 1
		elif (p[0]==m[0]) and (not cell_data):
			fid.write("CELL_DATA "+str(m[0])+"\n")
			cell_data = 1
		elif (not point_data) and (not cell_data):
			print "Data Matrix", count,"does not match size with nodes neither elements"
			count = count + 1
			continue
		if data[3*count+1] == '':
			data[3*count+1] = 'Data '+str(count+1)
		TYPE = 'double'
		fid.write(data[3*count]+' '+data[3*count+1]+' '+TYPE+'\n')
		fid.write("LOOKUP_TABLE defaul\n")
		np.savetxt(fid,data[3*count+2],fmt='%6.6f')
		count = count + 1
	fid.close()
	return 0

def typeMatrix(matrix):
	"""
		This function extract type of numpy array to use in VTK argument of datatype (int, float, double)
		IN DEVELOPMENT
	"""
	TYPE = str(matrix.dtype).split('\'')
	# regular expression to extract word below number 'int32'->'int' , 'float64'->'float' (in develop) 
	return TYPE
