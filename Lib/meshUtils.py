#! /usr/bin/python
## module meshUtils
# -*- coding: utf-8 -*- 
"""

    This module contains functions for the generation of structured meshes
    in 1D and 2D .
    
"""

__all__=['mesh1D','logmesh1D','meshtr2D','logmeshtr2D','linlogmeshtr2D','meshPlot']

__author__="Nicolas Guarin Z."

from numpy import *
import matplotlib.pyplot as plt

# --------------------- 1D meshes --------------------- #

def mesh1D(xmin,xmax,npoints):
    """

        Generate a 1D mesh where the points ares equally spaced.

        Parameters:
        -----------
        xmin:    coordinate of the beginning of the line segment
        xmax:    coordinate of the end of the line segment
        npoints: number of subdivisions


        Returns:
        --------
        coords:  numpy array like vector of the discretized domain with lenght N
        elems:  numpy array like vector of the relations between nodes.
    
    
        Raises:
	-------
    
   
        Last modification: date 21/10/2011
    
    """
    coords = linspace(xmin,xmax,npoints,endpoint=True)
    elems = zeros( (npoints-1,2) ,dtype=int )
    for i in range(0,npoints-1):
        elems[i,0] = i
        elems[i,1] = i+1
    return coords, elems
    
    


def logmesh1D(xmin,xmax,npoints,base,option):
    """

        Generate a 1D mesh where the points are not equally spaced....

        Parameters:
        -----------
        xmin:    coordinate of the beginning of the line segment
        xmax:    coordinate of the end of the line segment
        npoints: number of subdivisions
        base:    tells the base for the logarithm used to distribute coordinates.
                 for base --> 1 like 1.0001 the logmesh tends to a equally
                 distributed mesh
        option:  Is a string variable that will tell the kind
                of Mesh to be made.
		"L":     Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values on
                            the left.
		"R":    Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values
                            on the right side.
		"C":   Creates a 1D domain with values spaced evenly
                            on a log scale with more density of values on
                            the middle.        

        Returns:
        --------
        coords:  numpy array like vector of the discretized domain with lenght N
        elems:  numpy array like vector of the relations between nodes.
    
    
        Raises:
        -------
        Error if the user gives a wrong 'option' argument
        Numerical error if base is chosen equal or less than 1
   
        Last modification: date 21/10/2011
    
    """    
    if(option=="C"):
        if(mod(npoints,2)==0):
            coords = zeros((npoints+1))
        else:
            coords = zeros((npoints))
        xaux = logspace(0, 1, npoints/2+1, endpoint=True,base=base)
        xaux = (xaux-1.0)/(base-1.0)
        for i in range(0,npoints/2+1):
            coords[i] = -xaux[npoints/2-i]
            coords[i+npoints/2] = xaux[i]
        coords = coords*(xmax-xmin)/2.0 + (xmax+xmin)/2.0
    elif(option=="L"):
        coords = zeros((npoints))
        xaux = logspace(0, 1, npoints, endpoint=True,base=base)
        coords = (xaux-1.0)*(xmax-xmin)/(base-1.0) + xmin  
    elif(option=="R"):
        coords = zeros((npoints))
        xaux = logspace(0, 1, npoints, endpoint=True,base=base)
        for i in range(0,npoints):
            coords[i] = -xaux[-i-1]+base
        coords = coords*(xmax-xmin)/(base-1.0)+xmin
    else:
        print "Not known value of parameter 'option' ",option
    elems = zeros( (npoints-1,2) ,dtype=int )
    for i in range(0,npoints-1):
        elems[i,0] = i
        elems[i,1] = i+1
    return coords, elems



# --------------------- 2D meshes --------------------- #


def meshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints):
    """

        Generate a 2D mesh where the points are equally spaced.

        Parameters:
        -----------
        xmin:  initial value of the rectangular domain over x axis
        xmax:  final value of the rectangular domain over x axis
        ymin:  initial value of the rectangular domain over y axis
        ymax:  final value of the rectangular domain over y axis
        nxpoints:     Number of divisions over x
        nypoints:     Number of divisions over y


        Returns:
        --------
        coords:  numpy array like matrix of the discretized domain with shape Nx Ny
        elems:   numpy array like matrix of the relations between nodes.
    
    
        Raises:
        -------
    
        Last modification: date 21/10/2011
    
    """    
    coordx, elem = mesh1D(xmin,xmax,nxpoints)
    coordy, elem = mesh1D(ymin,ymax,nypoints)
    npoints = nxpoints*nypoints
    coords = zeros ( (npoints,2),dtype=float)
    cont = 0
    for j in range(0,nypoints):
        for i in range(0,nxpoints):
            coords[cont,0] = coordx[i]
            coords[cont,1] = coordy[j]
            cont = cont+1
    nelems = 2*(nxpoints-1)*(nypoints-1)
    elems = zeros ( (nelems,3) )
    cont = 0
    for j in range(0,nypoints-1):
        for i in range(0,nxpoints-1):
            elems[cont,0] = j*nxpoints+i
            elems[cont,1] = j*nxpoints+i+1
            elems[cont,2] = nxpoints + j*nxpoints+i+1
            elems[cont+1,0] = j*nxpoints+i
            elems[cont+1,1] = nxpoints + j*nxpoints+i+1
            elems[cont+1,2] = nxpoints + j*nxpoints+i
            cont = cont+2
    return coords, elems



def logmeshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints,optionx,optiony,basex,basey):
    """

        Generate a 2D mesh where the points are logarithmally spaced in both, x and
        y cordinates.

        Parameters:
        -----------
        xmin:  initial value of the rectangular domain over x axis
        xmax:  final value of the rectangular domain over x axis
        ymin:  initial value of the rectangular domain over y axis
        ymax:  final value of the rectangular domain over y axis
        nxpoints:     Number of divisions over x
        nypoints:     Number of divisions over y
        basex: the same base parameter but for the x coordinates of the rectangle
        basey: the same base parameter but for the y coordinates of the rectangle
        optionx:  Lets the user decide the orientation of the log distribution
                  for the custom mesh over the x axis. thi parameter can be
                  L for more points on the left
                  R for more points on the right
                  C for more points on the center
         optiony: Lets the user decide the orientation of the log distribution
                  for the custom mesh over the y axis. thi parameter can be
                  L for more points on the top
                  R for more points on the right
                  C for more points on the bottom


        Returns:
        --------
        coords:  numpy array like matrix of the discretized domain with shape Nx Ny
        elems:   numpy array like matrix of the relations between nodes.
    
    
        Raises:
        -------
        Error if the user gives a wrong 'option' argument
        Numerical error if base is chosen equal or less than 1
   
        Last modification: date 21/10/2011
    
    """    
    coordx, elem = logmesh1D(xmin,xmax,nxpoints,basex,optionx)
    coordy, elem = logmesh1D(ymin,ymax,nypoints,basey,optiony)
    npoints = nxpoints*nypoints
    coords = zeros ( (npoints,2),dtype=float)
    cont = 0
    for j in range(0,nypoints):
        for i in range(0,nxpoints):
            coords[cont,0] = coordx[i]
            coords[cont,1] = coordy[j]
            cont = cont+1
    nelems = 2*(nxpoints-1)*(nypoints-1)
    elems = zeros ( (nelems,3) )
    cont = 0
    for j in range(0,nypoints-1):
        for i in range(0,nxpoints-1):
            elems[cont,0] = j*nxpoints+i
            elems[cont,1] = j*nxpoints+i+1
            elems[cont,2] = nxpoints + j*nxpoints+i+1
            elems[cont+1,0] = j*nxpoints+i
            elems[cont+1,1] = nxpoints + j*nxpoints+i+1
            elems[cont+1,2] = nxpoints + j*nxpoints+i
            cont = cont+2
    return coords, elems




def linlogmeshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints,
                   optiony,basey):
    """

        Generate a 2D mesh where the points are equally spaced in x direction
        and logarithmally spaced in y direction.

        Parameters:
        -----------
        xmin:  initial value of the rectangular domain over x axis
        xmax:  final value of the rectangular domain over x axis
        ymin:  initial value of the rectangular domain over y axis
        ymax:  final value of the rectangular domain over y axis
        nxpoints:     Number of divisions over x
        nypoints:     Number of divisions over y
        basey: the same base parameter but for the y coordinates of the rectangle
         optiony: Lets the user decide the orientation of the log distribution
                  for the custom mesh over the y axis. thi parameter can be
                  L for more points on the top
                  R for more points on the right
                  C for more points on the bottom


        Returns:
        --------
        coords:  numpy array like matrix of the discretized domain with shape Nx Ny
        elems:   numpy array like matrix of the relations between nodes.
    
    
        Raises:
        -------
        Error if the user gives a wrong 'option' argument
        Numerical error if base is chosen equal or less than 1
    
   
        Last modification: date 21/10/2011
    
    """    
    coordx, elem = mesh1D(xmin,xmax,nxpoints)
    coordy, elem = logmesh1D(ymin,ymax,nypoints,basey,optiony)
    npoints = nxpoints*nypoints
    coords = zeros ( (npoints,2),dtype=float)
    cont = 0
    for j in range(0,nypoints):
        for i in range(0,nxpoints):
            coords[cont,0] = coordx[i]
            coords[cont,1] = coordy[j]
            cont = cont+1
    nelems = 2*(nxpoints-1)*(nypoints-1)
    elems = zeros ( (nelems,3) )
    cont = 0
    for j in range(0,nypoints-1):
        for i in range(0,nxpoints-1):
            elems[cont,0] = j*nxpoints+i
            elems[cont,1] = j*nxpoints+i+1
            elems[cont,2] = nxpoints + j*nxpoints+i+1
            elems[cont+1,0] = j*nxpoints+i
            elems[cont+1,1] = nxpoints + j*nxpoints+i+1
            elems[cont+1,2] = nxpoints + j*nxpoints+i
            cont = cont+2
    return coords, elems


# --------------------- Other functionalities --------------------- #



def meshPlot(coords,elems,facecolor,coordnum,elemnum):
    """
        Plot mesh nodes (points) and elements (conectivities between nodes).

        Parameters:
        -----------
        coords:    numpy array like matrix of the discretized domain with shape Nx Ny
        elems:     numpy array like matrix of the relations between nodes
        facecolor: element face color in the drawing (just used in 2D meshes)
        coordnum:  (STILL NOT USED) makes visible the nodes numbering
        elemnum:   (STILL NOT USED) makes visible the elements numbering


        Returns:
        --------
    
    
        Raises:
        -------
        Error if 'facecolor' doesn't correspond with a color name
    
   
        Last modification: date 21/10/2011
    
    """
    if(elems.shape[1]==2):
        y = coords*0
        plt.figure()
        plt.plot(coords,y,'-k')
        plt.plot(coords,y,'ob')
        plt.axis('equal')
        plt.show()
    else:
        plt.figure()
        plt.hold(True)
        for i in range(0,elems.shape[0]):
            xcoord = zeros( (size(elems[i,:])))
            ycoord = zeros( (size(elems[i,:])))
            for j in range(0,size(elems[i,:])):
                xcoord[j] = coords[elems[i,j],0]
                ycoord[j] = coords[elems[i,j],1]
            plt.fill(xcoord,ycoord,facecolor=facecolor,alpha=1.0, edgecolor='k')
        plt.plot(coords[:,0],coords[:,1],'ob')
        plt.axis('equal')
        plt.axis('off')
        plt.show()
    
