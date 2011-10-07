#! /usr/bin/python
## module Integral
"""

    This module contains functions for the generation of structured meshes
    in 1D and 2D.
    
"""

from numpy import *
import matplotlib.pyplot as plt

### LINEAR MESH 1D ###
"""

    Generate a 1D mesh where the points ares equally spaced.

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 1/Oct/2011
    
"""
def mesh1D(xmin,xmax,npoints):
    coords = linspace(xmin,xmax,npoints,endpoint=True)
    elems = zeros( (npoints-1,2) ,dtype=int )
    for i in range(0,npoints-1):
        elems[i,0] = i
        elems[i,1] = i+1
    return coords, elems


### LOGARITHMICAL MESH 1D ###
"""

    Generate a 1D mesh where the points are not equally spaced....

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 1/Oct/2011
    
"""
def logmesh1D(xmin,xmax,npoints,base,option):
    
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


### LINEAR MESH 2D ###
"""

    Generate a 2D mesh where the points are equally spaced.

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 2/Oct/2011
    
"""
def meshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints):
    
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


### LOG MESH 2D ###
"""

    Generate a 2D mesh where the points are not equally spaced.

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 2/Oct/2011
    
"""
def logmeshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints,optionx,optiony,basex,basey):
    
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


### LIN-LOG MESH 2D ###
"""

    Generate a 2D mesh where the points are not equally spaced.

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 2/Oct/2011
    
"""
def linlogmeshtr2D(xmin,xmax,ymin,ymax,nxpoints,nypoints,
                   optiony,basey):
    
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


### MESH PLOT ###
"""

    Plot a mesh....

    Input:


    Output:


    Created by: Nicolas Guarin
    Last modification: 2/Oct/2011
    
"""
def meshPlot(coords,elems,facecolor,coordnum,elemnum):

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
    
