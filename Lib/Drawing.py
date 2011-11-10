#!/usr/bin/python
## module Drawing
# -*- coding: utf8 -*-
"""

    This module contains functions for make easy the plotting work.
    
"""

__all__=['surface','surfaceM', 'contourM','curveV']
__author__="Edward Y. Villegas and Nicolas Guarin."

from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mp
import matplotlib.pyplot as plot
from numpy import *
mp.text.usetex=True # Ever use tex engine to write text

# --------------------- 2D plotting --------------------- #
def curveV(v1,v2,labels,path):
    """

        Plot a cuve....

        Parameters:
        -----------
        v1:     List or array for x data.
        v2:     List or array for y data.
        labels: List with strings for the axes labels.
        path:   Path where the file will be placed.
        


        Returns:
        --------

    
    
        Raises:
        -------

   
        Last modification: date 10/11/2011
    
    """
    plot.plot(v1,v2)
    plot.xlabel(labels[0])
    plot.ylabel(labels[1])
    plot.title(labels[2])
    plot.savefig(path)

# --------------------- 3D Plotting --------------------- #
def surface(domain,fun,labels):
    """

        Plot a surface...

        Parameters:
        -----------
        domain:     
        fun:     List or array for y data.
        labels: List with strings for the axes labels.


        Returns:
        --------

    
    
        Raises:
        -------

   
        Last modification: date 10/11/2011
    
    """
    if not (len(domain)==6):
        if len(domain)==4:
            dx = (domain[1]-domain[0])/20.
            dy = (domain[3]-domain[2])/20.
        else:
            print "Domain should be have 4 or 6 elements,"\
                  "but it have", len(domain)
        return
    else:
        dx = domain[4]*1.
        dy = domain[5]*1.
    xmin = domain[0]*1.
    xmax = domain[1]*1.
    ymin = domain[2]*1.
    ymax = domain[3]*1.
    figure = plot.figure()
    axes = Axes3D(figure)
    axes.set_xlabel(labels[0])
    axes.set_ylabel(labels[1])
    axes.set_zlabel(labels[2])
    axes.set_aspect('equal')
    x = arange(xmin,xmax,dx)
    y = arange(ymin,ymax,dy)
    x,y = meshgrid(x,y)
    z=eval(fun)
    axes.plot_surface(x, y, z, rstride=1, cstride=1, cmap=mp.cm.jet)
    plot.show()


def surfaceM(domain,Values,labels,path):
    """

        Plot a surface ...

        Parameters:
        -----------
        domain:     
        fun:     
        labels: List with strings for the axes labels.
        path:


        Returns:
        --------

    
    
        Raises:
        -------

   
        Last modification: date 10/11/2011
    
    """
    figure = plot.figure()
    axes = Axes3D(figure)
    axes.set_xlabel(labels[0])
    axes.set_ylabel(labels[1])
    axes.set_zlabel(labels[2])
    axes.set_aspect('equal')
    axes.plot_surface(domain[0],domain[1], Values.real, rstride=1, \
                      cstride=1, cmap=mp.cm.jet)
    plot.savefig(path)


def contourM(domain,values,labels,path):
    """

        Plot a surface ...

        Parameters:
        -----------
        domain:
        values:
        fun:     
        labels: List with strings for the axes labels.
        path:


        Returns:
        --------

    
    
        Raises:
        -------

   
        Last modification: date 10/11/2011
    
    """
    figure = plot.figure()
    plot.contourf(domain[0],domain[1],values.real)
    plot.xlabel(labels[0])
    plot.ylabel(labels[1])
    plot.title(labels[2])
    plot.savefig(path)
    

