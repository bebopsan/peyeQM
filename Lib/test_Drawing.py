#!/usr/bin/python
# -*- coding: utf8 -*-

from drawing import *
from numpy import *

# curve_vect
##x = linspace(-5,5,100,endpoint=True)
##y = exp(-x**2)
##curve_vect(x,y,['x','y','Title'],'D:/try.svg')



### 3D Plots ###

fun = '3.*(1-x)**2*exp(-x**2 - (y+1)**2) -\
10.*(x/5. - x**3 - y**5)*exp(-x**2-y**2)-\
1./3.*exp(-(x+1)**2 - y**2)'

domain = [-3,3,-3,3,10./50.,10./50.]

# surface

surface(domain,fun,['x','y','z'])


# surface_matrix_vect
x = linspace(-3,3,49,endpoint=True)
y = linspace(-3,3,49,endpoint=True)
x,y = meshgrid(x,y)
z = eval(fun)
surf_matrix_vect([x,y],z,['x','y','z'],'D:/surf.svg')


# contour_vect
contour_vect([x,y],z,['x','y','z'],'D:/contour.svg')


