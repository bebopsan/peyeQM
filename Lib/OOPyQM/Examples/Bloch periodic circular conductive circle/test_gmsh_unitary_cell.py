# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:52:13 2013

This creates a unitary cell with a circular inclusion placed at the midle.

Pseudo:
- Import stuff
- Instantiate corner point at origin
- Instantiate Unitary cell with the point from before and the characteristic 
  lenght for each side, In this case "a". 
  Aditionally, define physical tags for the regions outside and
  inside of the inclusion, these as well as the transfinite definition will
  act as optional parameters. If no Physical tags are defined gmsh wont 
  include them when the mesh file is saved.
- Define a Point() with placement relative to the empty unitary cell. And 
  tag taken from the last point of list of points inside it. 

@author: Santiago Echeverri Chacón
"""


import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)
from gmsh_library import *
from swampy import Lumpy
lumpy = Lumpy.Lumpy()
lumpy.make_reference()
point = Point(0,0,0)

a = 1
phys_tag_outside = 30
phys_tag_inside = 31

frame = Unitary_cell(0, point, a, a, transfinite = (2,1), phys_tag = phys_tag_outside)
circle_center = Point(frame.points[-1].id_tag, a/2.,a/2.)
frame.add_circular_inclussion(1,circle_center,0.2, phys_tag = phys_tag_inside)
lumpy.object_diagram()
for s in frame.surfaces[5:]:
    print s
f = open('unit_cell_1-02_two_reg.geo','w')
f.write(frame.__str__())
f.close