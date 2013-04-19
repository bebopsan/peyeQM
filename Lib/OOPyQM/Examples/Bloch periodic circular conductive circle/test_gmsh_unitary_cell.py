# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:52:13 2013

@author: Santiago Echeverri Chacón
"""

from gmsh_library import *

point = Point(0,0,0)
frame = Unitary_cell(0, point, 5, 5, transfinite = (2,1), phys_tag = 30)
circle_center = Point(frame.points[-1].id_tag, 2.5,2.5)
frame.add_circular_inclussion(1,circle_center,1)
for s in frame.surfaces[5:]:
    print s
f = open('unit_cell.geo','w')
f.write(frame.__str__())
f.close