# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 16:01:07 2013
This script creates a finite lattice of dielectric rods.

Pseudo:

- Reference path to libraries
- import stuff
- create bottom left corner point for grid
- Define lenght for unitary cell
- define relative point of circle center in unitary cell
- // //  tags for physical surfaces as to distinguish dielectric from air
- // //  Sketh of structure for the grid
- create a property bag to to assign into a dictionary. This property bag 
  will be identified with key '1' of the dictionary.
- Instantiate domain as a Grid.
- define boundary of the domain.

@author: santiago
"""

import os, sys
lib_path = os.path.abspath('../..')
sys.path.append(lib_path)
from gmsh_library import * 

blc = Point(0,0,0) #Bottom left corner

a = 2 # unit cell width
cc = Point(0, a/2.,a/2.) # circle center

phys_tag_outside = 100
# phys_tag_inside = 31

pt = [phys_tag_outside, None]#, phys_tag_inside] # physical tags

sketch = [['0','0','0','0','0','0','0','0','0']]
c = Properties_bag('circle', origin = cc, r = 0.1)
e = Properties_bag('empty')
d = {'1': c,'0':e}

domain =  Grid(blc, a, sketch, d, phys_tags = pt, transfinite = (2,1))
domain.define_boundary()

f = open('line.geo','w')
f.write(domain.__str__())
f.close

print 'Done'