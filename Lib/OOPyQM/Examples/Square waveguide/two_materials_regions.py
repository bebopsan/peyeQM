# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 11:01:40 2013

This script creates dictionaries that represent properties of regions
insidea a certain domain of a FEM simulation. Those dictionaries are
saved as pickle dumps inside a list and are to be read 
when runing the simulation.

Pseudo:
- Import stuff
- Define a filename root

@author: santiago
"""
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)
from Classes import Region
import pickle

#filename = 'unit_cell_1-04_two_reg_guide.reg'
filename = 'square_waveguide_2a.reg'
air_prop = {'mu':1, 'epsilon':1}
air = Region(tag = '8', name = 'Air', material_prop = air_prop)


#air1 = Region(tag = '31', name = 'Air', material_prop = air_prop)
dielectric_prop = {'mu':1, 'epsilon':8.9}
dielectric = Region(tag = '31', name = 'Dielectric', material_prop = dielectric_prop)

f = open(filename,'w')
regions = [air]
pickle.dump(regions, f)
f.close()
print 'done'