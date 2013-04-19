# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 11:01:40 2013

@author: santiago
"""
import os, sys
lib_path = os.path.abspath('..')
sys.path.append(lib_path)
from Classes import Region
import pickle

filename = 'unit_cell.reg'

air_prop = {'mu':1, 'epsilon':1}
air = Region(tag = '30', name = 'Air', material_prop = air_prop)

dielectric_prop = {'mu':1, 'epsilon':10}
dielectric = Region(tag = '17', name = 'Dielectric', material_prop = dielectric_prop)

f = open(filename,'w')
regions = [air]
pickle.dump(regions, f)
f.close()
print 'done'