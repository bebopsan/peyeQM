# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 10:13:04 2011

@author: santiago
"""

from read_mesh import read_mesh, read_bc
from vectors import  image_reference_bloch_vectors
from utils import bloch_multiplication
nodes, elements = read_mesh('test1.msh')

lines = elements[1]

bc = read_bc('test2.bc')
bloch = bc[2]

image_reference_bloch_vectors(lines, bloch)