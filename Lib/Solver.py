## module Solver
# -*- coding: utf-8 -*-
"""
Module Solver recieves data from the preprocessing stage and performs
the procesor stage of the Finite Element Method procedure.

It should be able to take a file as input or the matrices and options directly
from the main.
"""

__all__=['Srhoedinger','Dirichlet']
__author__='Santiago Echeverri Chacón'

import Read
import Write

def Schroedinger(File,Nodes=0,Elems=0,Params=[],dim=1,\
                 AnalysisType='stationary',AnalisisParam=[]):
    
