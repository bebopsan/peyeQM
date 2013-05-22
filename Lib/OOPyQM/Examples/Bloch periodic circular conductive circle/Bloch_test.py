# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 18:57:41 2013

@author: santiago
"""
from numpy import array, pi, asarray, set_printoptions
from Classes import Boundaries 
bloch = Boundaries() 
a = pi / 4.
b = pi / 4.
k_x = 1
k_y = 1
nodes = array([[0, 0],[a,0], [a,b], [0,b] ])
K = array([[-1, 0, 1, 0, 0, 0, 0, 0],\
           [ 0,-1, 0, 0, 0, 0, 0, 1],\
           [ 1, 0,-1, 0, 0, 0, 0, 1],\
           [ 0, 0, 0,-1, 0, 1, 0, 1],\
           [ 0, 0, 0, 0,-1, 0, 1, 0],\
           [ 0, 0, 0, 1, 0,-1, 0, 0],\
           [ 0, 0, 0, 0, 1, 0, -1, 0],\
           [ 0, 1, 0, 0, 0, 0, 0, -1]])
K = asarray(K,dtype=complex)
ref_im = array([[0, 2],[1,3],[0, 4],[1,5], [0,6],[1,7]])
print ref_im
K2 = bloch.bloch_multiplication(k_x,k_y,nodes,ref_im,K) 
set_printoptions(precision=1)
print 'K2[0]', K2[0]
print 'K2[0]-K2[0].conjugate().transpose()',K2[0]-K2[0].transpose().conjugate()
set_printoptions(precision=1)
K3 = bloch.bloch_sum(ref_im, K2[0])
print 'K3[0]', K3[0]
print 'K3[0]-K3[0].conjugate().transpose()',K3[0]-K3[0].conjugate().transpose()


