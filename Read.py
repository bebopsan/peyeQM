#! /usr/bin/python
## module Read
''' points,triang=ReadVRML(archivo).
    Reads a data file of file with VRML2 format, and
    returns two array that represents points and triangles 
    The input "archivo" is the name of the file to be opened
    (e.g "dosTriangulos.wrl")
'''
from numscan import numscan
import numpy as np
def ReadVRML(archivo):
    #Open the file
    f=open(archivo,'r')
    stop=0
    while stop==0:
        line=f.readline()
        if 'point [' in line:
            stop=1
    
    n=np.zeros(3)
    n=numscan(f)
    points=n
    while stop==1:
        
        back=f.tell()
        line=f.readline()
        f.seek(back)
        n=numscan(f)
        
        if n!=[]:
            points=np.vstack((points,n))  
        if ']'in line:
            stop=0
                      
    while stop==0:
        line=f.readline()
        if 'coordIndex' in line:
            stop=1
    n=numscan(f)
    triangles=n[0:3]
    while stop==1:
        back=f.tell()
        line=f.readline()
        f.seek(back)
        n=numscan(f)
        if n!=[]:
            triangles=np.vstack((triangles,n[0:3]))  
        if ']'in line:
            stop=0

    return points,triangles
    
