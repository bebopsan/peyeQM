## module ReadMesh
''' points,triang=ReadMSH(mesh,geom).
    Reads two files with GMSH format, a file with info about the
    mesh and a another one describing the geomety of the problem.
    
    Returns two arrays that represent points, lines and triangles.
    The input "archivo" is the name of the file to be opened
    (e.g "archivo.msh" and "geom.geo")
'''

from numscan import numscan
import numpy as np
def ReadMesh(mesh,geometry=0):
    #Open the file
    f=open(mesh,'r')
    stop=0
    while stop==0:
        line=f.readline()
        if '$Nodes' in line:
            stop=1
    n=numscan(f)
    nNodes=n[0]
    back=f.tell()
    n=np.array(numscan(f)).size
    if n==4:
        f.seek(back)
        Nodes=np.zeros((nNodes,3))
        
        for i in range(0,nNodes):
            Nodes[i,:]=numscan(f)[1:4]

        while stop==1:
            line=f.readline()
            if "$Elements" in line:
                stop=0
    elif n==3:
        f.seek(back)
        Nodes=np.zeros((nNodes,2))
        for i in range(0,nNodes):
            Nodes[i,:]=numscan(f)[1:3]

        while stop==1:
            line=f.readline()
            if "$Elements" in line:
                stop=0
    nElements=numscan(f)[0]
    

    while stop==0:
        back=f.tell()
        if numscan(f)[1]!=15:
            stop=1

    f.seek(back)
    if n==4:
        
        eLines=numscan(f)[4:8]
        while stop==1:
            
            eLines=np.vstack((eLines,numscan(f)[4:8]))
            back=f.tell()
            if numscan(f)[1]!=1or numscan(f)==[]:
                stop=0
            f.seek(back)
        
        
        eTriangles=numscan(f)[4:9]
        back=f.tell()
        if numscan(f) != []:
            while stop==0:
                eTriangles=np.vstack((eTriangles,numscan(f)[4:9]))
                back=f.tell()
                
                if "$End"in f.readline():
                    stop=1
                f.seek(back)
        f.seek(back)
    else:
        
        eLines=numscan(f)[4:7]
        while stop==1:
            
            eLines=np.vstack((eLines,numscan(f)[4:7]))
            back=f.tell()
            if numscan(f)[1]!=1 or numscan(f)==[]:
                stop=0
            f.seek(back)

        eTriangles=numscan(f)[4:8]

        while stop==0:
            
            eTriangles=np.vstack((eTriangles,numscan(f)[4:8]))
            back=f.tell()
            
            if "$End"in f.readline():
                stop=1
            f.seek(back)
    
##    f=open(geometry,'r')
##
##    while stop==1:
##        back=f.tell()
##        line=f.readline()
##        if "Point" in line:
##            stop=0
##            f.seek(back)
##    back=f.tell()
    ##    Points=numscan(f)[0:3]
    ##    while stop==0:
    ##        
    ##        Points=np.vstack((Points,numscan(f)[0:3]))
    ##        back=f.tell()
    ##        print Points
    ##        if "Line" in f.readline():
    ##            stop=0
    ##        f.seek(back)
 
    return Nodes, eLines#, eTriangles

def ReadSolverInput(Input):
    f=open(Input,'r')
    line=f.readline()
    while '$Solver' not in line:
        line=f.readline()
        if '$Solver' in line:
            
            Dimension=f.readline()
            BCType=f.readline()
            here=f.tell()
            if f.readline ==[]:
                params=[]
            else:
                f.seek(here)
                parameter=np.array(numscan(f))
                while True:
                    here=f.tell()
                    b=numscan(f)
                    if b==[]:
                        break
                    parameter=np.vstack((parameter,b))
            f.seek(here)
            Eq=f.readline()
            Type=f.readline()
            AnalisisParam=np.array(f.readline())
            SolverInput=[Dimension,BCType,parameter,Eq,Type,AnalisisParam]
            return SolverInput
            break
        if line=='':
            SolverInput=[]
            print 'Found no input for solver module' 
            return SolverInput
            break
            
        
        
        
    
    
    


























    
        
