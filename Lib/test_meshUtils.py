 #! /usr/bin/python

from meshUtils import*

##### 1D Test ###
##coords0, elems0=mesh1D(-1,1,20)
####meshPlot(coords0,elems0,'none',True,True)
##
##facecolor='none'
##coords1, elems1=logmesh1D(-1,1,20,3,"R")
##print coords1, elems1
##meshPlot(coords1,elems1,facecolor,True,True)



### 2D Test ###

coords3, elems3 = meshtr2D(-1,1,-1,1,10,10)
print coords3
print elems3
facecolor = 'none'
meshPlot(coords3,elems3,facecolor,True,True)

##coords4, elems4 = meshsq2D(-1,1,-1,1,10,10)
##print coords4
##print elems4
##facecolor = 'none'
##meshPlot(coords4,elems4,facecolor,True,True)

##coords5, elems5 = logmeshtr2D(-1,1,-1,1,20,20,'L','R',2,2)
##print coords5
##print elems5
##facecolor = 'none'
##meshPlot(coords5, elems5,facecolor,True,True)


##coords5, elems6 = linlogmeshtr2D(-1,1,-1,1,5,5,'L',3)
##print coords6
##print elems6
##facecolor = 'none'
##meshPlot(coords6,elems6,facecolor,True,True)
