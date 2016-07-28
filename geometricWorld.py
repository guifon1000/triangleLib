import pyevtk.hl as pv
import numpy as np
import triangles_to_VTK as tv
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import sys


name='zozozob'

Nu=1
Nv=1

#class Quad:

"""
      A -> B
     ^     |
    /      |
   /       |
  /        |
 /         v
D<---------C
           
"""
A=tl.Point(-0.05,0.05,10.)
B=tl.Point(0.,0.,10.)
C=tl.Point(0.,0.,0.)
D=tl.Point(-3.75,1.25,0.)

d={}
d['zero']=np.array([0.,0.])
tv.triangle_faces_to_VTK(name,
                      x=np.array([A.x,B.x,C.x,D.x]), y=np.array([A.y,B.y,C.y,D.y]), z=np.array([A.z,B.z,C.z,D.z]),
                      faces=np.array([[0,1,2],[2,0,3]]),
                      point_data=None,
                      cell_data= d)

pts=[]

ip=-1
for i in range(Nu):
    for j in range(Nv):
        ip+=1
        pts.append(tl.Point(,0.,0.))


print [p.x for p in pts]

tv.pointsToVTK(name,
                      x=np.array([p.x for p in pts]) , y=np.array([p.y for p in pts]), z=np.array([p.x for p in pts]),
                      data= None)

    
