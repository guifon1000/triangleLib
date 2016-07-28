import pyevtk.hl as pv
import numpy as np
import triangles_to_VTK as tv
import pyevtk.hl as vhl
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import sys


name='zozozob'

Nu=80
Nv=150

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
tv.triangle_faces_to_VTK(name+'all',
                      x=np.array([A.x,B.x,C.x,D.x]), y=np.array([A.y,B.y,C.y,D.y]), z=np.array([A.z,B.z,C.z,D.z]),
                      faces=np.array([[0,1,2],[2,0,3]]),
                      point_data=None,
                      cell_data= d)

pts=[]

for i in range(Nu):
    realI=float(i)/float(Nu-1)
    cprime=tl.Point(C.x+realI*(B.x-C.x),C.y+realI*(B.y-C.y),C.z+realI*(B.z-C.z))
    dprime=tl.Point(D.x+realI*(A.x-D.x),D.y+realI*(A.y-D.y),D.z+realI*(A.z-D.z))
    for j in range(Nv):
        realJ=float(j)/float(Nv-1)
        pts.append(tl.Point(cprime.x+realJ*(dprime.x-cprime.x),cprime.y+realJ*(dprime.y-cprime.y),cprime.z+realJ*(dprime.z-cprime.z)))
     
ip=0
faces=[]
for i in range(Nu-1):
    for j in range(Nv-1):
        ip+=1
        if np.mod(ip,Nv)==0 : ip+=1
        print ip,ip+1,ip+Nu+1,ip+Nu
        print ip,ip+1,ip+Nu+1
        print ip,ip+Nu+1,ip+Nu
        print '----------------------'
        faces.append([ip-1,ip,ip+Nv])
        faces.append([ip-1,ip+Nv,ip+Nv-1])
print faces
tv.triangle_faces_to_VTK(name+'rect',
                      x=np.array([p.x for p in pts]) , y=np.array([p.y for p in pts]), z=np.array([p.z for p in pts]),
                      faces=np.array(faces),
                      point_data=None,
                      cell_data= None)


