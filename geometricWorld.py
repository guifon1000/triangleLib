import pyevtk.hl as pv
import numpy as np
import triangles_to_VTK as tv
import pyevtk.hl as vhl
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import sys
import vlm0 as vlm
import asciiArt

name='sail'

Nu=50
Nv=10


asciiArt.openFile('minusEtCortex2')
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

A=tl.Point(-1.,1.,0.01)
B=tl.Point(1.,1.,0.)
C=tl.Point(1.,-1.,0.)
D=tl.Point(-1.,-1.,0.01)


A=tl.Point(-1.,1.,0.08)
B=tl.Point(1.,0.4,0.)
C=tl.Point(1.,-1.,0.)
D=tl.Point(-1.,-1.,-0.08)


fs=vlm.freeStream(1.,0.,0.)

q0=tl.bigQuad(A,B,C,D)
q0.subQuad(Nu,Nv)
q0.facetize()




tab=[]
rht=[]
for q in q0.quads:
    tab.append(vlm.dipolePanel(q.p0,q.p1,q.p2,q.p3))

sfce=vlm.Surface(tab,q0)
sfce.dipoleMatrix(sfce.tab)

print sfce.M
plt.matshow(sfce.M)
plt.colorbar()
plt.show()

for q in sfce.tab:
    rht.append(-tl.dot(fs.v,q.pan.ng))
rht=np.array(rht)
print " right hand term"
print rht
print ''
print ''


elem=[]
mu = np.dot(np.linalg.inv(sfce.M),np.transpose(rht))

sfce.addCellData(rht,'rht')
sfce.addCellData(mu,'mu')
sfce.writeVTK('mumu')

#print mu


for i,p in enumerate(tab):
    p.setIntensity(mu[i])
    elem.append(p) 
print '-----------------------------------------------------------'
vel=[]
normal=tl.Vector(0.,-1.,0.)
elem.append(fs)
#elem.append(vlm.Source(1.,0.,-0.1,0.))
#elem.append(vlm.Doublet(1.,-1.,1.,0.,normal))
#elem.append(vlm.Source(-100.,0.,0.1,0.))

f=1.


#mesh = vlm.Mesh(10,10,10,10.,10.,10.)
#mesh.compute(elem)

#mesh.writeVTK()
asciiArt.openFile('snoopy')
