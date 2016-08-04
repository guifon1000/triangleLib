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

Nu=10
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

A=tl.Point(1.,1.,0.)
B=tl.Point(-1.,1,0.)
C=tl.Point(-1.,-1.,0.)
D=tl.Point(1.,-1.,0.)


#q0=tl.Quad(A,B,C,D)
#q0.subQuad(Nu,Nv)
#q0.facetize()

pan=tl.Panel(A,B,C,D)

vel=[]
elem=[]
fs=vlm.freeStream(10.,0.,0.)
normal=tl.Vector(0.,-1.,0.)
elem.append(fs)
#elem.append(vlm.Source(1.,0.,-0.1,0.))
#elem.append(vlm.Doublet(1.,-1.,1.,0.,normal))
#elem.append(vlm.Source(-100.,0.,0.1,0.))

f=1.


elem.append(vlm.lineVortex(A,B,f))
elem.append(vlm.lineVortex(B,C,f))
elem.append(vlm.lineVortex(C,D,f))
elem.append(vlm.lineVortex(D,A,f))

mesh = vlm.Mesh(50,50,50,100.,100.,100.)
mesh.compute(elem)
mesh.writeVTK()
asciiArt.openFile('snoopy')
