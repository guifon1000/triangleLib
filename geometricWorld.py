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
import dynamics as dyn

name='sail'

Nu=5
Nv=5


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

A=tl.Point(-0.25,2.,0.15)
B=tl.Point(0.3,2.,-0.45)
C=tl.Point(0.5,-2.,-0.3)
D=tl.Point(-1.,-2.,0.35)
fs=vlm.freeStream(1.,0.,0.)

q0=tl.bigQuad(A,B,C,D)
q0.subQuad(Nu,Nv)
q0.facetize()


mat=dyn.materialPoint(0.,0.,0.,1.)

#A=vlm.fluidPoint(-42.,-50.,-100.)
#B=vlm.fluidPoint(42.,50.,100.)
#lv=vlm.lineVortex(A,B,1000.)


#plt.clf()
#leg=[]
#for x in np.arange(0.5,2.,0.05):
#    z=0.
#    absc=[]
#    ordo=[]
#    leg.append('x = '+str(x))
#    for y in np.arange(-5.,5.,0.01) :
#        plo = lv.getVelocity(x,y,z)
#        vp = tl.Vector(plo[0],plo[1],plo[2]).norm0
#        absc.append(y)
#        ordo.append(vp)
#    plt.plot(absc,ordo,'.')
#    plt.xlabel('y')
#    plt.ylabel('normV')
#    #plt.legend(leg)
#plt.show()
#


#plt.clf()
#Np=1000
#d=8.
#mat=np.zeros((Np,Np))
#for i in range(Np):
#    for j in range(Np):
#        z = 0.
#        x = -d+(2.*d)*i/float(Np)
#        y = -d+(2.*d)*j/float(Np)
#        v = lv.getVelocity(x,y,z)
#        val = tl.Vector(v[0],v[1],v[2]).norm0
#        mat[i][j]=val
#plt.pcolor(mat)
#plt.colorbar()
#plt.show()




tab=[]
rht=[]
for q in q0.quads:
    tab.append(vlm.dipolePanel(q.p0,q.p1,q.p2,q.p3))

sfce=vlm.Surface(tab,q0)
sfce.dipoleMatrix(sfce.tab)

plt.clf()
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
mu = np.dot(np.linalg.inv(sfce.M),rht)

sfce.addCellData(rht,'rht')
sfce.addCellData(mu,'mu')
sfce.writeVTK('mumu')

#print mu


for i,p in enumerate(tab):
    p.setIntensity(mu[i])
    elem.append(p)

#elem=[]
elem.append(fs)




A=tl.Point(-1.,-100.,0.)
B=tl.Point(-1.,100.,0.)


C=tl.Point(1.,-100.,0.)
D=tl.Point(1.,100.,-0.)



l0=vlm.lineVortex(A,B,1.)
l1=vlm.lineVortex(C,D,1.)



#elem.append(l0)
#elem.append(l1)










plt.clf()
Np=3
d=4.
mat=np.zeros((Np,Np))
for i in range(Np):
    for j in range(Np):
        y = 0.
        x = -d+(2.*d)*j/float(Np)
        z = -d+(2.*d)*i/float(Np)
        v=tl.Vector()
        for k,p in enumerate(elem):
            v += p.getVelocity(x,y,z)
        val = tl.Vector(v[0],v[1],v[2]).norm0
        mat[i][j]=val
plt.pcolor(mat)
plt.colorbar()
plt.show()


print '-----------------------------------------------------------'
vel=[]
normal=tl.Vector(0.,-1.,0.)

#elem.append(lv)
#elem.append(vlm.Source(1.,0.,-0.1,0.))
#elem.append(vlm.Doublet(1.,-1.,1.,0.,normal))
#elem.append(vlm.Source(-100.,0.,0.1,0.))

f=1.


mesh = vlm.Mesh(50,50,50,5.,5.,5.)
mesh.compute(elem)

mesh.writeVTK()
asciiArt.openFile('snoopy')
