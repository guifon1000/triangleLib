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
import MPLout as mpl
#asciiArt.openFile('minusEtCortex2')

testMSH = tl.quadsPanel('testGMSH')
plt.clf()
mpl.plot3(testMSH)


plt.show()

fs = vlm.freeStream(1.0,0.,0.)

tab=[]
rht=[]
for q in testMSH.quads:
    tab.append(vlm.dipolePanel(q.p0,q.p1,q.p2,q.p3))
sfce=vlm.Surface(tab,testMSH)
sfce.dipoleMatrix(sfce.tab)

plt.clf()
plt.matshow(sfce.M)
plt.colorbar()
plt.show()

for q in sfce.tab:
    rht.append(-tl.dot(fs.v,q.pan.ng))
rht=np.array(rht)

elem=[]
mu = np.dot(np.linalg.inv(sfce.M),rht)
sfce.addCellData(mu,'mu')
sfce.writeVTK('mumu')


for i,p in enumerate(tab):
    p.setIntensity(mu[i])
    elem.append(p)

elem.append(fs)
mesh = vlm.Mesh(50,30,20,5.,3.,2.)
mesh.compute(elem)

mesh.writeVTK()
asciiArt.openFile('snoopy')
