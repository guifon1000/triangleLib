import numpy as np
import sys
sys.path.append('../triangleLib/')
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
sys.path.append('../VortexLatticeMethod/')
import vlm
import sys
import pygmsh as pg
import meshio as meshio
import mesh




#asciiArt.openFile('./asciiArt/minusEtCortex')

print "zobz"
P = [tl.Point(-1.0,10.0,0.1),\
        tl.Point(1.0,10.0,-0.3),\
        tl.Point(1.0,-10.0,-0.7),\
        tl.Point(-1.0,-10.0,0.4)]

Nu = 5
Nv = 5



ob = mesh.surfaceMesh(P,Nu,Nv)

ob.MPLout()
fs = vlm.freeStream(1.0,0.,0.)

tab=[]
rht=[]
for q in ob.quads:
    p0 = ob.points[q[0]]
    p1 = ob.points[q[1]]
    p2 = ob.points[q[2]]
    p3 = ob.points[q[3]]
    t = vlm.dipolePanel(p0,p1,p2,p3)
    tab.append(t)
sfce=vlm.Surface(tab,ob)
sfce.dipoleMatrix(sfce.tab)

plt.matshow(sfce.M)
plt.colorbar()
plt.show()

for q in sfce.tab:
    rht.append(-tl.dot(fs.v,q.pan.ng))
rht=np.array(rht)
elem=[]
mu = np.dot(np.linalg.inv(sfce.M),rht)
sfce.addCellData(mu,'mu')
meshio.write('test.vtu', ob.points, ob.cells, cell_data={'mu':mu})# ells['quad'])    #{'mu':mu}

M = mesh.volumicMesh(4, 4, 4, 50 ,10 , 50)
for i,p in enumerate(tab):
    p.setIntensity(mu[i])
    elem.append(p)
elem.append(fs)
M.compute(elem)
#asciiArt.openFile('./asciiArt/snoopy')
