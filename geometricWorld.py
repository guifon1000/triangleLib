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
import MPLout as mpl
import sys
import pygmsh as pg
import meshio as meshio
#asciiArt.openFile('minusEtCortex2')

geom = pg.Geometry()

# Draw a cross.

#    [-0.5,  0.0, 0.0],
#    [-0.1, -0.1, 0.0],
#    [0.0,  -0.5, 0.0],
#    [0.1,  -0.1, 0.0],
#    [0.5,   0.0, 0.0],
#    [0.1,   0.1, 0.0]
#    ],
#    lcar=0.05
#    )

p0=geom.add_point([-0.5,  -0.25, 0.0],lcar=0.1)
p1=geom.add_point([0.0,  0.0 , 0.0],lcar=0.1)
p2=geom.add_point([0.5,  0.0 , 0.0],lcar=0.1)

l0 = geom.add_line(p0,p1)
l1 = geom.add_line(p1,p2)
#geom._GMSH_CODE.append('Recombine Surface {%s};' % poly)

axis = [0, 0, 1]

a = geom.extrude(
    'Line{%s}' % l0,
    translation_axis=axis,
    rotation_axis=axis,
    point_on_axis=[0, 0, 0],
    angle=0.8 / 6.0 * np.pi
    )

geom._GMSH_CODE.append('Recombine Surface {%s[]};' % a)
a = geom.extrude(
    'Line{%s}' % l1,
    translation_axis=axis,
    rotation_axis=axis,
    point_on_axis=[0, 0, 0],
    angle=0.8 / 6.0 * np.pi
    )

geom._GMSH_CODE.append('Recombine Surface {%s[]};' % a)

fg = open('test0.geo','w')
for l in geom.get_code():
    fg.write(l)
points, cells = pg.generate_mesh(geom)
meshio.write('test.vtu', points, cells)
plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#for p in points:
#    ax.scatter(p[0], p[1], p[2]) 
for t in cells['quad']:
    p0=points[t[0]]
    p1=points[t[1]]
    p2=points[t[2]]
    p3=points[t[3]]
    ax.plot([p0[0],p1[0],p2[0],p3[0],p0[0]],\
            [p0[1],p1[1],p2[1],p3[1],p0[1]],\
            [p0[2],p1[2],p2[2],p3[2],p0[2]],'k')
plt.show()



testMSH = tl.quadsPanel(points,cells['quad'])
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
