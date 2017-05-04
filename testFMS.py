import mesh 
import triangleLib as tl
import matplotlib.pyplot as plt
import numpy as np
import bidLib as l2d
import sys
sys.path.append('./profiles/')
from splineProfileMultiParam import Profile
l = []
N = 100
facdom = 2.
r = 0.2
for i in range(N):
    teta = float(i)*2.*np.pi/float(N)
    p = (r*np.cos(teta) , r*np.sin(teta))


    #l.append(p)


f = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.04,0.029],npt = 50)

revX = f.x[::-1]
revextra = f.extra[::-1]
for i in range(len(revX)):
    p = (r*revX[i], r*revextra[i])
    l.append(p)
for i in range(1,len(f.x)):
    p = (r*f.x[i], r*f.intra[i])
    l.append(p)
pl = mesh.Polyline2D(l)
pl.box_2d()
