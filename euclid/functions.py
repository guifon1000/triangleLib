import sys
sys.path.append('./classes/')
from classes.Point import Point
from classes.Vector import Vector
from classes.Frame import Frame
from scipy import interpolate
import numpy as np

def cross(u,v):
    pv=[]
    pv.append(u[1]*v[2]-u[2]*v[1])
    pv.append(u[2]*v[0]-u[0]*v[2])
    pv.append(u[0]*v[1]-u[1]*v[0])
    vp=Vector(pv)
    return vp


def dot(u,v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]



def matrix_to_quaternion(m):
    diag=np.diag(m)
    np.append(diag,1.)
    
    tr= np.trace(m)+1.
    if tr>0.:
        s=0.5/np.sqrt(tr)
        x=(m[2,1]-m[1,2])*s
        y=(m[0,2]-m[2,0])*s
        z=(m[1,0]-m[0,1])*s
        w=0.25/s
        return Quaternion(w,x,y,z)
         


def parameter_frame(tck, s, mode = 'frenet'):
    from classes.Frame import Frame
    basis = []
    t =  interpolate.splev(s , tck, der = 0)
    orig = Point([t[0], t[1], t[2]])
    if mode == 'frenet':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        t =  interpolate.splev(s , tck, der = 2)
        xpp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T.unit())
        B = cross(xp, xpp)
        basis.append(B.unit())
        N = cross(B, T)
        basis.append(N.unit())

    if mode == 'Xnat':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T)
        B = cross((1.,0.,0.), T)
        basis.append(B)
        N = cross(T, B)
        basis.append(N)
    matrix = np.zeros((3,3), dtype =float)
    for i in range(3) : matrix[i] = basis[i]
    return Frame((orig, matrix))

