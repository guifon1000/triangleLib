import triangleLib as tl
import numpy as np


def angleX(p):
    A = tl.Point(p[0][0],p[0][1])
    B = tl.Point(p[1][0],p[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(1.,0.,0.)
    return tl.angle(v,u)


def normal(p):
    A = tl.Point(p[0][0],p[0][1])
    B = tl.Point(p[1][0],p[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(0.,0.,1.)
    pv = tl.cross(v,u) 
    return [pv[0],pv[1]]


