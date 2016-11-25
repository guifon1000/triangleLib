import triangleLib as tl
import numpy as np


def angleX(p):
    A = tl.Point(p[0][0],p[0][1])
    B = tl.Point(p[1][0],p[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(1.,0.,0.)
    return tl.angle(u,v)


def normal(p):
    A = tl.Point(p[0][0],p[0][1])
    B = tl.Point(p[1][0],p[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(0.,0.,1.)
    return [tl.cross(u,v)[0],tl.cross(u,v)[1]]


