import triangleLib as tl
import numpy as np



def angle(u,v):
    return tl.angle(tl.Vector(u[0],u[1]),tl.Vector(v[0],v[1]))


def cross(u,v):
    return tl.cross(tl.Vector(u[0],u[1]),tl.Vector(v[0],v[1]))[0:-1]


def angleX(panel):
    A = tl.Point(panel[0][0],panel[0][1])
    B = tl.Point(panel[1][0],panel[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(1.,0.,0.)
    return angle(v,u)


def normal(p):
    A = tl.Point(p[0][0],p[0][1])
    B = tl.Point(p[1][0],p[1][1])
    u = tl.Vector()
    u.fromPoints(A,B)
    v = tl.Vector(0.,0.,1.)
    pv = tl.cross(v,u) 
    return [pv[0],pv[1]]


