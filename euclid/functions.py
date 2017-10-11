import sys
sys.path.append('./classes/')
from classes.Vector import Vector




def cross(u,v):
    pv=[]
    pv.append(u[1]*v[2]-u[2]*v[1])
    pv.append(u[2]*v[0]-u[0]*v[2])
    pv.append(u[0]*v[1]-u[1]*v[0])
    vp=Vector(pv)
    return vp


def dot(u,v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]
