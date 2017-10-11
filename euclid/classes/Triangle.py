import sys
import numpy as np
from Vector import Vector
from Point import Point
sys.path.append('../')
from functions import cross

lcar = 0.1

class Triangle(list):
    """    
               
                self[2]
                   ^
                 /  \
                /    \
           s3  /      \ s0
              /        \
             /     0    \
            v            \
     self[0] -----------> self[1]    (  0 = out normal )
                   s1
    """
    def __init__(self, *largs):
        super(Triangle, self).__init__(*largs)

    
    def __neg__(self):
        return Triangle((self[0], self[2], self[1]))



    @property
    def cg(self):
        return Point(np.mean(self   ,axis=0))

    @property
    def normal(self):
        u = Vector((self[1][0]-self[0][0],\
    			self[1][1]-self[0][1],\
    			self[1][2]-self[0][2]))
        v = Vector((self[2][0]-self[1][0],\
    			self[2][1]-self[1][1],\
     			self[2][2]-self[1][2]))
        return cross(u,v)

    def pop_to_geom(self, geom):
        datas = []
        datas.append(geom.add_point(self[0], lcar))
        datas.append(geom.add_point(self[1], lcar))
        datas.append(geom.add_point(self[2], lcar))
        datas.append(geom.add_line(datas[0], datas[1]))
        datas.append(geom.add_line(datas[1], datas[2]))
        datas.append(geom.add_line(datas[2], datas[0]))
        datas.append(geom.add_line_loop((datas[3], datas[4], datas[5])))
        datas.append(geom.add_plane_surface(datas[6]))
