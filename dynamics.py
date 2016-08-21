import triangleLib as tl
import matplotlib.pyplot as plt
import sympy as sp
import numpy as np




class materialPoint(tl.Point):
    def __init__(self,x,y,z,m,dof='0'):
        super(materialPoint,self).__init__(x,y,z)
        self.mass = m
        self.dof = dof
class Spring:
    def __init__(self,p1,p2,k):
        self.l0=tl.distance(p1,p2)
        self.k = k



if __name__ == '__main__':
    A = materialPoint(-1.,0.,0.,1.,'0')
    B0 = materialPoint(1.,0.,0.,1.,'xyz')
    s = Spring(A,B0,10.)
    B = tl.Point(10.,10.,10.)
    Fab = tl.Vector()
    Fab.fromPoints(B,A)
    Fab*=s.k*(tl.distance(A,B)-s.l0)
    matK=np.eye(3)
    print s.k*matK
    print 'ok'


    # creation of a set of Points

