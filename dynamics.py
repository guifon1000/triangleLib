import triangleLib as tl
import matplotlib.pyplot as plt
import sympy as sp





class materialPoint(tl.Point):
    def __init__(self,x,y,z,m):
        super(materialPoint,self).__init__(x,y,z)
        self.mass = m

class Spring:
    def __init__(self,p1,p2,k):
        self.l0=tl.distance(p1,p2)
        self.k = k



if __name__ == '__main__':
    A = materialPoint(-1.,0.,0.,1.)
    B = materialPoint(1.,0.,0.,1.)
    s = Spring(A,B,1.)
    print 'ok'
