import numpy as np
import bidLib as bid



class PointElement(list):
    def __init__(self,x,y):
        super(PointElement,self).__init__([x,y])

class LinearElement(list):
    # the local y is always the outward normal
    def __init__(self,p0,p1):
        super(LinearElement,self).__init__([p0,p1])
        self.angleX = bid.angleX(self)
        self.normal = bid.normal(self)
        pmat = np.zeros((2,2),dtype = float)
        pmat[0][0] = np.cos(self.angleX) 
        pmat[0][1] = np.sin(self.angleX) 
        pmat[1][0] = -np.sin(self.angleX) 
        pmat[1][1] = np.cos(self.angleX) 
        self.pmat = pmat
    def locXY(self,x,y):
        ori = np.array([self[0][0],self[0][1]])
        return np.dot(self.pmat,np.array([x,y]))-ori

A = PointElement(1.,0.)
B = PointElement(0.,1.)
p = LinearElement(A,B)
print p.locXY(1.,1.)
print '--------'

