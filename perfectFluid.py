import numpy as np
import bidLib as bid
import matplotlib.pyplot as plt


class PointElement(list):
    def __init__(self,x,y):
        super(PointElement,self).__init__([x,y])

class LinearElement(list):
    # the local y is always the outward normal
    def __init__(self,p0,p1):
        super(LinearElement,self).__init__([p0,p1])
        self.angleX = bid.angleX(self)
        self.normal = [-bid.normal(self)[i] for i in range(2)]
        pmat = np.zeros((2,2),dtype = float)
        pmat[0][0] = np.cos(self.angleX) 
        pmat[0][1] = np.sin(self.angleX) 
        pmat[1][0] = -np.sin(self.angleX) 
        pmat[1][1] = np.cos(self.angleX)
        self.pmat = pmat
        self.xa = p0[0]
        self.ya = p0[1]
        self.xb = p1[0]
        self.yb = p1[1]
    def locXY(self,x,y):
        ori = np.array([self[0][0],self[0][1]])
        prod = np.dot(self.pmat,np.array([x-ori[0],y-ori[1]]))
        out = prod
        return out
    def globXY(self,x,y):
        ori = np.array([self[0][0],self[0][1]])
        prod = np.dot(np.linalg.inv(self.pmat),np.array([x,y]))
        out = [prod[0]+ori[0],prod[1]+ori[1]]
        return out

def loadXfoilFile(name):
    f = open(name,'r').readlines()
    points = []
    for i in range(1,len(f)):
        l = f[i].split()
        p = PointElement(float(l[0]),float(l[1]))
        points.append(p)
    panels = []
    for i in range(len(points)-1):
        p0 = points[i]
        p1 = points[i+1]
        panels.append(LinearElement(p0,p1))
    plt.clf()
    for p in panels:
        plt.plot([p.xa,p.xb],[p.ya,p.yb])
        pt = p.globXY(0.01,0.01)
        plt.scatter(pt[0],pt[1])
    plt.axis('equal')
    plt.show() 
        


b = []
panels = []
A = PointElement(1.,1.)
B = PointElement(2.,-1)
pan = LinearElement(A,B)
panels.append(pan)
vinf = np.array([1.,0.])
for p in panels:
    b.append(np.dot(vinf,p.normal))

print b
print '--------'
#loadXfoilFile('profiles/NACAcamber0012.dat')
