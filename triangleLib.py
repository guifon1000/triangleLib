import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import random
from time import sleep
import math



class Point:
    def __init__(self,x=0.,y=0.,z=0.):
        self.x=x
        self.y=y
        self.z=z
        self.index=-1
    def setPos(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
    def ran(self):
        self.x=random()
        sleep(0.1)
        self.y=random()
        sleep(0.1)
        self.z=random()
        sleep(0.1)
    def addi(self,e):
        pp=Point()
        pp.setPos(self.x+e.x,self.y+e.y,self.z+e.z)
        return pp

def distance(p1,p2):
    return math.sqrt((p2.x-p1.x)**2+(p2.y-p1.y)**2+(p2.z-p1.z)**2)


class Segment:
    def __init__(self):
        self.p1=Point()
        self.p2=Point()
    def setPoints(self,p1,p2):
        self.p1=p1
        self.p2=p2
    def middle(self):
        self.mid=Point()
        self.mid.x=0.5*(self.p1.x+self.p2.x)
        self.mid.y=0.5*(self.p1.y+self.p2.y)
        self.mid.z=0.5*(self.p1.z+self.p2.z)



class Vector:
    def __init__(self):
        self.x=0.
        self.y=0.
        self.z=0.
    def fromPoints(self,A,B):
        self.x=B.x-A.x
        self.y=B.y-A.y
        self.z=B.z-A.z
        self.calcNorm0()
    def calcNorm0(self):
        self.norm0=math.sqrt(self.x**2+self.y**2+self.z**2)
    def set(self,x,y,z):
        self.x=x
        self.y=y
        self.z=z
        self.calcNorm0()
    def normalize(self):
        self.calcNorm0()
        self.x=self.x/self.norm0
        self.y=self.y/self.norm0
        self.z=self.z/self.norm0
    def draw(self,ax,p):
        pp=p.addi(self)
        ax.plot([p.x,pp.x],[p.y,pp.y],[p.z,pp.z])
        
    def projectOnCoordinatesSystem(self,u0,v0,w0):
        cx=u0.x*self.x+u0.y*self.y+u0.z*self.z
        cy=v0.x*self.x+v0.y*self.y+v0.z*self.z
        cz=w0.x*self.x+w0.y*self.y+w0.z*self.z
        v=Vector()
        v.set(cx,cy,cz)
        return v


        
def cross(u,v, norm=True):
    pv=[]
    pv.append(u.y*v.z-u.z*v.y)
    pv.append(u.z*v.x-u.x*v.z)
    pv.append(u.x*v.y-u.y*v.x)
    vp=Vector()
    vp.set(pv[0],pv[1],pv[2])
    if norm:vp.normalize()
    return vp

class Plane:
    """
    defines a plane 
    """
    def __init__(self):
        self.a=0.
        self.b=0.
        self.c=0.
        self.d=0.
    def set(self,x):
        self.a=x[0]
        self.b=x[1]
        self.c=x[2]
        self.d=x[3]
    def setMediator(self,A,B):
        s=Segment()
        s.setPoints(A,B)
        m=s.middle()
        v=Vector()
        v.fromPoints(s.mid,B)
        v.normalize()
        self.a=v.x
        self.b=v.y
        self.c=v.z
        self.d=-(v.x*s.mid.x+v.y*s.mid.y+v.z*s.mid.z)
    def set3Points(self,A,B,C):
        v0=Vector()
        v1=Vector()
        v0.fromPoints(A,B)
        v1.fromPoints(A,C)
        v2=cross(v0,v1)
        v=v2
        self.a=v.x
        self.b=v.y
        self.c=v.z
        self.d=-(v.x*A.x+v.y*A.y+v.z*A.z)
    def draw(self,ax):
        xx, yy = np.meshgrid(range(10), range(10))
        z=(-self.d-self.a*xx-self.b*yy)/self.c
        ax.plot_surface(xx, yy, z)


def intersect3Planes(p1,p2,p3):
    A=np.array([[p1.a,p1.b,p1.c],[p2.a,p2.b,p2.c],[p3.a,p3.b,p3.c]])
    b=np.array([-p1.d,-p2.d,-p3.d]).transpose()
    sol=np.linalg.solve(A, b)
    pvor=Point()
    pvor.setPos(sol[0],sol[1],sol[2])
    return pvor



 
class Triangle:
    """    
              
               3
              
             / ^
            /  |
       s3  /   |
          / 0  |s2
         /     |
        v          
       1  ---> 2    (  0 = out normal )
          s1

    """
    def __init__(self):
        self.p1=Point()
        self.p2=Point()
        self.p3=Point()
        self.s1=Segment()
        self.s2=Segment()
        self.s3=Segment()
        self.u=Vector()
        self.v=Vector()
        self.w=Vector()
        self.s1.setPoints(self.p1,self.p2)
        self.s2.setPoints(self.p2,self.p3)
        self.s3.setPoints(self.p3,self.p1)

    def setPoints(self,p1,p2,p3):
        self.p1=p1
        self.p2=p2
        self.p3=p3
    def createSys(self):
        """
        creates a coordinate system (u,v,w) for which u is aligned with the first
        edge of the triangle
        """
        v0=Vector()
        v1=Vector()
        v0.set(self.p2.x-self.p1.x,self.p2.y-self.p1.y,self.p2.z-self.p1.z)
        v1.set(self.p3.x-self.p2.x,self.p3.y-self.p2.y,self.p3.z-self.p2.z)
        v0.normalize()
        v1.normalize()
        v2=cross(v0,v1)
        v3=cross(v2,v0)
        v2.normalize()
        v3.normalize()
        self.u=v0
        self.v=v3
        self.w=v2

    def createLocalSys_X(self):
        """
        creates a coordinate system (u,v,w) for which w is the normal of the triangle
        (as calculated by createSys) and we look for u to be as aligned as possible with
        the global X direction
        """
        #if self.w.norm0 : print self.w.norm0
        x=Vector()
        x.set(1.,0.,0.)
        try :
            y2=cross(self.w,x)
            y2.normalize()
        except :
            y2=Vector()
            y2.set(0.,1.,0.)
        x2=Vector()
        x2=cross(y2,self.w)
        x2.normalize()
        
        self.localSys=[x2,y2,self.w]
        #print x2.norm0, y2.norm0, self.w.norm0
        #print y2.norm0


    def middleIze(self):
        t.s1.middle()
        t.s2.middle()
        t.s3.middle()

    def circumCenter(self):
        self.createSys()
        A1=Plane()
        A2=Plane()
        A3=Plane()
        A1.setMediator(self.p1,self.p2)
        A2.setMediator(self.p2,self.p3)
        A3.set3Points(self.p1,self.p2,self.p3)
        vor=intersect3Planes(A1,A2,A3)
        self.circumcenter=vor
        return vor
    
    def gravityCenter(self):
	pcg=Point()
	xg=(self.p1.x+self.p2.x+self.p3.x)/3.
	yg=(self.p1.y+self.p2.y+self.p3.y)/3.
	zg=(self.p1.z+self.p2.z+self.p3.z)/3.
	pcg.setPos(xg,yg,zg)
	self.cg=pcg
	return pcg


    


    def draw(t,ax):
        x=[t.p1.x,t.p2.x,t.p3.x]
        y=[t.p1.y,t.p2.y,t.p3.y]
        z=[t.p1.z,t.p2.z,t.p3.z]
        x.append(t.p1.x)
        y.append(t.p1.y)
        z.append(t.p1.z)
        t.middleIze()
        ax.scatter(t.s1.mid.x,t.s1.mid.y,t.s1.mid.z)
        ax.scatter(t.s2.mid.x,t.s2.mid.y,t.s2.mid.z)
        ax.scatter(t.s3.mid.x,t.s3.mid.y,t.s3.mid.z)
        ax.plot(x,y,z)
        


if __name__=='__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    



    t=Triangle()
    t.p1.ran()
    t.p2.ran()
    t.p3.ran()
    vor=t.circumCenter()
    ax.scatter(vor.x,vor.y,vor.z)
    ax.axis('equal')
    plt.show()
    
