import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import random,randint
from time import sleep
import math
import mesh
import os
#import meshio








class Point(list):
    def __init__(self,*args,**kwargs):
        """
            class Point -> array of three floats (the coordinates), only cartesian
            also has the attributes .x , .y and .z --> same as self[0], self[1], self[2]
            can have an index
        """
        super(Point,self).__init__(args)
        #self.append([float(x),float(y),float(z)])
        self.x=self[0]
        self.y=self[1]
        self.z=self[2]
        if kwargs.has_key('name'):
            self.name = kwargs['name']
        else:
            self.name = None
        

    def setName(self,st):
        self.name = st

    def setIndex(self,i):
        """ 
            set the index of the point
        """
        self.index = i

    def ran(self):
        """
            Shitty random point. Should it really call sleep ?
        """
        self.x=random()
        sleep(0.1)
        self.y=random()
        sleep(0.1)
        self.z=random()
        sleep(0.1)

    def addi(self,e):
        """
            returns self + e (e is a Point too)
        """
        pp=Point()
        pp.setPos(self.x+e.x,self.y+e.y,self.z+e.z)
        return pp






def distance(p1,p2):
    """
        the distance between 2 Points
    """
    return math.sqrt((p2.x-p1.x)**2+(p2.y-p1.y)**2+(p2.z-p1.z)**2)




class Segment(list):
    def __init__(self,p1,p2):
        """
            class Segment -> array of 2 points
        """
        super(Segment,self).__init__([p1,p2])
    def middle(self):
        """
            finds the middle of the Segment
        """
        self.mid=Point()
        self.mid.x=0.5*(self.p1.x+self.p2.x)
        self.mid.y=0.5*(self.p1.y+self.p2.y)
        self.mid.z=0.5*(self.p1.z+self.p2.z)


class Vector(list):
    def __init__(self, *args):
        """
            class Vector   
        """
        if len(args)==3 :
            for i in range(3) : self.append(args[i])
        if len(args)==2 :
            for i in range(3) : 
                self.append(args[1][i]-args[0][i])

        self.x=self[0]
        self.y=self[1]
        self.z=self[2]
        self.norm0 = math.sqrt(self.x**2+self.y**2+self.z**2)

    def __add__(self,other):
        """
            Adds an other Vector
        """
        return Vector(self.x+other.x,self.y+other.y,self.z+other.z)

    def __mul__(self,alpha):
        """
            multiplies the vector by a real
        """
        return Vector(self.x*alpha,self.y*alpha,self.z*alpha)

    def __rmul__(self,alpha):
        """
            right multiplication by a real
        """
        return Vector(self.x*alpha,self.y*alpha,self.z*alpha)

    def calcNorm0(self):
        """
            computes the norm of the vector
        """
        self.norm0=math.sqrt(self.x**2+self.y**2+self.z**2)

    def set(self,x,y,z):
        """
            sets the coordinates of the Vector
        """
        self.x=x
        self.y=y
        self.z=z
        self.calcNorm0()

    def normalize(self):
        """
            normalize the vector
        """
        self.calcNorm0()
        self.x=self.x/self.norm0
        self.y=self.y/self.norm0
        self.z=self.z/self.norm0

    def draw(self,ax,p):
        """
            adds the point to a 3d existing matplotlib view, ax
        """
        pp=p.addi(self)
        ax.plot([p.x,pp.x],[p.y,pp.y],[p.z,pp.z])


    def projectOnCoordinatesSystem(self,u0,v0,w0):
        """
            returns the Vector given in the 3 Vectors u0,v0,w0 basis
        """
        cx=u0.x*self.x+u0.y*self.y+u0.z*self.z
        cy=v0.x*self.x+v0.y*self.y+v0.z*self.z
        cz=w0.x*self.x+w0.y*self.y+w0.z*self.z
        v=Vector()
        v.set(cx,cy,cz)
        return v

class Polyline(list):
    def __init__(self,*args):
        for p in args:
            self.append([float(p[i]) for i in range(len(p)) ])
                    


class Panel(object):
    """
    
    p0 ----> p1                  p0 <------ p1                  p1                                                                                                     
     ^ \   / |                      \   0  ^                   ^ |                                                                                                  
     |  \ /  |    ____|\             \    /                   /  |                                                                                              
     |   X c0|          \             \  /          +        /   |
     |  / \  |    ____  /              v                   c0    |                                                                                        
     | /   \ v        |/               c0                    ^   | 
    p3 <--- p2                                                \  | 
                                                               \ v
                                                                p2
                                   t:                             

    """                                                         
    def __init__(self,p0,p1,p2,p3):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        x0=0.25*(p0[0]+p1[0]+p2[0]+p3[0]) 
        y0=0.25*(p0[1]+p1[1]+p2[1]+p3[1]) 
        z0=0.25*(p0[2]+p1[2]+p2[2]+p3[2]) 
        cent0=Point(x0,y0,z0)
        self.cg = cent0
        
        d0=Vector(cent0,p0)
        d1=Vector(cent0,p1)
        d2=Vector(cent0,p2)
        d3=Vector(cent0,p3)
        

        n0 = cross(d3,d0,norm=False)
        n1 = cross(d0,d1,norm=False)
        n2 = cross(d1,d2,norm=False)
        n3 = cross(d2,d3,norm=False)
        
        self.ng=Vector(0.25*(n0.x+n1.x+n2.x+n3.x),\
                    0.25*(n0.y+n1.y+n2.y+n3.y),\
                    0.25*(n0.z+n1.z+n2.z+n3.z))
        self.ng.normalize()

def cross(u,v, norm=True):
    pv=[]
    pv.append(u.y*v.z-u.z*v.y)
    pv.append(u.z*v.x-u.x*v.z)
    pv.append(u.x*v.y-u.y*v.x)
    vp=Vector(pv[0],pv[1],pv[2])
    if norm:vp.normalize()
    return vp

def dot(u,v):
    return u.x*v.x+u.y*v.y+u.z*v.z

def angle(u,v):
    
    angle = np.arctan2(cross(u,v,norm=False).norm0,dot(u,v))
    if angle == 0.:return 0.
    if dot(cross(u,v),Vector(0.,0.,1.0))>=0.:
        return angle
    else :
        return -angle



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
    def set3Points(self,*args):
        v0=Vector(args[0],args[1])
        v1=Vector(Args[0],args[2])
        v2=cross(v0,v1)
        v=v2
        self.a=v[0]
        self.b=v[1]
        self.c=v[2]
        self.d=-(v[0]*args[0][0]+v[1]*args[0][1]+v[2]*args[0][2])
    def setPointNormal(self,p,normal):
        print 'zob'
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


class Triangulation(object):
    def __init__(self,**kwargs):
        if (kwargs.has_key('file')) and (kwargs['file'].endswith('.msh') or \
					kwargs['file'].endswith('.MSH')):
            f = open('./'+str(kwargs['file']),'r').readlines()
            _vertices = []
            app = None
            nv = None
            for ist in range(len(f)):
                l = f[ist]
                if l.startswith('$Nodes'):
                    nv = int(f[ist+1])
                    app = ist
                    break
            else:
 		print 1/0
            app +=1
            print 'there are '+str(nv)+' vertices strating at line '+str(app+1)+' :'+f[app+1]
            print 'ending at line '+str(app+nv)+' : '+f[app+nv]
            for i in range(app+1,app+nv+1):
                l=f[i].split()
                _vertices.append((float(l[1]), float(l[2]), float(l[3])))
            self.vertices = _vertices
            _faces = []
            for i in range(ist+nv+1,len(f)):
                lp = f[i].split()
                if (len(lp) > 3) and (lp[1] == '2') :
                        _faces.append((int(lp[5]),\
                                          int(lp[6]),\
                                          int(lp[7])))
            self.faces = _faces
            _normals = []
            for fa in self.faces :
                p0 = Point(*self.vertices[fa[0]-1])
                p1 = Point(*self.vertices[fa[1]-1])
                p2 = Point(*self.vertices[fa[2]-1])
                v1 = Vector(float(p1[0])-float(p0[0]),float(p1[1])-float(p0[1]),float(p1[2])-float(p0[2]))
                v2 = Vector(float(p2[0])-float(p1[0]),float(p2[1])-float(p1[1]),float(p2[2])-float(p1[2]))
                n = cross(v1,v2,norm=False)
                _normals.append((n[0],n[1],n[2]))
            self.normals = _normals

    def set_position(self):
        for ve in self.vertices:
            ve[0] += self.posg[0]
            ve[1] += self.posg[1]
            ve[2] += self.posg[2]
 
class Triangle(list):
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
    def __init__(self,*args):
        super(Triangle,self).__init__(args)
        self.s1=Segment(self[0], self[1])
        self.s2=Segment(self[1], self[2])
        self.s3=Segment(self[2], self[3])
        p1 = self[0]
        p2 = self[1]
        p3 = self[2]
        self.u=Vector(p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2])
        self.v=Vector(p3[0]-p2[0],p3[1]-p2[1],p3[2]-p2[2])
        self.w=Vector(p1[0]-p3[0],p1[1]-p3[1],p1[2]-p3[2])
        
        #self.s1.setPoints(self.p1,self.p2)
        #self.s2.setPoints(self.p2,self.p3)
        #self.s3.setPoints(self.p3,self.p1)

    def setPoints(self,p1,p2,p3):
        self.p1=p1
        self.p2=p2
        self.p3=p3
    

    def computeNormal(self):    # wrong
        return cross(self.u,self.v,norm = False) 
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
        self.circumcenter=vor   # vor ? really ?
        return vor
    
    def gravityCenter(self):
	pcg=Point()
	xg=(self.p1.x+self.p2.x+self.p3.x)/3.
	yg=(self.p1.y+self.p2.y+self.p3.y)/3.
	zg=(self.p1.z+self.p2.z+self.p3.z)/3.
	pcg.setPos(xg,yg,zg)
	self.cg=pcg
	return pcg


    


    def draw(self,ax):
        x=[self.p1.x,self.p2.x,self.p3.x]
        y=[self.p1.y,self.p2.y,self.p3.y]
        z=[self.p1.z,self.p2.z,self.p3.z]
        x.append(self.p1.x)
        y.append(self.p1.y)
        z.append(self.p1.z)
        #self.middleIze()
        #ax.scatter(self.s1.mid.x,self.s1.mid.y,self.s1.mid.z)
        #ax.scatter(self.s2.mid.x,self.s2.mid.y,self.s2.mid.z)
        #ax.scatter(self.s3.mid.x,self.s3.mid.y,self.s3.mid.z)
        ax.plot(x,y,z)
        


if __name__=='__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # 2D polygone  
    N = 4
    index=-1
    points = []
    points.append(Point(-1.,1.,0.))
    points.append(Point(1.,1.,0.))
    points.append(Point(1.,-1.,0.))
    points.append(Point(-1.,-1.,0.))
    zob = Vector(points[0],points[1])
    triangulation = Triangulation(file='sphere.msh')
    links = []
    links.append([0,1])
    links.append([1,2])
    links.append([2,3])
    links.append([3,0])
    stl = mesh.Triangulation(stl = os.getcwd()+'/meshing/Bobomb_HD.stl')
    tt = []
    pts = stl.points
    cd = []
    for tri in stl.tris:
        p0 = Point(pts[tri[0]][0],pts[tri[0]][1],pts[tri[0]][2])
        p1 = Point(pts[tri[1]][0],pts[tri[1]][1],pts[tri[1]][2])
        p2 = Point(pts[tri[2]][0],pts[tri[2]][1],pts[tri[2]][2])
        tri0 = Triangle(p0,p1,p2)
        norm = tri0.computeNormal()
        cd.append(norm)
        tt.append([norm[0],norm[1],norm[2]])
        #tri0.draw(ax)
    #meshio.write('test.vtu', stl.points,stl.cells, cell_data = {'normal' : np.array(tt)})# ells['quad'])    #{'mu':mu}
    meshio.write('test.vtu', stl.points,stl.cells, cell_data = None)# ells['quad'])    #{'mu':mu}
    ax.axis('equal')
    plt.show()
    
