import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from random import random,randint
from time import sleep
import math
import mesh
import os
import json
#import meshio








class Point(list):
    def __init__(self,*largs, **kwargs):
        """
            class Point -> array of three floats (the coordinates), only cartesian
            also has the attributes .x , .y and .z --> same as self[0], self[1], self[2]
            can have an index
        """
        super(Point,self).__init__(*largs)
        #self.append([float(x),float(y),float(z)])
        self.x=self[0]
        self.y=self[1]
        self.z=self[2]

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
            returns self + e (e is a Vector)
        """
        pp=Point([self[0]+e[0],self[1]+e[1],self[2]+e[2]])
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
    def __init__(self, *largs):
        """
            class Vector -> no inheritance and it should be an array of floats !!!    
        """
        super(Vector,self).__init__(*largs)
        self.norm0 = math.sqrt(self[0]**2. + self[1]**2. + self[2]**2. )

    def __add__(self,other):
        """
            Adds an other Vector
        """
        return Vector(self[0]+other[0],self[1]+other[1],self[2]+other[2])

    def __mul__(self,alpha):
        """
            multiplies the vector by a real
        """
        return Vector(self[0]*alpha,self[1]*alpha,self[2]*alpha)

    def __rmul__(self,alpha):
        """
            right multiplication by a real
        """
        return Vector(self[0]*alpha,self[1]*alpha,self[2]*alpha)

    def fromPoints(self,A,B):
        """
            defines the Vector with 2 points, sets the norm
        """
        self[0]=B[0]-A[0]
        self[1]=B[1]-A[1]
        self[2]=B[2]-A[2]
        self.calcNorm0()

    def calcNorm0(self):
        """
            computes the norm of the vector
        """
        self.norm0=math.sqrt(self[0]**2+self[1]**2+self[2]**2)

    def set(self,x,y,z):
        """
            sets the coordinates of the Vector
        """
        self[0]=x
        self[1]=y
        self[2]=z
        self.calcNorm0()

    def normalize(self):
        """
            normalize the vector
        """
        self.calcNorm0()
        self[0]=self[0]/self.norm0
        self[1]=self[1]/self.norm0
        self[2]=self[2]/self.norm0
        self.calcNorm0()

    def draw(self,ax,p):
        """
            adds the point to a 3d existing matplotlib view, ax
        """
        pp = p.addi(Point([self[0], self[1], self[2]]))
        ax.plot([p[0],pp[0]],[p[1],pp[1]],[p[2],pp[2]],c='k')

    def projectOnCoordinatesSystem(self,u0,v0,w0):
        """
            returns the Vector given in the 3 Vectors u0,v0,w0 basis
        """
        cx=u0[0]*self[0]+u0[1]*self[1]+u0[2]*self[2]
        cy=v0[0]*self[0]+v0[1]*self[1]+v0[2]*self[2]
        cz=w0[0]*self[0]+w0[1]*self[1]+w0[2]*self[2]
        v=Vector((cx,cy,cz))
        return v



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
        
        d0=Vector()
        d0.fromPoints(cent0,p0)
        d1=Vector()
        d1.fromPoints(cent0,p1)
        d2=Vector()
        d2.fromPoints(cent0,p2)
        d3=Vector()
        d3.fromPoints(cent0,p3)
        

        n0 = cross(d3,d0,norm=False)
        n1 = cross(d0,d1,norm=False)
        n2 = cross(d1,d2,norm=False)
        n3 = cross(d2,d3,norm=False)
        
        self.ng=Vector()
        self.ng.set(0.25*(n0.x+n1.x+n2.x+n3.x),\
                    0.25*(n0.y+n1.y+n2.y+n3.y),\
                    0.25*(n0.z+n1.z+n2.z+n3.z))
        self.ng.normalize()

def cross(u,v, norm=True):
    pv=[]
    pv.append(u[1]*v[2]-u[2]*v[1])
    pv.append(u[2]*v[0]-u[0]*v[2])
    pv.append(u[0]*v[1]-u[1]*v[0])
    vp=Vector(pv)
    #print vp
    if norm:vp.normalize()
    return vp

def dot(u,v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]

def angle(u,v):
    angle = np.arctan2(cross(u,v,norm=False).norm0,dot(u,v))
    #
    if angle == 0.:return 0.
    # warning : only in 2D case
    if dot(cross(u,v), Vector((0.,0.,1.0)))>0.:
        return angle
    else :
    #    print "gne"
        return -angle

def unit_vector_like(vec):
    v = vec
    v.normalize()
    return v

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
    def __init__(self,p1,p2,p3):
        self.p0 = p1
        self.p1 = p2
        self.p2 = p3


    def computeNormal(self):    
        u = Vector((self.p1[0]-self.p0[0],\
			self.p1[1]-self.p0[1],\
			self.p1[2]-self.p0[2]))
        v = Vector((self.p2[0]-self.p1[0],\
			self.p2[1]-self.p1[1],\
			self.p2[2]-self.p1[2]))
        return cross(u,v,norm=True)


              
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
	xg=(self.p0.x+self.p1.x+self.p2.x)/3.
	yg=(self.p0.y+self.p1.y+self.p2.y)/3.
	zg=(self.p0.z+self.p1.z+self.p2.z)/3.
	pcg = Point((xg, yg, zg))
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
        

class Triangulation(dict):
    """
    dictionary
    * KEYS : 
        'vertices' : array of points, 3-real arrays (x,y,z)
        'faces'    : array of faces, 3 int arrays (3 index in vertices list)
    """
    def __init__(self, points = None, faces = None, **kwargs):
        self.vertices = points
        self.faces = faces

    def load_file(self,name):
        f = json.load(open(name,'r'))
        try:
            self.vertices = f['vertices']
            self.faces = f['faces']
        except:
            print 'not enough info in '+name
            return 0

    @property
    def cg(self):
        return np.mean(np.array(self.vertices)   ,axis=0)

    def _dict(self):
        d = {}
        d['vertices'] = self.vertices
        d['faces'] = self.faces
        return d


    def write_obj_file(self,name):
        if '.obj' in name :
            with open('./samples/'+name, 'w') as f:
                f.write('# OBJ file\n')
                for v in self.vertices:
                    f.write('v %.4f %.4f %.4f\n' % (v[0],v[1],v[2]))
                for t in self.faces:
                    f.write('f')
                    for p in t :
                        f.write(" %d" % p)
                    f.write('\n')

        elif '.json' in name :
            f = open('./samples/'+name, 'w')
            json.dump(self._dict(), f, indent =4)
    
    def reorient_convex(self):
        # reorient the faces of the icosahedron
        _new_faces = []
        for f in self.faces:
            p0 = Point(self.vertices[f[0]-1])
            p1 = Point(self.vertices[f[1]-1])
            p2 = Point(self.vertices[f[2]-1])
            tr = Triangle(p0, p1, p2)
            cg = tr.gravityCenter()
            n0 = tr.computeNormal()
            n1 = Vector(cg)
            if dot(n0,n1) < 0:
                fi = (f[0], f[2], f[1])
            else:
                fi = f
            _new_faces.append(fi)
        self.faces = _new_faces
        # ok the icosahedron is correctly oriented in the dict d

    def translate(self, vec):
        _vertices = []
        for p in self.vertices:
            _vertices.append([\
                    p[0] + vec[0],\
                    p[1] + vec[1],\
                    p[2] + vec[2]])
        self.vertices = _vertices

    def refine_on_sphere(self, radius = 1.) :
        d2 = {}
        _vertices = []
        _faces = []
        for f in self.faces:
            p0 = self.vertices[f[0]-1]
            p1 = self.vertices[f[1]-1]
            p2 = self.vertices[f[2]-1]
            i0 = None
            i1 = None
            i2 = None
            i3 = None
            i4 = None
            i5 = None
            if p0 not in _vertices:
                _vertices.append(p0) 
                i0 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p0:
                        i0 = i
                        break
            if p1 not in _vertices:
                _vertices.append(p1) 
                i1 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p1:
                        i1 = i
                        break
            if p2 not in _vertices:
                _vertices.append(p2)
                i2 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p2:
                        i2 = i
                        break
            p3 = [0.5*(p1[j]+p2[j]) for j in range(3)]  
            p4 = [0.5*(p2[j]+p0[j]) for j in range(3)]  
            p5 = [0.5*(p1[j]+p0[j]) for j in range(3)] 
            if p3 not in _vertices:
                _vertices.append(p3) 
                i3 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p3:
                        i3 = i
                        break
            if p4 not in _vertices:
                _vertices.append(p4) 
                i4 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p4:
                        i4 = i
                        break
            if p5 not in _vertices:
                _vertices.append(p5)
                i5 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p5:
                        i5 = i
                        break
            i0 += 1 
            i1 += 1 
            i2 += 1 
            i3 += 1 
            i4 += 1 
            i5 += 1
            _faces.append((i4, i3, i2))
            _faces.append((i5, i1, i3))
            _faces.append((i5, i3, i4))
            _faces.append((i0, i5, i4))
        _vertices2 = []
        for p in _vertices:
            d = np.sqrt( (p[0] - self.cg[0])**2. +\
                         (p[1] - self.cg[1])**2. +\
                         (p[2] - self.cg[2])**2. )
            ps =[ radius * (p[i] - self.cg[i])/d + self.cg[i] for i in range(3) ]
            _vertices2.append(ps)
        self.vertices = _vertices2
        self.faces = _faces


class Sphere(Triangulation):
    def __init__(self, radius = 1., center = (0.,0.,0.), refin = 1 ):
        d = Triangulation()
        d.load_file('./samples/icosahedron.json')
        d.translate(center)
        for i in range(refin):
            d.refine_on_sphere(radius)
        self.vertices = d.vertices
        self.faces = d.faces
        self.plane_pos = []
        self.radius = float(radius)
        plt.clf()
        for j,p in enumerate(self.vertices):
            adipos = [(p[i] - self.cg[i])/radius for i in range(3) ]
            lat = np.arctan2(adipos[1], adipos[0])
            if abs(adipos[2]) < 1. : 
                lon = np.arccos(adipos[2])
            else:
                lon = 0.
            self.plane_pos.append((lat,lon))
            #plt.scatter(lat,lon,c= 'k', s = 1)
        #plt.show()






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
    
