from kivy3 import Vector3
from kivy3.core.geometry import Geometry
from kivy3.core.face3 import Face3
import pygmsh as pg
import triangleLib as tl


class Triangle(Geometry):

    _tri_vertices = [(-1, 1, -1), (1, 1, -1),
                  (1, -1, -1)]

    _tri_faces = [(0, 1, 2)]

    _tri_normals = [(0, 0, 1)]



    def __init__(self, width, height, depth, **kw):
        name = kw.pop('name', '')
        super(Triangle, self).__init__(name)
        self.width_segment = kw.pop('width_segment', 1)
        self.height_segment = kw.pop('height_segment', 1)
        self.depth_segment = kw.pop('depth_segment', 1)

        self.w = width
        self.h = height
        self.d = depth

        self._build_tri()

    def _build_tri(self):

        for v in self._tri_vertices:
            v = Vector3(0.5 * v[0] * self.w,
                        0.5 * v[1] * self.h,
                        0.5 * v[2] * self.d)
            self.vertices.append(v)

        n_idx = 0
        for f in self._tri_faces:
            face3 = Face3(*f)
            normal = self._tri_normals[n_idx / 2]
            face3.vertex_normals = [normal, normal, normal]
            n_idx += 1
            self.faces.append(face3)



class Msh(Geometry):
    def __init__(self, fi, **kwargs):
        name = kwargs.pop('name', '')
        super(Msh, self).__init__(name)
        self.width_segment = kwargs.pop('width_segment', 1)
        self.height_segment = kwargs.pop('height_segment', 1)
        self.depth_segment = kwargs.pop('depth_segment', 1)
        f = open(fi,'r').readlines()
        _vertices = []
        app = False
        ist = 0
        nv = None
        while app == False :
            l = f[ist] 
	    if l.startswith('$Nodes'):
                nv = int(f[ist+1])
                app = True
            ist+=1
        print 'there are '+str(nv)+' vertices strating at line '+str(ist+1)+' :'+f[ist+1]
        print 'ending at line '+str(ist+nv)+' : '+f[ist+nv]
        for i in range(ist+1,ist+nv+1):
            l=f[i].split()
            _vertices.append((float(l[1]), float(l[2]), float(l[3])))
        self.msh_vertices = _vertices
        _faces = []
        for i in range(ist+nv+1,len(f)):
            lp = f[i].split()
            if (len(lp) > 3) and (lp[1] == '2') :
                    _faces.append((int(lp[5]),\
                                      int(lp[6]),\
                                      int(lp[7])))
        self.msh_faces = _faces
        _normals = []
        for fa in self.msh_faces :
            p0 = tl.Point(*self.msh_vertices[fa[0]-1])
            p1 = tl.Point(*self.msh_vertices[fa[1]-1])
            p2 = tl.Point(*self.msh_vertices[fa[2]-1])
            v1 = tl.Vector(float(p1[0])-float(p0[0]),float(p1[1])-float(p0[1]),float(p1[2])-float(p0[2]))
            v2 = tl.Vector(float(p2[0])-float(p1[0]),float(p2[1])-float(p1[1]),float(p2[2])-float(p1[2]))
            n = tl.cross(v1,v2,norm=False)
            _normals.append((n[0],n[1],n[2]))
        self.msh_normals = _normals
         
        self._build_msh()


    def _build_msh(self):
        for v in self.msh_vertices:
            v0 = Vector3(v[0]*0.4 ,
                        v[1]*0.4 ,
                        v[2]*0.4 )
            self.vertices.append((v0[0],v0[1],v0[2]))
        print len(self.vertices)
        n_idx = 0
        for f in self.msh_faces:
            face3 = Face3(*f)
            normal = self.msh_normals[n_idx]
            #face3.vertex_normals = [normal, normal, normal]
            face3.normal = normal
            n_idx += 1
            self.faces.append(face3)
 
