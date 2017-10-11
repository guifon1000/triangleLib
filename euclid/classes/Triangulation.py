import sys
sys.path.append('../')
from functions import dot
import json
import numpy as np
from Point import Point
from Vector import Vector
from Triangle import Triangle

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

 
    def write_obj_file(self,name):
        if '.obj' in name :
            with open('../generated_geometries/'+name, 'w') as f:
                f.write('# OBJ file\n')
                for v in self.vertices:
                    f.write('v %.4f %.4f %.4f\n' % (v[0],v[1],v[2]))
                for t in self.faces:
                    f.write('f')
                    for p in t :
                        f.write(" %d" % p)
                    f.write('\n')

        elif '.json' in name :
            f = open('../generated_geometries/'+name, 'w')
            json.dump(self._dict(), f, indent =4)
   
    @property
    def triangle_list(self):
        _triangles = []
        for f in self.faces:
            p0 = Point(self.vertices[f[0]-1])
            p1 = Point(self.vertices[f[1]-1])
            p2 = Point(self.vertices[f[2]-1])
            tr = Triangle((p0, p1, p2))
            _triangles.append(tr)
        return _triangles

     
    def reorient_convex(self):
        # reorient the faces IF CONVEX
        _new_faces = []
        for f in self.faces:
            p0 = Point(self.vertices[f[0]-1])
            p1 = Point(self.vertices[f[1]-1])
            p2 = Point(self.vertices[f[2]-1])
            tr = Triangle((p0, p1, p2))
            n0 = tr.normal.unit()
            n1 = Vector(tr.cg)
            if dot(n0,n1) < 0:
                fi = (f[0], f[2], f[1])
            else:
                fi = f
            _new_faces.append(fi)
        self.faces = _new_faces
        # ok the convex is correctly oriented in the dict d

    def translate(self, vec):
        _vertices = []
        for p in self.vertices:
            _vertices.append([\
                    p[0] + vec[0],\
                    p[1] + vec[1],\
                    p[2] + vec[2]])
        self.vertices = _vertices




    def refine_2(self) :
        d2 = {}
        _vertices = []
        _faces = []



        for f in self.faces:
            p0 = self.vertices[f[0]-1]
            p1 = self.vertices[f[1]-1]
            p2 = self.vertices[f[2]-1]
            if p0 in _vertices:
                i0 = _vertices.index(p0)
            else:
                _vertices.append(p0)
                i0 = len(_vertices)-1
            if p1 in _vertices:
                i1 = _vertices.index(p1)
            else:
                _vertices.append(p1)
                i1 = len(_vertices)-1
            if p2 in _vertices:
                i2 = _vertices.index(p2)
            else:
                _vertices.append(p2)
                i2 = len(_vertices)-1


            p3 = Point([0.5*(p1[j]+p2[j]) for j in range(3)])
            p4 = Point([0.5*(p2[j]+p0[j]) for j in range(3)])
            p5 = Point([0.5*(p1[j]+p0[j]) for j in range(3)])

            if p3 in _vertices:
                i3 = _vertices.index(p3)
            else:
                _vertices.append(p3)
                i3 = len(_vertices)-1
            if p4 in _vertices:
                i4 = _vertices.index(p4)
            else:
                _vertices.append(p4)
                i4 = len(_vertices)-1
            if p5 in _vertices:
                i5 = _vertices.index(p5)
            else:
                _vertices.append(p5)
                i5 = len(_vertices)-1



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
        return Triangulation(points = _vertices, faces = _faces)


    
