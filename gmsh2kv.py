from kivy3 import Vector3
from kivy3.core.geometry import Geometry
from kivy3.core.face3 import Face3
import pygmsh as pg
import triangleLib as tl

class Msh(Geometry):
    def __init__(self, **kwargs):
        tri0 = tl.Triangulation(file = kwargs['file'])
        name = kwargs.pop('name', '')
        super(Msh, self).__init__(name)
        self.width_segment = kwargs.pop('width_segment', 1)
        self.height_segment = kwargs.pop('height_segment', 1)
        self.depth_segment = kwargs.pop('depth_segment', 1)
        self._vertices = tri0.vertices
        self._faces = tri0.faces
        self._normals = tri0.normals
        self._build_msh()


    def _build_msh(self):
        for v in self._vertices:
            v0 = Vector3(v[0],
                        v[1],
                        v[2] )
            self.vertices.append((v0[0],v0[1],v0[2]))
        n_idx = 0
        for f in self._faces:
            face3 = Face3(*f)
            normal = self._normals[n_idx]
            face3.vertex_normals = [normal, normal, normal]
            face3.normal = normal
            n_idx += 1
            self.faces.append(face3)
 
