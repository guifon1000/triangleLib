from Triangulation import Triangulation
from Point import Point
from Vector import Vector
import numpy as np



class Sphere(Triangulation):
    def __init__(self, radius = 1., center = (0.,0.,0.), refin = 1 , **kwargs):
        d = Triangulation()
        d.load_file('../samples/icosahedron.json')
        d.translate(center)
        for i in range(refin):
            d = d.refine_2()
        _vertices = []
        for p in d.vertices:
            dist = np.sqrt(np.sum([ (p[i] - d.cg[i])**2. for i in range(3)]))
            _vertices.append(Point([d.cg[i] + (p[i]-d.cg[i])*radius/dist for i in range(3)]))
        self.vertices = _vertices
        self.faces = d.faces
        self.plane_pos = []
        self.radius = float(radius)
