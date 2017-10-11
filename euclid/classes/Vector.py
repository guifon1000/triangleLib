import numpy as np
import pygmsh as pg


class Vector(list):
    def __init__(self, *largs):
        super(Vector,self).__init__(*largs)


    @property
    def norm(self):
        return np.sqrt(self.norm2)

    @property
    def norm2(self):
        return np.sum([c**2. for c in self])

    def unit(self):
        return Vector( [ self[i]/self.norm for i in range(3)])

