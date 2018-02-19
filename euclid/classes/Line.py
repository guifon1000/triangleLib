from classes.Point import Point
from classes.Vector import Vector
import sys
sys.path.append('../')
from functions import cross
class Line(list):
    """
    defines a line 
    """
    def __init__(self, *largs):
        if ([type(largs[0][i]) for i in range(len(*largs))] == [Point,Vector]):
            pt = largs[0][0]
            vec = largs[0][1].unit()
            super(Line, self).__init__((pt, vec))

    def parameter_point(self, par):
        pt = [self[0][i] + par * self[1][i] for i in range(3) ]
        return Point(pt)
