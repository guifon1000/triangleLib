from classes.Point import Point
from classes.Vector import Vector
from classes.Line import Line
import sys
sys.path.append('../')
from functions import cross
class Line(list):
    """
    defines a segment 
    """
    def __init__(self, *largs):
        if ([type(largs[0][i]) for i in range(len(*largs))] == [Point, Point]):
            super(Segment, self).__init__(*largs)
