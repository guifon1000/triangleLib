
from scipy import interpolate
from Vector import Vector
from Point import Point
import sys 
sys.path.append('../')
from functions import cross

class Frame(list):
    def __init__(self, s, tck, **kwargs):
        super(Frame,self).__init__()
        if kwargs.has_key('type'):
            t =  interpolate.splev(s , tck, der = 0)
            orig = Point([t[0], t[1], t[2]])
            self.append(orig)
            if kwargs['type'] == 'frenet':
                t =  interpolate.splev(s , tck, der = 1)
                xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                t =  interpolate.splev(s , tck, der = 2)
                xpp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                T = xp.unit()
                self.append(T)
                B = cross(xp, xpp)
                self.append(B)
                N = cross(B, T)
                self.append(N)

            if kwargs['type'] == 'Xnat':
                t =  interpolate.splev(s , tck, der = 1)
                xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                T = xp.unit()
                self.append(T)
                B = cross((1.,0.,0.), T)
                self.append(B)
                N = cross(T, B)
                self.append(N)


