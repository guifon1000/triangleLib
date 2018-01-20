import numpy as np
from Point import Point
from Polyline3D import Polyline3D


class Polyline2D(list):    #always closed
    def __init__(self, *largs,**kwargs):
        super(Polyline2D,self).__init__(*largs)
        if kwargs.has_key('closed'):
            if kwargs['closed'] : 
                if not self.is_closed:
                    self.append(self[0])
        if kwargs.has_key('z'):
            self.z=kwargs['z']
        else:
            self.z = 0.
        self.pt3d = []
        for i in range(len(self)):
            p = Point(self[i])
            self.pt3d.append( Point([0. , 0., 0.]))


    @property
    def is_closed(self):
        if self[0] == self[-1]:
            return True
        else:
            return False

    def close(self):
        if self.is_closed:
            print 'the polyline is already closed'
        else:
            self.append(self[0])

    def to_frame(self, frame, **kwargs):
        try:
            fac_scale = kwargs['scale']
        except:
            fac_scale = 1.

        pol3d = []
        for p in self:
            loc  =  np.dot([0., p[0], p[1]], frame[1])
            pol3d.append(Point([frame[0][i] + fac_scale * loc[i] for i in range(3) ]))
        return Polyline3D(pol3d)


    def pop_to_geom(self, geom):
        pts = []
        lns = []
        pol = self.pt3d
        if self.check_closed:
            for i,p in enumerate(pol[:-1]):
                p = geom.add_point(p,0.1)
                pts.append(p)
            for i in range(len(pts)-1):
                l = geom.add_line(pts[i],pts[i+1])
                lns.append(l)
           
            l = geom.add_line(pts[-1], pts[0])
            lns.append(l)
        return lns
