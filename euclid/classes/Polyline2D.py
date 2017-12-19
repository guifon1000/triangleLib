import numpy as np
from Point import Point



class Polyline2D(list):    #always closed
    def __init__(self, *largs,**kwargs):
        super(Polyline2D,self).__init__(*largs)
        if kwargs.has_key('closed'):
            if kwargs['closed'] : self.append(self[0])
        if kwargs.has_key('z'):
            self.z=kwargs['z']
        else:
            self.z = 0.
        self.pt3d = []
        if self[0] == self[-1] :
            self.check_closed = True
        else : 
            self.check_closed = False
        for i in range(len(self)):
            p = Point(self[i])
            self.pt3d.append( Point([0. , 0., 0.]))


    def to_frame(self, frame, **kwargs):
        try:
            fac_scale = kwargs['scale']
            fac_scale = 1.
        except:
            fac_scale = 1.


        for i,p in enumerate(self):
            loc  =  np.dot([0., p[0], p[1]], frame[1])
            self.pt3d[i][0] = frame[0][0] + loc[0]
            self.pt3d[i][1] = frame[0][1] + loc[1]
            self.pt3d[i][2] = frame[0][2] + loc[2]



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
