
from Point import Point
class Polyline2D(list):    #always closed
    def __init__(self, *largs,**kwargs):
        super(Polyline2D,self).__init__(*largs)
        if kwargs.has_key('z'):
            self.z=kwargs['z']
        else:
            self.z = 0.
        if kwargs.has_key('closed'):
            self.check_closed = kwargs['closed']
        else : 
            self.check_closed = False
        if self.check_closed : 
            p = self[0]
            self.append(p)
        self.pt3d = []
        for i in range(len(self)):
            p = Point(self[i])
            self.pt3d.append( Point([p[0] , p[1], 0.]))


    def to_frame(self, f, **kwargs):
        try:
            fac_scale = kwargs['scale']
        except:
            fac_scale = 1.
        for i in range(len(self.pt3d)) :
            self.pt3d[i][0] *= fac_scale * (f[2][0] + f[3][0] )
            self.pt3d[i][0] += f[0][0] 
            self.pt3d[i][1] *= fac_scale * (f[2][1] + f[3][1] )
            self.pt3d[i][1] += f[0][1] 
            self.pt3d[i][2] *= fac_scale * (f[2][2] + f[3][2] )
            self.pt3d[i][2] += f[0][2] 

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
