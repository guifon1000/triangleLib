import numpy as np
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import triangleLib as tl

class Frame(list):
    def __init__(self, s, tck, **kwargs):
        super(Frame,self).__init__()
        if kwargs.has_key('type'):
            t =  interpolate.splev(s , tck, der = 0)
            orig = tl.Point([t[0], t[1], t[2]])
            self.append(orig)
            if kwargs['type'] == 'frenet':
                t =  interpolate.splev(s , tck, der = 1)
                xp  = tl.Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                t =  interpolate.splev(s , tck, der = 2)
                xpp  = tl.Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                T = tl.unit_vector_like(xp)
                self.append(T)
                B = tl.cross(xp,xpp, norm = True)
                self.append(B)
                N = tl.cross(B,T, norm = True )
                self.append(N)

            if kwargs['type'] == 'Xnat':
                t =  interpolate.splev(s , tck, der = 1)
                xp  = tl.Vector( [float(t[0]), float(t[1]), float(t[2]) ])
                T = tl.unit_vector_like(xp)
                self.append(T)
                B = tl.cross((1.,0.,0.), T, norm = True)
                self.append(B)
                N = tl.cross(T,B, norm = True )
                self.append(N)


if __name__ == '__main__':
    x = [0.0, 0.1, 0.4, 0.6, 0.9]
    y = [0.0, 0.05, 0.09, 0.11, 0.4]
    z = [0.0, 0.9, 1.8, 5.2, 7.9]

    tck, u = interpolate.splprep([x,y,z], s=2)
    #x_knots, y_knots, z_knots = interpolate.splev(tck[0], tck)
    t = np.linspace(0, 1, 100)
    fig2 = plt.figure(2)
    ax3d = fig2.add_subplot(111, projection='3d')
    for s in t:
        f = Frame(s, tck, type = 'Xnat')
        f[1].draw(ax3d, f[0])
        f[2].draw(ax3d, f[0])
        f[3].draw(ax3d, f[0])

    plt.axis('equal')
    plt.show()

