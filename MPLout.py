import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt




def plot3(MSH):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for f in MSH.pts:
        ax.scatter(f.x, f.y, f.z,marker='s')

    for f in MSH.quads:
        s0 = [f.p0.x,f.p1.x,f.p2.x,f.p3.x,f.p0.x]
        s1 = [f.p0.y,f.p1.y,f.p2.y,f.p3.y,f.p0.y]
        s2 = [f.p0.z,f.p1.z,f.p2.z,f.p3.z,f.p0.z]
        ax.plot(s0,s1,s2,'k')


    for f in MSH.faces:
        p0=f[0]
        p1=f[1]
        p2=f[2]


        s0 = [MSH.pts[p0].x,MSH.pts[p1].x,MSH.pts[p2].x]
        s1 = [MSH.pts[p0].y,MSH.pts[p1].y,MSH.pts[p2].y]
        s2 = [MSH.pts[p0].z,MSH.pts[p1].z,MSH.pts[p2].z]
        ax.plot(s0,s1,s2,'k')
    plt.show()
