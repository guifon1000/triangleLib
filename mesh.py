import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
sys.path.append('./profiles/')
import splineProfileMultiParam as prf

"""
    mesh lib 
"""

def mesh_cylinder2D():
    geom = pg.Geometry()
    diam = 1.0
    domfac = 10. 
    cyl = geom.add_circle((0,0,0),diam,0.04,num_sections = 3)
    p0 = geom.add_point((-domfac*diam,domfac*diam,0.),0.1)
    p1 = geom.add_point((-domfac*diam,-domfac*diam,0.),0.1)
    p2 = geom.add_point((2.*domfac*diam,-domfac*diam,0.),0.1)
    p3 = geom.add_point((2.*domfac*diam,domfac*diam,0.),0.1)
    transfinite_frontiers = 10
    west = geom.add_line(p0,p1)
    geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (west,transfinite_frontiers) )
    south = geom.add_line(p1,p2)
    geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (south,transfinite_frontiers) )
    east = geom.add_line(p2,p3)
    geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (east,transfinite_frontiers) )
    north = geom.add_line(p3,p0)
    geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (north,transfinite_frontiers) )
    lns = []
    lns.append(west)
    lns.append(south)
    lns.append(east)
    lns.append(north)
    we = geom.add_physical_line(west)
    ea = geom.add_physical_line(east)
    no = geom.add_physical_line(north)
    so = geom.add_physical_line(south)
    for l in cyl:
        lns.append(l)
        #li = geom.add_line(l)
        geom.add_physical_line(l)
    ll = geom.add_line_loop(lns)
    sf0 = geom.add_plane_surface(ll)
    #ext = geom.extrude("Surface{"+sf0+"}",translation_axis=(0.,0.,1.))
    psf = geom.add_physical_surface(sf0)
    points, cells = pg.generate_mesh(geom)

    epex = 1.
    ls = cells['line']
    ge2 = pg.Geometry()
    for s in ls:
        seg = [int(st) for st in s]
        p0 = points[seg[0]]
        p1 = points[seg[1]]
        p2 = [p1[0],p1[1],p1[2]+epex]
        p3 = [p0[0],p0[1],p0[2]+epex]
        ge2.add_point(p0,0.1)
        ge2.add_point(p1,0.1)
        ge2.add_point(p2,0.1)
        ge2.add_point(p3,0.1)

    f2  = open('test33.geo','w')
    for l in ge2.get_code():
        f2.write(l)
    f2.close()
    print 1/0


    X = np.array([p[0] for p in points])
    Y = np.array([p[1] for p in points])
    Z = np.array([p[2] for p in points])
    tris = cells['triangle']
    plt.clf()
    print cells
    for p in points:plt.scatter(p[0],p[1])
    for t in tris:
        p0 = points[t[0]]
        p1 = points[t[1]]
        p2 = points[t[2]]
        plt.plot((p0[0],p1[0]),(p0[1],p1[1]))
        plt.plot((p1[0],p2[0]),(p1[1],p2[1]))
        plt.plot((p2[0],p0[0]),(p2[1],p0[1]))
    plt.axis('equal')
    plt.show()

    #self.points, self.cells = pg.generate_mesh(geom)




class Triangulation(list):
    def __init__(self,stl = None):
        print "triangulation creation"
        if stl is not None : 
            geom = pg.Geometry()
            geom._GMSH_CODE.append('Merge \''+stl+'\';\n')
            self.points, self.cells = pg.generate_mesh(geom)


        fg = open('test0.geo','w')
        for l in geom.get_code():
            fg.write(l)
        self.points, self.cells = pg.generate_mesh(geom)

        self.X = np.array([p[0] for p in self.points])
        self.Y = np.array([p[1] for p in self.points])
        self.Z = np.array([p[2] for p in self.points])
        self.tris = self.cells['triangle']



class thickWing(object):
    def __init__(self,pf,name):
        geom = pg.Geometry()
        ptsExtra = []
        lns = []
        for (i,x) in enumerate(pf.x):
            p = geom.add_point((x,pf.extra[i],1.),0.)
            ptsExtra.append(p)
        for i in range(len(ptsExtra)-1):
            l = geom.add_line(ptsExtra[i],ptsExtra[i+1])
            lns.append(l)
        ptsIntra = []
        for (i,x) in enumerate(pf.x):
            p = geom.add_point((x,pf.intra[i],1.),0.)
            ptsIntra.append(p)
        for i in range(len(ptsIntra)-1):
            l = geom.add_line(ptsIntra[i],ptsIntra[i+1])
            lns.append(l)
        lnTe = geom.add_line(ptsExtra[-1],ptsIntra[-1])
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (lnTe,2) )
        lns.append(lnTe)
        lloop = geom.add_line_loop(lns)
        sfStart = geom.add_ruled_surface(lloop) 
        geom.add_physical_surface(sfStart) 
                #geom._GMSH_CODE.append('Recombine Surface {%s};' % sfIntra)
        fg = open(name+'.geo','w')
        for l in geom.get_code():
            fg.write(l)


def loadXfoilFile(name):
    f = open(name,'r').readlines()
    points = []
    for i in range(1,len(f)):
        l = f[i].split()
        p = [float(l[0]),float(l[1])]
        points.append(p)
    panels = []
    for i in range(len(points)-1):
        p0 = points[i]
        p1 = points[i+1]
        panels.append([p0,p1])
    return points,panels



class meshAround2dObject(object):
    def __init__(self,name):
        points,panels = loadXfoilFile(name) 
        geom = pg.Geometry()
        for i,p in enumerate(panels):
            p0 = p[0]
            p1 = p[1]
            p_0 = geom.add_point((p0[0],p0[1],0.),0.)
            p_1 = geom.add_point((p1[0],p1[1],0.),0.)
            l = geom.add_line(p_0,p_1)
            plt.scatter([p0[0],p1[0]],[p0[1],p1[1]])
        fg = open(name+'2D.geo','w')
        for l in geom.get_code():
            fg.write(l)
        plt.axis('equal')
        plt.show() 
        






class surfaceMesh(object):
    def __init__(self,points,Nu,Nv):
        self.Nu = Nu
        self.Nv = Nv
        MSHpts = []
        geom = pg.Geometry()
        for i,p in enumerate(points):
            p = geom.add_point([p.x,p.y,p.z],lcar=0.1)
            MSHpts.append(p)

        l0 = geom.add_line(MSHpts[0],MSHpts[1])
        l1 = geom.add_line(MSHpts[1],MSHpts[2])
        l2 = geom.add_line(MSHpts[2],MSHpts[3])
        l3 = geom.add_line(MSHpts[3],MSHpts[0])



        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l0,Nu) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l1,Nv) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l2,Nu) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l3,Nv) )
        ll = geom.add_line_loop((l0,l1,l2,l3))
        sf = geom.add_ruled_surface(ll)
        geom._GMSH_CODE.append('Transfinite Surface {%s};' % sf)
        geom.add_physical_surface(sf, label='zob')

        for l in geom.get_code().split('\n'):
            print l


        geom._GMSH_CODE.append('Recombine Surface {%s};' % sf)

        fg = open('test0.geo','w')
        for l in geom.get_code():
            fg.write(l)
        self.points, self.cells = pg.generate_mesh(geom)
         
        self.X = np.array([p[0] for p in self.points])
        self.Y = np.array([p[1] for p in self.points])
        self.Z = np.array([p[2] for p in self.points])
        self.quads = self.cells['quad']
        
    def MPLout(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for t in self.cells['quad']:
            p0=self.points[t[0]]
            p1=self.points[t[1]]
            p2=self.points[t[2]]
            p3=self.points[t[3]]
            ax.plot([p0[0],p1[0],p2[0],p3[0],p0[0]],\
                [p0[1],p1[1],p2[1],p3[1],p0[1]],\
                [p0[2],p1[2],p2[2],p3[2],p0[2]],'k')
        plt.show()


class volumicMesh(object):
    def __init__(self, Lx, Ly, Lz, N1, N2, N3, center=[0.,0.,0.]):
        self.N1 = N1
        self.N2 = N2
        self.N3 = N3
        geom = pg.Geometry()
        p0 = geom.add_point([center[0]-Lx/2.,\
                            center[1]+Ly/2.,\
                            center[2]-Lz/2.],\
                            lcar=0.1)


        p1 = geom.add_point([center[0]+Lx/2.,\
                            center[1]+Ly/2.,\
                            center[2]-Lz/2.],\
                            lcar=0.1)


        p2 = geom.add_point([center[0]+Lx/2.,\
                            center[1]-Ly/2.,\
                            center[2]-Lz/2.],\
                            lcar=0.1)


        p3 = geom.add_point([center[0]-Lx/2.,\
                            center[1]-Ly/2.,\
                            center[2]-Lz/2.],\
                            lcar=0.1)


        l0 = geom.add_line(p0,p1)
        l1 = geom.add_line(p1,p2)
        l2 = geom.add_line(p2,p3)
        l3 = geom.add_line(p3,p0)



        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l0,N1) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l1,N2) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l2,N1) )
        geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l3,N2) )

        ll = geom.add_line_loop((l0,l1,l2,l3))
        sf = geom.add_ruled_surface(ll)
        geom._GMSH_CODE.append('Transfinite Surface {%s};' % sf)
        geom._GMSH_CODE.append('Recombine Surface {%s};' % sf)
        geom._GMSH_CODE.append('vol[] = Extrude {0,0,%s} { Surface {%s};Layers{%s};Recombine;};' % (Lz,sf,N3))
        geom.add_physical_volume('vol',label='zobi')

        fg = open('volMesh.geo','w')
        for l in geom.get_code():
            fg.write(l)
        fg.close()
        self.points, self.cells = pg.generate_mesh(geom)
        self.hex=self.cells['hexahedron']
        meshio.write('volumic0.vtu', self.points, self.cells)# ells['quad'])    #{'mu':mu}

    def compute(self,elem):
        fp=[]
        for p in self.points:
            #p2 = vlm.fluidPoint(p[0],p[1],p[2])
            p2.computeVelocity(elem)
            fp.append(p2)

        self.d={'ux' : np.array([p.vel[0] for p in fp]),\
                'uy' : np.array([p.vel[1] for p in fp]),\
                'uz' : np.array([p.vel[2] for p in fp])}
        meshio.write('volumic1.vtu', self.points, self.cells, point_data=self.d)
