import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
sys.path.append('./VortexLatticeMethod/')
import vlm as vlm


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

class surfaceMesh(object):
    def __init__(self,points,Nu,Nv):
        self.Nu = Nu
        self.Nv = Nv

        geom = pg.Geometry()
        MSHpts=[]
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
            p2 = vlm.fluidPoint(p[0],p[1],p[2])
            p2.computeVelocity(elem)
            fp.append(p2)

        self.d={'ux' : np.array([p.vel[0] for p in fp]),\
                'uy' : np.array([p.vel[1] for p in fp]),\
                'uz' : np.array([p.vel[2] for p in fp])}
        meshio.write('volumic1.vtu', self.points, self.cells, point_data=self.d)
