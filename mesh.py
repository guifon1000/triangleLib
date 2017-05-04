import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
import bidLib as l2d
"""
    mesh lib 
"""

def write_geo(name,geom):
    fg = open(name+'.geo','w')
    for l in geom.get_code():
        fg.write(l)
    fg.close()


class Polyline2D(list):    #always closed
    def __init__(self, *largs,**kwargs):
        super(Polyline2D,self).__init__(*largs)
        if kwargs.has_key('z'):
            self.z=kwargs['z']
        else:
            self.z = 0.
        self.pt3d = []
        for i in range(len(self)):
            p = l2d.point(self[i],z=self.z)
            self.pt3d.append(p)

    def extrude_segments(self,ep):
        pass

    def box_2d(self,**kwargs):
        """
          3  <------------- 2
          |                 ^
          v                 |
          0  -------------> 1
        """

        geom = pg.Geometry()
        try:
            facdom = kwargs['facdom']
        except:
            facdom=0.5
        try:
            self.trax = (0.,0.,kwargs['thick'])
        except:
            self.trax = (0.,0.,0.1)
        try:
            self.name = kwargs['name']
        except:
            print('I need a \'name\' parameter !')
            return 
        lcar = 1.0#2.*self.trax[2]

        lineloop = []
        pbox = []
        p0 = l2d.point((-facdom,-facdom),z=self.z)
        p1 = l2d.point((facdom,-facdom),z=self.z) 
        p2 = l2d.point((facdom,facdom),z=self.z) 
        p3 = l2d.point((-facdom,facdom),z=self.z) 
        p0 = geom.add_point(p0,lcar)
        p1 = geom.add_point(p1,lcar)
        p2 = geom.add_point(p2,lcar)
        p3 = geom.add_point(p3,lcar)
        pbox.append(p0)
        pbox.append(p1)
        pbox.append(p2)
        pbox.append(p3)
        south = geom.add_line(p0,p1)
        east = geom.add_line(p1,p2) 
        north = geom.add_line(p2,p3) 
        west = geom.add_line(p3,p0) 
        lineloop.append(south)
        lineloop.append(east)
        lineloop.append(north)
        lineloop.append(west)
        ext = geom.extrude(south,translation_axis=self.trax)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
        geom.add_physical_surface(ext[1],label='Ymin')
        ext = geom.extrude(east,translation_axis=self.trax)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
        geom.add_physical_surface(ext[1],label='Xmax')
        ext = geom.extrude(north,translation_axis=self.trax)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
        geom.add_physical_surface(ext[1],label='Ymax')
        ext = geom.extrude(west,translation_axis=self.trax)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
        geom.add_physical_surface(ext[1],label='Xmin')
        print '---------------------------------'
        opl = []
        lobj = []
        for p in self.pt3d:
            pp = geom.add_point(p,1.)
            opl.append(pp)
        for i,p in enumerate(opl[:-1]):
            li = geom.add_line(opl[i],opl[i+1])
            lobj.append(li)
            
            lineloop.append(li)
        #close the polyline
        li = geom.add_line(opl[-1],opl[0])
        lobj.append(li)
        lineloop.append(li)


        
        sfobj = []
        for li in lobj:
            ext = geom.extrude(li,translation_axis=self.trax)
            #geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
            sfobj.append(ext[1])
        geom.add_physical_surface(sfobj,label='cylinder')


        lloop = geom.add_line_loop(lineloop)
        sf0 = geom.add_plane_surface(lloop)

        #geom.add_physical_surface(sf0,label='Zmin')


        ll2 = []
        opl2 = []
        bpl2 = []
        for p in pbox:
            p2 = l2d.point((p.x[0],p.x[1]),z=self.z+self.trax[2])
            p2 = geom.add_point(p2,lcar)
            bpl2.append(p2)

        for i,p2 in enumerate(bpl2[:-1]):
            l2 = geom.add_line(p2,bpl2[i+1])
            ll2.append(l2)

        l2 = geom.add_line(bpl2[-1],bpl2[0])
        ll2.append(l2)


        for p in opl:
            p2 = l2d.point((p.x[0],p.x[1]),z=self.z+self.trax[2])
            p2 = geom.add_point(p2,lcar)
            opl2.append(p2)


        for i,p2 in enumerate(opl2[:-1]):
            l2 = geom.add_line(p2,opl2[i+1])
            ll2.append(l2)

        l2 = geom.add_line(opl2[-1],opl2[0])
        ll2.append(l2)

        #for p in pbox:
            #p2 = l2d.point((p.x[0],p.x[1]),z=self.z+self.trax[2])

        lloop2 = geom.add_line_loop(ll2)

        sf1 = geom.add_plane_surface(lloop2)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % sf1.id )

        #geom.add_physical_surface(sf1,label='Zmax')
              


        write_geo(self.name,geom)
        print('2D Geometry successfully created')
        #points,cells,pt_data, cell_data, field_data = pg.generate_mesh(geom)
        #meshio.write('volumic0.vtu', points, cells)# ells['quad'])    #{'mu':mu}




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
