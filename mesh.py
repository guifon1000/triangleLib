import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
sys.path.append('./profiles/')
import splineProfileMultiParam as prf
from scipy.interpolate import splev

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
        span = 0.5
        corde = 0.05
        fac = 0.0
        beta = 0.
        n=5
        forwardT = 0.75
        spanWing = np.linspace(0.,span,num=n,endpoint=False)
        def prandtl(x,R,corde,fac,beta):
            expo = np.exp(-fac*(1.-x/R)/((x/R)*np.cos(beta)))
            return corde*(2./np.pi)*np.arccos(expo)
        chordSpan = prandtl(spanWing,span,corde,fac,beta)
        geom = pg.Geometry()
        extrados = []
        intrados = []
        tip = []
        base = []
        for ispan in range(len(chordSpan)-1):
            # we work between ispan and ispan + 1
            ptsExtra0 = []
            ptsExtra1 = []
            ptsIntra0 = []
            ptsIntra1 = []
            for (i,x) in enumerate(pf.x):
                p = geom.add_point(((x-forwardT)*chordSpan[ispan],pf.extra[i]*chordSpan[ispan],spanWing[ispan]),0.)
                ptsExtra0.append(p)
                p = geom.add_point(((x-forwardT)*chordSpan[ispan+1],pf.extra[i]*chordSpan[ispan+1],spanWing[ispan+1]),0.)
                ptsExtra1.append(p)
                p = geom.add_point(((x-forwardT)*chordSpan[ispan],pf.intra[i]*chordSpan[ispan],spanWing[ispan]),0.)
                ptsIntra0.append(p)
                p = geom.add_point(((x-forwardT)*chordSpan[ispan+1],pf.intra[i]*chordSpan[ispan+1],spanWing[ispan+1]),0.)
                ptsIntra1.append(p)
            for i in range(len(ptsExtra0)-1):
                l0 = geom.add_line(ptsExtra0[i+1],ptsExtra0[i])
                l1 = geom.add_line(ptsExtra1[i],ptsExtra1[i+1])
                lf = geom.add_line(ptsExtra0[i],ptsExtra1[i])
                lr = geom.add_line(ptsExtra1[i+1],ptsExtra0[i+1])
                lloop = geom.add_line_loop([lf,l1,lr,l0])
                extra = geom.add_ruled_surface(lloop)
                extrados.append(extra)
                geom._GMSH_CODE.append('Recombine Surface {%s};' % extra)
                geom.add_physical_surface(extra)
                l0 = geom.add_line(ptsIntra0[i],ptsIntra0[i+1])
                l1 = geom.add_line(ptsIntra1[i+1],ptsIntra1[i])
                lf = geom.add_line(ptsIntra1[i],ptsIntra0[i])
                lr = geom.add_line(ptsIntra0[i+1],ptsIntra1[i+1])
                lloop = geom.add_line_loop([l0,lr,l1,lf])
                intra = geom.add_ruled_surface(lloop)
                geom._GMSH_CODE.append('Recombine Surface {%s};' % intra)
                geom.add_physical_surface(intra)
                intrados.append(intra)
            l0 = geom.add_line(ptsIntra1[-1],ptsIntra0[-1])
            l1 = geom.add_line(ptsExtra1[-1],ptsIntra1[-1])
            geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l1,3) )
            l2 = geom.add_line(ptsExtra0[-1],ptsExtra1[-1])
            l3 = geom.add_line(ptsIntra0[-1],ptsExtra0[-1])
            geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (l3,3) )
            lloop = geom.add_line_loop([l0,l1,l2,l3])
            tr_ed = geom.add_ruled_surface(lloop)
            geom._GMSH_CODE.append('Recombine Surface {%s};' % tr_ed) 
            geom.add_physical_surface(tr_ed)
            if ispan == 0:
                linesEx = []
                linesIn = []
                for i in range(len(ptsExtra0)-1):
                    l0 = geom.add_line(ptsExtra0[i],ptsExtra0[i+1])
                    linesEx.append(l0)
                lte = geom.add_line(ptsExtra0[-1],ptsIntra0[-1])
                geom._GMSH_CODE.append('Transfinite Line{%s} = %s Using Progression 1;' % (lte,3) )
                for i in range(len(ptsExtra0)-1):
                    l1 = geom.add_line(ptsIntra0[::-1][i],ptsIntra0[::-1][i+1])
                    linesIn.append(l1)
                #lloop = geom.add_line_loop(zob)
                #base = geom.add_plane_surface(lines)
                #geom.add_physical_surface(base)
            #
        fg = open(name+'.geo','w')
        for l in geom.get_code():
            fg.write(l)

        self.points, self.cells = pg.generate_mesh(geom)

        self.X = np.array([p[0] for p in self.points])
        self.Y = np.array([p[1] for p in self.points])
        self.Z = np.array([p[2] for p in self.points])
        self.quads = self.cells['quad']
        meshio.write(name+'.vtu', self.points, self.cells)# ells['quad'])    #{'mu':mu}

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
            #p2 = vlm.fluidPoint(p[0],p[1],p[2])
            p2.computeVelocity(elem)
            fp.append(p2)

        self.d={'ux' : np.array([p.vel[0] for p in fp]),\
                'uy' : np.array([p.vel[1] for p in fp]),\
                'uz' : np.array([p.vel[2] for p in fp])}
        meshio.write('volumic1.vtu', self.points, self.cells, point_data=self.d)
