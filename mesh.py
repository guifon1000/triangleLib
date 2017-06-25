import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
import bidLib as l2d
import json

"""
    mesh lib 
"""







def square(L,x=0.):
    first_polyline = []
    first_polyline.append(l2d.point((-L,-L),z=0.))
    first_polyline.append(l2d.point((L,-L),z=0.))
    first_polyline.append(l2d.point((L,L),z=0.))
    first_polyline.append(l2d.point((-L,L),z=0.))
    pol2 = first_polyline 
    for p in pol2:
        p[0]+=x
    return pol2

#def read_polyline_file(name):
    #d = json.load(open(name+'.json','r'))
    #return d

#def create_closed_object(**kwargs):


def soluz_box(L,H,D,ld):
    first_polyline = []
    first_polyline.append(l2d.point((0.,0.5*float(H)+D),z=0.))
    first_polyline.append(l2d.point((0.5 * float(L)+ ld , 0.5*float(H)+D ),z=0.))
    first_polyline.append(l2d.point((0.5 * float(L)+ ld , 0.5*float(H)),z=0.))
    first_polyline.append(l2d.point((0.5 * float(L) , 0.5*float(H)),z=0.))
    first_polyline.append(l2d.point((0.5 * float(L) , 0.),z=0.))
    pol = first_polyline
    for p in first_polyline[::-1][1:]:
        pol.append(l2d.point((p[0],-p[1]),z=0.))
    pol2 = pol
    for p in pol[::-1][1:len(pol)-1]:
        pol2.append(l2d.point((-p[0],p[1]),z=0.))
    return pol2

def rectangle(lx,ly,xc) : #xc is a 2 dim tuple, position of the center
    first_polyline = []
    x0 = xc[0]
    y0 = xc[1]
    first_polyline.append(l2d.point(  (x0 + 0.5*lx , y0 - 0.5 * ly)  ,  z=0.))
    first_polyline.append(l2d.point(  (x0 - 0.5*lx , y0 - 0.5 * ly)  ,  z=0.))
    first_polyline.append(l2d.point(  (x0 - 0.5*lx , y0 + 0.5 * ly)  ,  z=0.))
    first_polyline.append(l2d.point(  (x0 + 0.5*lx , y0 + 0.5 * ly)  ,  z=0.))
    return first_polyline






def motif0(L,l,e,a,thickness,x=0.):
    first_polyline = []
    first_polyline.append(l2d.point((0.,0.5*float(L)+thickness),z=0.))
    first_polyline.append(l2d.point((0.5 * float(l)+ 0.5* (e-a),0.5*float(L)+thickness),z=0.))
    first_polyline.append(l2d.point((0.5 * float(l)+ 0.5* (e-a),0.5*float(L)),z=0.))
    first_polyline.append(l2d.point((0.5 * float(l) , 0.5*float(L)),z=0.))
    first_polyline.append(l2d.point((0.5 * float(l) , 0.),z=0.))

    pol = first_polyline
    for p in first_polyline[::-1][1:]:
        pol.append(l2d.point((p[0],-p[1]),z=0.))
        
    pol2 = pol
    for p in pol[::-1][1:len(pol)-1]:
        pol2.append(l2d.point((-p[0],p[1]),z=0.))
    for p in pol2:
        p[0]+=x
    return pol2

def revolve(pol):
    geom = pg.Geometry()
    pts = []
    lns = []
    lcar = 0.01
    for i,p in enumerate(pol.pt3d):
        pts.append(geom.add_point(p,lcar))
    for i,p in enumerate(pts[:-1]):
        lns.append(geom.add_line(p,pts[i+1]))
        
    write_geo('revolve',geom)
    


def simple_square(name):
    objects = []
    pol2 = rectangle( 1., 1. , (0. , 0.)   )
    out = Polyline2D(pol2)
    geom = box(out,objects,name)
    return geom



def echanger(name):
    thickness = 0.002
    Larg = 3.
    L = 1.
    l = 0.0025
    e = 0.005
    a = 0.001
    
    thickness = 0.002
    H = 1.
    D = 0.025
    L = 2.
    ld = D

    No = 200
    a = (L+(1-No)*thickness)/No
    objects = []
    for i in range(No-1):
        #pol2 = motif0(L,l,e,a,thickness,x=-0.5*Larg+Larg*float(i)*(1./float(No-1)))
        pol2 = rectangle(thickness, H , (-0.5*L+a+0.5*thickness + float(i) * (a+thickness)   , 0.)  )
        objects.append(Polyline2D(pol2))
    #pol2 = motif0(1.05*L,1.05*Larg,0.05,0.00025,0.05)
    pol2 = soluz_box(L,H,D,ld)
    out = Polyline2D(pol2)
    geom = box(out,objects,name)
    #poly2d.box_2d(name = 'vol0')
    #geom.add_line(pt0,pt1)
    return geom


def square_in_box(name):
    gmp = []
    L = 0.5
    l = 0.1
    objects = []
    No = 1
    Larg = 3.
    for i in range(No):
        pol2 = square(l)
        objects.append(Polyline2D(pol2))

    pol2 = square(L)
    out = Polyline2D(pol2)
    geom = box(out,objects,name)
    #poly2d.box_2d(name = 'vol0')
    #geom.add_line(pt0,pt1)
    return geom

def write_geo(name,geom):
    fg = open(name+'.geo','w')
    for l in geom.get_code():
        fg.write(l)
    fg.close()


def box(outBox,inObjects,name):
    print str(len(inObjects))+' polylines !'
    geom = pg.Geometry()
    boxPoints = []
    boxLines = []
    boxPoints2 = []
    boxLines2 = []
    lcar = 100.
    for pto in outBox.pt3d:
        boxPoints.append(geom.add_point(pto,lcar))
        p2 = (pto[0],pto[1],pto[2]+0.1)
        boxPoints2.append(geom.add_point(p2,lcar))
    for i in range(len(boxPoints)-1):
        boxLines.append(geom.add_line(boxPoints[i],boxPoints[i+1]))
        boxLines2.append(geom.add_line(boxPoints2[i],boxPoints2[i+1]))
    boxLines.append(geom.add_line(boxPoints[-1],boxPoints[0]))
    boxLines2.append(geom.add_line(boxPoints2[-1],boxPoints2[0]))

    boxloop = []
    boxloop2 = []
    objloop = []
    objloop2 = []
    bcBox = {}
    bcObj = {}
    for il,l in enumerate(boxLines):
        boxloop.append(l)
        ext = geom.extrude(l,translation_axis=(0.,0.,0.1))
        geom._GMSH_CODE.append('Reverse Surface{%s};' % ext[1].id )
        geom.add_physical_surface(ext[1],label='BC.'+str(il))
        #bcBox['BC.'+str(il)] = outBox.types[il]

    for il,l in enumerate(boxLines2):
        boxloop2.append(l)
    for io,ob in enumerate(inObjects):
        opl = []
        opl2 = []
        lobj = []

        for p in ob.pt3d:
            pp = geom.add_point(p,1.)
            opl.append(pp)
        for i,p in enumerate(opl[:-1]):
            li = geom.add_line(opl[i],opl[i+1])
            lobj.append(li)
            objloop.append(li)
        #close the polyline
        li = geom.add_line(opl[-1],opl[0])
        lobj.append(li)
        objloop.append(li)

        for p in ob.pt3d:
            pp = geom.add_point((p[0],p[1],p[2]+0.1),1.)
            opl2.append(pp)
        for i,p in enumerate(opl2[:-1]):
            li = geom.add_line(opl2[i],opl2[i+1])
            objloop2.append(li)
        #close the polyline
        li = geom.add_line(opl2[-1],opl2[0])
        objloop2.append(li)

    lengthes = {}
    for l in objloop:
        le0 = np.sqrt((l.points[0].x[0]-l.points[1].x[0])**2.+(l.points[0].x[1]-l.points[1].x[1])**2.)
        if (str(le0) not in lengthes.keys()):
            lengthes[str(le0)] = []
            lengthes[str(le0)].append(l)
        else:
            lengthes[str(le0)].append(l)

    idl = 0
    for kl in lengthes.keys():
        sfOb = []
        idl += 1
        lab = 'wall.'+str(idl)
        idli = 1
        for li in lengthes[kl]:
            labi = lab
            labi += '.'+str(idli)
            ext = geom.extrude(li,translation_axis=(0.,0.,0.1))
            sfOb.append(ext[1])
        geom.add_physical_surface(sfOb,label=lab)
            #idli += 1 
    #all_lines = objloop+boxloop
    #all_lines2 = objloop2+boxloop2
    #lzmin = geom.add_line_loop(all_lines)
    #lzmax = geom.add_line_loop(all_lines2)
    #zmax = geom.add_plane_surface(lzmax)
    #zmax = geom.add_plane_surface(lzmax)
    #geom._GMSH_CODE.append('Reverse Surface{%s};' % zmax.id )
    #Zmin = geom.add_physical_surface(zmin,label='Zmin')
    #Zmax = geom.add_physical_surface(zmax,label='Zmax')
    write_geo(name,geom)
    return geom


class Polyline2D(list):    #always closed
    def __init__(self, *largs,**kwargs):
        super(Polyline2D,self).__init__(*largs)
        if kwargs.has_key('z'):
            self.z=kwargs['z']
        else:
            self.z = 0.
        if kwargs.has_key('types'):
            self.types = kwargs['types']



        self.pt3d = []
        for i in range(len(self)):
            p = l2d.point(self[i],z=self.z)
            self.pt3d.append( (p[0] , p[1], self.z))
        


    def write_xfoil(self, **kwargs):
        print('=============  write XFOIL file==============')
        name = kwargs['name']
        import subprocess
        subprocess.call(['mkdir',name])
        f = open(name+'/'+name+'.dat','w')
        f.write('\n')
        for doublon in self:
            f.write(str(doublon[0])+' '+str(doublon[1])+'\n')
        f.close()
        return name+'/'+name+'.dat'
            
            
            
    def extrude_segments(self,ep):
        pass

    def box_2d(self,**kwargs):
        """
          3  <------------- 2
          |                 ^
          v                 |
          0  -------------> 1

          creates a geometry of the boxed geometry, with a thickness
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
        # 4 corner box
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
        # N corner box
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
        lloop2 = geom.add_line_loop(ll2)
        sf1 = geom.add_plane_surface(lloop2)
        geom._GMSH_CODE.append('Reverse Surface{%s};' % sf1.id )

        write_geo(self.name,geom)
        print('2D Geometry successfully created')
        points,cells,pt_data, cell_data, field_data = pg.generate_mesh(geom)
        meshio.write('volumic0.vtu', points, cells)# ells['quad'])    #{'mu':mu}




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



if __name__=='__main__':
    points = [(0.,1.),(1.,1.),(0.,0.)]
    pol = Polyline2D(points)
    revolve(pol)
