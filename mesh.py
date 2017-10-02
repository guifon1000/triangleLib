import numpy as np
import meshio
import pygmsh as pg
import matplotlib.pyplot as plt
import triangleLib as tl
import sys
import bidLib as l2d
import json
from frames import Frame
"""
    mesh lib 
"""





def read_msh_triangulation(fi, **kwargs):
    f = open(fi,'r').readlines()
    _vertices = []
    app = False
    ist = 0
    nv = None
    while app == False :
        l = f[ist] 
        if l.startswith('$Nodes'):
            nv = int(f[ist+1])
            app = True
        ist+=1
    for i in range(ist+1,ist+nv+1):
        l=f[i].split()
        _vertices.append((float(l[1]), float(l[2]), float(l[3])))
    _faces = []
    for i in range(ist+nv+1,len(f)):
        lp = f[i].split()
        if (len(lp) > 3) and (lp[1] == '2') :
                _faces.append((int(lp[5]),\
                                  int(lp[6]),\
                                  int(lp[7])))
    d = Triangulation(_vertices, _faces)
    return d




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

def circle(rad,xc,n) : #xc is a 2 dim tuple, position of the center
    first_polyline = []
    x0 = xc[0]
    y0 = xc[1]
    for i in range(n):
        teta = i*2.*np.pi/n
        first_polyline.append(l2d.point(  (x0 + rad * np.cos(teta) , y0 + rad * np.sin(teta) )  ,  z=0.))
    return first_polyline
        
def circle_in_box(name):
    out = Polyline2D(rectangle(25.,5.,(5.,0.)))
    objects = [Polyline2D(circle(0.1,(0.,0.),200))]
    geom = box(out,objects,name)
    return geom

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
    geom = box2D(out,objects,name)
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
    geom = box2D(out,objects,name)
    #poly2d.box_2d(name = 'vol0')
    #geom.add_line(pt0,pt1)
    return geom

def write_geo(name,geom):
    fg = open(name+'.geo','w')
    for l in geom.get_code():
        fg.write(l)
    fg.close()






def box3D(Lx,Ly,Lz,inObjects,**kwargs):
    print('----------------- 3D box construction ------------------------')
     




def box2D(outBox,inObjects,name):
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
        if kwargs.has_key('closed'):
            self.check_closed = kwargs['closed']
        else : 
            self.check_closed = False
        if self.check_closed : 
            p = self[0]
            self.append(p)
        self.pt3d = []
        for i in range(len(self)):
            p = l2d.point(self[i],z=self.z)
            self.pt3d.append( tl.Point([p[0] , p[1], self.z]))
    
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



    def add_to_geom(self, geom):
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

    def translate(self,vec):
        for i in range(len(self)) :
            self.pt3d[i] = [self.pt3d[i][j]+vec[j] for j in range(3)]

    def center_of_gravity(self):
        s = np.array(self.pt3d)
        self.cg = np.mean(s, axis = 0)
            
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


class Triangulation(dict):
    """
    dictionary
    * KEYS : 
        'vertices' : array of points, 3-real arrays (x,y,z)
        'faces'    : array of faces, 3 int arrays (3 index in vertices list)
    """
    def __init__(self, points, faces, **kwargs):
        self['vertices'] = points
        self['faces'] = faces
   

    def write_obj_file(self,name):
        with open(name+'.obj','w') as f:
            f.write('# OBJ file\n')
            for v in d['vertices']:
                f.write('v %.4f %.4f %.4f\n' % (v[0],v[1],v[2]))
            for t in d['faces']:
                f.write('f')
                for p in t :
                    f.write(" %d" % p)
                f.write('\n')




 
    def reorient_convex(self):
        # reorient the faces of the icosahedron
        _new_faces = []
        for f in self['faces']:
            p0 = tl.Point(self['vertices'][f[0]-1])
            p1 = tl.Point(self['vertices'][f[1]-1])
            p2 = tl.Point(self['vertices'][f[2]-1])
            tr = tl.Triangle(p0, p1, p2)
            cg = tr.gravityCenter()
            n0 = tr.computeNormal()
            n1 = tl.Vector(cg)
            if tl.dot(n0,n1) < 0:
                fi = (f[0], f[2], f[1])
            else:
                fi = f
            _new_faces.append(fi)
   
        self['faces'] = _new_faces 
        # ok the icosahedron is correctly oriented in the dict d


    def refine_1(self) :
        d2 = {}
        _vertices = []
        _faces = []
        for f in self['faces']:
            p0 = self['vertices'][f[0]-1]
            p1 = self['vertices'][f[1]-1]
            p2 = self['vertices'][f[2]-1]
            i0 = None
            i1 = None
            i2 = None
            i3 = None
            i4 = None
            i5 = None
            if p0 not in _vertices:
                _vertices.append(p0) 
                i0 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p0:
                        i0 = i
                        break
            if p1 not in _vertices:
                _vertices.append(p1) 
                i1 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p1:
                        i1 = i
                        break
            if p2 not in _vertices:
                _vertices.append(p2)
                i2 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p2:
                        i2 = i
                        break
            p3 = [0.5*(p1[j]+p2[j]) for j in range(3)]  
            p4 = [0.5*(p2[j]+p0[j]) for j in range(3)]  
            p5 = [0.5*(p1[j]+p0[j]) for j in range(3)] 
            if p3 not in _vertices:
                _vertices.append(p3) 
                i3 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p3:
                        i3 = i
                        break
            if p4 not in _vertices:
                _vertices.append(p4) 
                i4 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p4:
                        i4 = i
                        break
            if p5 not in _vertices:
                _vertices.append(p5)
                i5 = len(_vertices)-1
            else:
                for i,v in enumerate(_vertices):
                    if v==p5:
                        i5 = i
                        break
            i0 += 1 
            i1 += 1 
            i2 += 1 
            i3 += 1 
            i4 += 1 
            i5 += 1
            _faces.append((i4, i3, i2))
            _faces.append((i5, i1, i3))
            _faces.append((i5, i3, i4))
            _faces.append((i0, i5, i4))
        _vertices2 = []
        for p in _vertices:
            d = np.sqrt(p[0]**2. + p[1]**2. + p[2]**2.)
            ps = (p[0]/d ,\
                p[1]/d ,\
                p[2]/d)
            _vertices2.append(ps) 
        self['vertices'] = _vertices2
        self['faces'] = _faces


def thick_wing(pf,name):
    geom = pg.Geometry()
    pts = []
    lns = []
    pol = pf.polyline().pt3d
    for (i,p) in enumerate(pol[:-1]):
        p = geom.add_point(p,0.1)
        pts.append(p)
    for i in range(len(pts)-1):
        l = geom.add_line(pts[i],pts[i+1])
        lns.append(l)
    l = geom.add_line(pts[-1],pts[0])
    lns.append(l)
    lloop = geom.add_line_loop(lns)
    sfStart = geom.add_plane_surface(lloop)
    geom.add_physical_surface(sfStart) 
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






def wing():

    # creation of a volumic wing with the extrusion of a 2d profile on a 
    # generatrix defined by a 3d-spline

    sys.path.append('./profiles')
    from splineProfileMultiParam import  Profile
    from scipy import interpolate
    from mpl_toolkits.mplot3d import Axes3D


    # creation of the geometry pygmsh
    geom = pg.Geometry()
    
    # control points of the generatrix
    x = [0.0, 0.1, 0.4, 0.6, 0.9]
    y = [0.0, 0.05, 0.09, 0.11, 0.4]
    z = [0.0, 0.3, 0.8, 1.2, 2.9]

    # tck, u represent the parametric 3d curve
    tck, u = interpolate.splprep([x,y,z], s=2)
    name = 'wingZero'
    Nslices = 100 # number of slices
    npt = 63 # points of the profile
    t = np.linspace(0., 1., Nslices) # parametric space
    pf = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = npt) # creation of the 2d profile
    fi = Frame(t[0], tck, type = 'Xnat')
    pol = pf.polyline()
    pol.to_frame(fi, scale = 0.5)
    li = pol.add_to_geom(geom)
    lloop = []
    for l in li : lloop.append(l)
    ll = geom.add_line_loop(lloop)
    li0 = li
    sf = geom.add_plane_surface(ll)
    phys = []
    geom.add_physical_surface(sf)
    phys.append(sf)
    for i in range(Nslices-1):
        si = t[i]
        sip1 = t[i+1]
        #pf = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = npt)
        fip1 = Frame(sip1, tck, type = 'Xnat')
        pol = pf.polyline()
        pol.to_frame(fip1, scale = 0.5*np.cos(sip1*0.4*np.pi))
        lip1 = pol.add_to_geom(geom)
        for j in range(len(li)):
            lij = li0[j]
            lip1j = lip1[j]
            lti = geom.add_line(lij.points[0], lip1j.points[0])
            ltip1 = geom.add_line(lij.points[1], lip1j.points[1])
            lloop = geom.add_line_loop([lti, lip1j, -ltip1, -lij])
            sf = geom.add_ruled_surface(lloop)
            phys.append(sf)
        li0 = lip1
    lloop = []
    for l in lip1 : lloop.append(-l)
    ll = geom.add_line_loop(lloop)
    li0 = li
    sf = geom.add_plane_surface(ll)
    phys.append(sf)
    physS = geom.add_surface_loop(phys)
    geom.add_physical_surface(phys, label = 'profile')
    write_geo(name, geom)
    import subprocess
    exe_gmsh = '/home/fon/gmsh-2.16.0-Linux/bin/gmsh'
    subprocess.call([exe_gmsh,name+'.geo','-2','-o',name+'.msh','>',name+'.mshlog'])
    X, cells, pt_data, cell_data, field_data = meshio.read(name+'.msh')
    box3D(10.,10., 10., [0])
    print X
    print cells['triangle']
    import readMSH as rmsh
    d = rmsh.read_msh_file(name)
    rmsh.write_fms_file(name,**d)
    d = read_msh_triangulation(name+'.msh')
    d.write_obj_file(name)


if __name__ == '__main__1':
    from random import random
    Npoints = 5
    scale = 10.
    name = 'zob'
    geom = pg.Geometry()
    points = []
    lines = []
    thicknesses = [0.1*scale*np.cos(float(i)) for i in range(Npoints) ]

    p0 = tl.Point((0., 0., 0. ))
    p = geom.add_point(p0, scale)
    points.append(p)

    for i in range(Npoints):
        p0 = (random() - 0.5, random() - 0.5, 0.)
        p = geom.add_point(p0, scale)
        l = geom.add_line(points[0], p)
        points.append(p)
        lines.append(l)
    sector = []
    A = points[0].x
    B = lines[0].points[1].x
    first = tl.Vector((B[i]-A[i] for i in range(3)))
    for l in lines[1:]:
        B = l.points[1].x

        vec = tl.Vector((B[i]-A[i] for i in range(3)))
        print '------------------------------'
        print ' first = ' + str(first)
        print ' vec = '+ str(vec)
        angle = tl.angle(first,vec)
        print ' angle = '+ str(angle)
        sector.append((l.id, angle))
    def get_angle(tup):
        return tup[1]
    sort_sect = sorted(sector, key = get_angle)
    print sort_sect
    print ' - - - - - - - - - - - - '
    sectors = []
    sectors.append((lines[0].id, sort_sect[0][0])) 
    for i in range(len(sort_sect)-1):
        sectors.append((sort_sect[i][0], sort_sect[i+1][0]))
    sectors.append((sort_sect[-1][0], lines[0].id)) 
    print sectors
    # at this point, sectors contains the ordered couples of lines 
    thick = 0.1
    
    for sec in sectors :
        # we look for a point at distance = thickness of each line
        print sec
        l0 = sec[0]
        l1 = sec[1]
        for l in lines :
            if l.id == sec[0] : l0 =l
            if l.id == sec[1] : l1 =l
        print l0.points
    #    print 'line 1 : '+str(lines[sec[0]].points[:].x)
    #    print 'line 2 : '+str(lines[sec[1]].points[:].x)
    write_geo(name, geom)
    d = geom.__dict__
    print d






if __name__ == '__main__3':
    # find what is the closest of an icosahedron
    d= read_msh_triangulation('./icosahedron.msh')
    d.reorient_convex()
    for i in range(1):
        d.refine_1() 
    d.write_obj_file('icosahedron')


