import sys
sys.path.append('./classes/')
from classes.Point import Point
from classes.Vector import Vector
from classes.Frame import Frame
from scipy import interpolate
import numpy as np

def write_geo(name,geom):
    fg = open(name+'.geo','w')
    for l in geom.get_code():
        fg.write(l)
    fg.close()


def cross(u,v):
    pv=[]
    pv.append(u[1]*v[2]-u[2]*v[1])
    pv.append(u[2]*v[0]-u[0]*v[2])
    pv.append(u[0]*v[1]-u[1]*v[0])
    vp=Vector(pv)
    return vp


def distance(p1, p2):
    return Vector([p2[i] - p1[i] for i in range(3) ]).norm
    


def dot(u,v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]


def angle(u,v, plane_normal = None):
    x = u.unit()
    y = v.unit()
    angle = np.arctan2(cross(x,y).norm,dot(x,y))
    #
    if plane_normal:
        if dot(cross(u,v), plane_normal)>0.:
            return angle
        else :
            return -angle
    else:
        return angle
  
    # warning : only in 2D case
    #if dot(cross(u,v), Vector((0.,0.,1.0)))>0.:
    #    return angle
    #else :
    #    print "gne"
    #    return -angle


def intersect_3_planes(p1,p2,p3):
    A=np.array([p1[:3], p2[:3], p3[:3]])
    b=np.array([-p1[3],-p2[3],-p3[3]]).transpose()
    sol=np.linalg.solve(A, b)
    return Point([sol[i] for i in range(3)])


def get_parameter(pt, line):
    if is_on_line(pt, line):
        imax = np.argmax([np.abs(c) for c in line[1]])
        return (pt[imax]-line[0][imax])/line[1][imax]
    else:
        return None




def intersect_2_lines(l1, l2):
    from classes.Plane import Plane
    v1 = l1[1]
    v2 = l2[1]
    normal_vect = cross(l1[1], l2[1])
    intersection = None
    if normal_vect.norm > 1.e-16:
        unit_normal = normal_vect.unit()
        components = [c**2. for c in unit_normal]
        # np.argpartition(components, -2) gives the array of indices of growing component 
        # (from smaller to biggest)
        # ind = np.argpartition(components, -2)[:2] gives the indices of the 2 smallest 
        # components of the vector
        ind = np.argpartition(components, -2)[:2]

        plane_1 = Plane((l1[0], unit_normal))

        if is_on_plane(l2[0], plane_1):
            i1 = ind[0]
            i2 = ind[1]
            A = np.array([[-v1[i1], v2[i1]], [-v1[i2], v2[i2]]])
            b = np.array([l1[0][i1] - l2[0][i1], l1[0][i2] - l2[0][i2]]).transpose()
            #print A
            #print b
            sol=np.linalg.solve(A, b)
            #print sol
            intersection = l1.parameter_point(sol[0]) 
    return intersection


def intersect_2_segments(seg1, seg2):
    from classes.Line import Line
    v1 = Vector([seg1[1][i] - seg1[0][i] for i in range(3)])
    v2 = Vector([seg2[1][i] - seg2[0][i] for i in range(3)])
    l1 = Line([seg1[0], v1])
    seg1_parameters = sorted([get_parameter(seg1[0], l1), get_parameter(seg1[1], l1)])
    l2 = Line([seg2[0], v2])
    seg2_parameters = sorted([get_parameter(seg2[0], l2), get_parameter(seg2[1], l2)])
    intersection = intersect_2_lines(l1,l2)
    out = False
    if intersection and ( seg1_parameters[0] < get_parameter(intersection, l1) < seg1_parameters[1]) and \
            ( seg2_parameters[0] < get_parameter(intersection, l2) < seg2_parameters[1]):
        out = True
    return out



def is_on_line(pt, line):
    c_ref = None
    s= None
    for i,coord in enumerate(line[1]):
        if coord != 0.:
            s = (pt[i] - line[0][i])/coord
            break
    pt_test = Point([line[0][i] + s * line[1][i] for i in range(3)])
    return distance(pt_test, pt)<1.e-10

def is_on_plane(pt, plane):
    return (plane[0]*pt[0] + plane[1]*pt[1] + plane[2]*pt[2] + plane[3])**2. < 1.e-20

def matrix_to_quaternion(m):
    from classes.Quaternion import Quaternion
    diag=np.diag(m)
    np.append(diag,1.)

    tr= np.trace(m)+1.
    if tr>0.:
        s=0.5/np.sqrt(tr)
        x=(m[2,1]-m[1,2])*s
        y=(m[0,2]-m[2,0])*s
        z=(m[1,0]-m[0,1])*s
        w=0.25/s
        return Quaternion((w,x,y,z))
         
def vtk_visu(*largs, **kwargs):
    import vtk
    print 'VTK VISU'
    points = []
    polylines = []
    frames = []
    idpoint = 0
    nbpoint = 0 
    for v in largs[0] : 
	if type(v).__name__ == 'Polyline2D':
           first = len(points)
           nbpoint = len(v)
           pol = []
           for p in v.pt3d:
               points.append(p)
	       pol.append(len(points) - 1)
           pol.append(first)
           polylines.append(pol) 
	if type(v).__name__ == 'Frame':
           first = len(points)
           nbpoint = len(v)
           pol = []
           points.append(v[0])
           i_orig_frame = len(points) - 1
           points.append(       [v[0][i] + v[1][0][i] for i in range(3)] )
           i_x = len(points) - 1
           points.append(       [v[0][i] + v[1][1][i] for i in range(3)] )
           i_y = len(points) - 1
           points.append(       [v[0][i] + v[1][2][i] for i in range(3)] )
           i_z = len(points) - 1
	   polx = [ i_orig_frame, i_x ]
           polylines.append(polx)
	   poly = [ i_orig_frame, i_y ]
           polylines.append(poly) 
	   polz = [ i_orig_frame, i_z ]
           polylines.append(polz)  
    pts = vtk.vtkPoints() 
    lns = vtk.vtkCellArray()
    pts.SetNumberOfPoints(len(points))
    for i,p in enumerate(points):
        pts.SetPoint(i, p[0], p[1], p[2])

    for pol in polylines:
        lns.InsertNextCell(len(pol)) 
        for p in pol:
            lns.InsertCellPoint(p)
        

    polygon = vtk.vtkPolyData()
    polygon.SetPoints(pts)
    polygon.SetLines(lns)



    polygonMapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        polygonMapper.SetInputConnection(polygon.GetProducerPort())
    else:
        polygonMapper.SetInputData(polygon)
        polygonMapper.Update()
    
    polygonActor = vtk.vtkActor()
    polygonActor.SetMapper(polygonMapper)        


    ren1 = vtk.vtkRenderer()
    ren1.AddActor(polygonActor)
    ren1.SetBackground(0.1, 0.2, 0.4)

    ren1.ResetCamera()



    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren1)
    renWin.SetSize(300, 300)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    iren.Initialize()
    iren.Start()


def parameter_frame(tck, s, mode = 'frenet'):
    from classes.Frame import Frame
    basis = []
    t =  interpolate.splev(s , tck, der = 0)
    orig = Point([t[0], t[1], t[2]])
    if mode == 'frenet':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        t =  interpolate.splev(s , tck, der = 2)
        xpp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T.unit())
        B = cross(xp, xpp)
        basis.append(B.unit())
        N = cross(B, T)
        basis.append(N.unit())

    if mode == 'Xnat':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T)
        B = cross((1.,0.,0.), T)
        basis.append(B)
        N = cross(T, B)
        basis.append(N)

    if mode == 'Ynat':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T)
        B = cross((0.,1.,0.), T)
        basis.append(B)
        N = cross(T, B)
        basis.append(N)


    if mode == 'Znat':
        t =  interpolate.splev(s , tck, der = 1)
        xp  = Vector( [float(t[0]), float(t[1]), float(t[2]) ])
        T = xp.unit()
        basis.append(T)
        B = cross((0.,0.,1.), T)
        basis.append(B)
        N = cross(T, B)
        basis.append(N)

    matrix = np.zeros((3,3), dtype =float)
    for i in range(3) : matrix[i] = basis[i]
    return Frame((orig, matrix))




def write_fms_file(name,kwargs):
    f = open(name+('.fms'),'w')
    d=kwargs
    print '------------------------------------------'
    print '------------------ FMS PART --------------'
    print '------------------------------------------'
    print d['faces']
    f.write('// patch names and types\n')
    f.write(str(len(d['faces']))+'\n')
    f.write('(\n')
    for p in d['faces']:
        f.write(str(p)+'\n')
        if 'cylinder' in p:
            f.write('wall\n')
        else:
            f.write('patch\n')
        f.write('\n')
    f.write(')\n\n')
    f.write('// coordinates of surface points\n')
    f.write(str(len(d['vertices']))+'\n')
    f.write('(\n')
    for p in d['vertices']:
        f.write( '( %f %f %f) '%(p[0], p[1], p[2] )   )
    f.write(')\n\n')
    f.write('// list of triangles\n')
    ntri = 0
    for k in d['faces']:
        ntri += len(d['faces'][k])
    f.write(str(ntri)+'\n')
    f.write('(\n')
    for i,group in enumerate(d['faces']):
        for j,t in enumerate(d['faces'][group]):
            tp = '( ( %d %d %d )  %d ) '% (t[0] - 1 , t[1] - 1 , t[2] - 1 , i   )
            f.write(tp)
    f.write(')\n\n')

    for i in range(4):f.write('0()\n')                                                                     


