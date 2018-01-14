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


def dot(u,v):
    return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]



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
    matrix = np.zeros((3,3), dtype =float)
    for i in range(3) : matrix[i] = basis[i]
    return Frame((orig, matrix))





