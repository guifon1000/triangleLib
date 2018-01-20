import sys
sys.path.append('../')
from classes.Point import Point
from classes.Triangle import Triangle
from classes.Triangulation import Triangulation, read_msh_file, merge
from classes.Vector import Vector
from classes.Sphere import Sphere
from classes.Frame import Frame
from classes.Quaternion import Quaternion
from classes.Extrusion import Extrusion
from classes.Box import Box
import pygmsh as pg
from modelers.profiles.splineProfileMultiParam import  Profile
from modelers.planet.Planet import  Planet
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from functions import parameter_frame, matrix_to_quaternion, vtk_visu, write_geo, write_fms_file


def pretty_print(string, symbol = '-'):
    s = ''
    for i in range(79):
        s += symbol
    s += '\n'
    lstr = len(string) + 2 
    s2 = ''
    if lstr >= 79:
        s2 = string + '\n'
    else:
        si = ' ' + string + ' '
        n1 = int((79 - lstr )/2)
        n2 = 79 - len(si) - n1
        s2 = ''
        for i in range(n1): s2 += symbol
        s2 += si
        for i in range(n2): s2 += symbol 
        s2 += '\n'
    print( s + s2 + s)







idtest = 1
pretty_print('TEST n.'+str(idtest) + ': TRIANGLE')


# create three points
p1 = Point((0.001,0.004,0.))
p2 = Point((0.23,1.,0.))
p3 = Point((1.,0.4,0.7))

points = [p1, p2, p3]

geo_points = []
# create a triangle
t = Triangle(points)

print ' ----- A TRIANGLE ------ '
print t
print 'attributes'
print 'cog = '+str(t.cg)
print 'normal = '+str(t.normal)

print ' REVERSE ---------> '
t = -t
print t
print 'attributes'
print 'cog = '+str(t.cg)
#t.cg.pop_to_geom(geom)
print 'normal = '+str(t.normal)

idtest += 1
pretty_print('TEST n.' + str(idtest) + ': VECTOR')
v = t.normal
print "raw vector :"+str(v)
print "norm 2 : "+str(v.norm2)
print "norm 1 : "+str(v.norm)
print "unit vector : "+str(v.unit())



idtest += 1
pretty_print('TEST n.' + str(idtest) + ': IMPORT JSON TRIANGULATION')


d = Triangulation()
d.load_file('../samples/icosahedron.json')
print d.cg


idtest += 1
pretty_print('TEST n.' + str(idtest) + ': SPHERES')


d.reorient_convex()
d.translate((1.,0.,-2.))
d = d.refine_2()
d = d.refine_2()
d = d.refine_2()
d.write_obj_file('refined_icosahedron.obj')
s1 = Sphere(refin = 0, center = (1.,1.,1.), radius = 0.3)
s2 = Sphere(refin = 1, center = (0.,-4.,2.), radius = 0.1)
s3 = Sphere(refin = 2, center = (0.,0.,0.), radius = 0.7)
p1 = Planet(refin = 1, center = (-3.,-5.,-3.), radius = 0.8)
s1.write_obj_file('sphere1.obj')
write_fms_file('sphere', s1)
s2.write_obj_file('sphere2.obj')
s3.write_obj_file('sphere3.obj')
p1.write_obj_file('planet1.obj')

idtest += 1
pretty_print('TEST n.' + str(idtest) + ': EXTRUSION ALONG 3D CURVE')

geom = pg.built_in.Geometry()
vtk_elements = []

pf = Profile(typ = 'fon',par = [0.99,0.1,0.018,0.035,0.001],npt = 21) # creation of the 2d profile
pol = pf.polyline(closed = True)  # TODO : Profile should herit of Polyline_2D

print type(pf)
# control points of the generatrix
global_scale = 0.09

x = [-6., -5., 0., 5., 6.]
y = [0., -0.5, -1.5, -0.5, 0.]
z = [0.0, 0.,0., 0., -0.00]


x = [v * global_scale for v in x]
y = [v * global_scale for v in y]
z = [v * global_scale for v in z]


# tck, u represent the parametric 3d curve
tck, u = interpolate.splprep([x,y,z], s=3)

frames = []
scales = []

for s in np.linspace(0., 1., num = 20):
    frame = parameter_frame(tck, s, mode = 'Znat')
    scale =  global_scale * (0.5 + 5.* s * (1. - s))
    scales.append(scale)
    vtk_elements.append(frame)
    frames.append(frame)

extrusion = Extrusion(pol, frames, scales = scales, close_caps = True, mode = 'profile')
print extrusion['physical']
write_fms_file('extrusion', extrusion)
box = Box(Point((0.,0.,0.)), 20.,20.,20.)
write_fms_file('box', box)
mer = merge( box, extrusion.reverse())
mer.write_obj_file('merge.obj')
print 'writing the obj file for extrusion'
extrusion.write_obj_file('extrusion.obj')
write_fms_file('merge', mer)
vtk_visu(vtk_elements)
write_geo('test_'+str(idtest), geom)

idtest += 1
pretty_print('TEST n.' + str(idtest) + ': QUATERNION')

q1 = Quaternion((1., 0.5, 0.6, 0.9))
print 'Q1 : '
print q1
print 'Q1 CONJUGATE : '
print q1.conjugate
print 'Q1 NORM : '
print abs(q1)
print 'Q1 INVERSE : '
print q1.inverse()
print 'Q1 * Q1^(-1)'
print q1 * q1.inverse()
print 'Q1 -> matrix'
print q1.to_matrix()

mat = np.zeros((3,3))
mat = np.array([ [0., 1., 0. ] , [-1, 0., 0.] , [0., 0., 1.] ])

print '=========================='
print 'MATRIX :'
print mat
print '=========================='
print 'MATRIX -> QUATERNION :'
q =  matrix_to_quaternion(mat)
print q 
print '=========================='
print 'QUATERNION -> MATRIX :'
print q.to_matrix()

q2 = Quaternion((4., 0.25, 0.26, 0.49))
print q1 * q2
