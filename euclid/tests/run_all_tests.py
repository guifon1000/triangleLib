import sys
sys.path.append('../')
from classes.Point import Point
from classes.Triangle import Triangle
from classes.Triangulation import Triangulation
from classes.Vector import Vector
from classes.Sphere import Sphere
from classes.Frame import Frame
from classes.Frame import Frame0
import pygmsh as pg
from modelers.profiles.splineProfileMultiParam import  Profile
from modelers.planet.Planet import  Planet
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from functions import parameter_frame

def write_geo(name,geom):
    fg = open(name+'.geo','w')
    for l in geom.get_code():
        fg.write(l)
    fg.close()


geom = pg.built_in.Geometry()



# create three points
p1 = Point((0.001,0.004,0.))
p2 = Point((0.23,1.,0.))
p3 = Point((1.,0.4,0.7))


idtest = 1

points = [p1, p2, p3]

geo_points = []
#for p in points: 
#    geo_points.append(p.pop_to_geom(geom))


#print [g.id for g in geo_points]


print "##################################################################################"
print "########################### TEST n."+str(idtest)+": TRIANGLE           #########################"
print "##################################################################################"
# create a triangle
t = Triangle(points)

triangles = [t]
#for tr in triangles: tr.pop_to_geom(geom)
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

#t = -t

idtest += 1

print "##################################################################################"
print "############################# TEST n."+str(idtest)+": VECTOR           #########################"
print "##################################################################################"

v = t.normal
print "raw vector :"+str(v)
print "norm 2 : "+str(v.norm2)
print "norm 1 : "+str(v.norm)
print "unit vector : "+str(v.unit())



idtest += 1
print "##################################################################################"
print "############################# TEST n."+str(idtest)+": IMPORT JSON TRIANGULATION       #############"
print "##################################################################################"


d = Triangulation()
d.load_file('../samples/icosahedron.json')
print d.cg


idtest += 1
print "##################################################################################"
print "############################# TEST n."+str(idtest)+": SPHERES       #############"
print "##################################################################################"


d.reorient_convex()
d.translate((1.,0.,-2.))
d = d.refine_2()
d = d.refine_2()
d = d.refine_2()
d.write_obj_file('refined_icosahedron.obj')
s1 = Sphere(refin = 0, center = (1.,1.,1.), radius = 0.3)
s2 = Sphere(refin = 1, center = (0.,-4.,2.), radius = 0.1)
s3 = Sphere(refin = 2, center = (-3.,1.,-1.), radius = 0.7)
p1 = Planet(refin = 1, center = (-3.,-5.,-3.), radius = 0.8)
s1.write_obj_file('sphere1.obj')
s2.write_obj_file('sphere2.obj')
s3.write_obj_file('sphere3.obj')
p1.write_obj_file('planet1.obj')


print 'ok'

idtest += 1
print "##################################################################################"
print "############################# TEST n."+str(idtest)+":VECTOR EXTRUSION       #############"
print "##################################################################################"

start_point = Point([0.,0.,0.])
end_point = Point([1.,1.,0.])

mat = np.zeros((3, 3), dtype = float)

alpha = np.pi/6


xf = Vector((np.cos(alpha), np.sin(alpha), 0.))
yf = Vector((-np.sin(alpha), np.cos(alpha), 0.))
zf = Vector((0., 0., 0.))

mat[0] = xf
mat[1] = yf
mat[2] = zf



print mat

# control points of the generatrix
x = [0.0, 0.69, 1.2, 2.7]
y = [0.0, 0.4, 0.8, 1.4]
z = [0.0, 2.9, 3.8, 8.8]

# tck, u represent the parametric 3d curve
tck, u = interpolate.splprep([x,y,z], s=3)
for s in np.linspace(0., 1., num = 50):
    frame = parameter_frame(tck, s, mode = 'frenet')
    pf = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = 12) # creation of the 2d profile
    pol = pf.polyline()
    #pol.pt3d = [Point([0., 0., 0.]) for p in pol  ]
    for i,p in enumerate(pol):
        loc  =  np.dot([0., p[0], p[1]], frame[1])
        print loc
        print 'tt'
        pol.pt3d[i][0] = frame[0][0] + loc[0]
        pol.pt3d[i][1] = frame[0][1] + loc[1]
        pol.pt3d[i][2] = frame[0][2] + loc[2]
        #print str(p)+' --> '+str(pol.pt3d[i])
    pol.pop_to_geom(geom)

write_geo('tests', geom)


1/0
#pf = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = 40) # creation of the 2d profile
#pol = pf.polyline()

#pol.to_frame(end_frame, scale = 5.)
#pol.pop_to_geom(geom)


#pf2 = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = 40) # creation of the 2d profile
#pol2 = pf2.polyline()

#pol2.to_frame(start_frame, scale = 5.)
#pol2.pop_to_geom(geom)



idtest += 1
print "##################################################################################"
print "############################# TEST n."+str(idtest)+":          PROFILE 3D       #############"
print "##################################################################################"

# control points of the generatrix
x = [0.0, 0.1, 0.4, 0.6, 0.9]
y = [0.0, 0.05, 0.09, 0.11, 0.4]
z = [0.0, 0.3, 0.8, 1.2, 2.9]

# tck, u represent the parametric 3d curve
tck, u = interpolate.splprep([x,y,z], s=2)
#ex = Extrusion()


name = 'wingZero'
Nslices = 100 # number of slices
npt = 63 # points of the profile
t = np.linspace(0., 1., Nslices) # parametric space
pf = Profile(typ = 'fon',par = [0.82,0.21,0.13,0.08,0.029],npt = npt) # creation of the 2d profile
fi = Frame0(t[0], tck, type = 'Xnat')
pol = pf.polyline()
pol.to_frame(fi, scale = 0.5)
li = pol.pop_to_geom(geom)

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
    fip1 = Frame0(sip1, tck, type = 'Xnat')
    pol = pf.polyline()
    pol.to_frame(fip1, scale = 0.5*np.cos(sip1*0.4*np.pi))
    lip1 = pol.pop_to_geom(geom)
    for j in range(len(li)):
        lij = li0[j]
        lip1j = lip1[j]
        lti = geom.add_line(lij.points[0], lip1j.points[0])
        ltip1 = geom.add_line(lij.points[1], lip1j.points[1])
        lloop = geom.add_line_loop([lti, lip1j, -ltip1, -lij])
        sf = geom.add_surface(lloop)
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


write_geo('tests', geom)
