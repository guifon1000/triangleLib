from Frame import Frame
from Triangulation import Triangulation
from functions import write_geo
import pygmsh as pg

class Extrusion(Triangulation):
    def __init__(self, polyline_2D, frames, **kwargs):
        geom = pg.built_in.Geometry()
        Nslices = len(frames)
        # get the first polyline
        first = polyline_2D.to_frame(frames[0])
        for i in range(Nslices - 1):
            second = polyline_2D.to_frame(frames[i+1])
            print '---------------'
            print '---------------'
            print first
            print 'vvvvvvvvvvvvvvv'
            print second
            print '---------------'
            print '---------------'
            #for j in range(len(first)):
            #    print str(first[j])+'   <--->   '+str(second[j])
            first = second
        write_geo('test_extrusion', geom)
    #frames.append(pol)

#for pol in frames :
#    pol.pop_to_geom(geom)
#    vtk_elements.append(pol)
