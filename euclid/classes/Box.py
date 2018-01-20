from Point import Point
from Triangulation import Triangulation


"""
                                p0             p1
                              +-------------+
                             /.            /|
                            / .           / |
                           /  .          /  |
     ^                    /   .         /   |
     |                   /    .        /    |
     |              p3  +-------------+ p2  |
     |                  |     .   x C |     |
     |                  |  p4 +.......|.....+ p5
     |         ^        |    .        |    /
 Z   |        /         |   .         |   /
     |       /          |  .          |  /
     |      /           | .           | /
     |     / Y          |.            |/
     |    /         p7  +-------------+  p6
     |   /
     |  /
     | /
     |/
     +---------------------------->
                  X
 """   

class Box(Triangulation):

    def __init__(self, center, Lx, Ly, Lz, **kwargs):
        _points = []
        p0 =  Point( (center[0] - 0.5 * Lx , center[1] + 0.5 * Ly , center[2] + 0.5 * Lz) )
        _points.append(p0)
        p1 =  Point( (center[0] + 0.5 * Lx , center[1] + 0.5 * Ly , center[2] + 0.5 * Lz) ) 
        _points.append(p1)
        p2 =  Point( (center[0] + 0.5 * Lx , center[1] - 0.5 * Ly , center[2] + 0.5 * Lz) ) 
        _points.append(p2)
        p3 =  Point( (center[0] - 0.5 * Lx , center[1] - 0.5 * Ly , center[2] + 0.5 * Lz) ) 
        _points.append(p3)
        p4 =  Point( (center[0] - 0.5 * Lx , center[1] + 0.5 * Ly , center[2] - 0.5 * Lz) ) 
        _points.append(p4)
        p5 =  Point( (center[0] + 0.5 * Lx , center[1] + 0.5 * Ly , center[2] - 0.5 * Lz) ) 
        _points.append(p5)
        p6 =  Point( (center[0] + 0.5 * Lx , center[1] - 0.5 * Ly , center[2] - 0.5 * Lz) ) 
        _points.append(p6)
        p7 =  Point( (center[0] - 0.5 * Lx , center[1] - 0.5 * Ly , center[2] - 0.5 * Lz) ) 
        _points.append(p7)


        _faces = [ [5,4,1], [5,8,4], \
                   [3,6,2], [7,6,3], \
                   [7,3,4], [4,8,7], \
                   [6,1,2], [6,5,1], \
                   [7,8,5], [7,5,6], \
                   [4,3,1], [3,2,1] ]


 
        _physical = [  [1, 'Xmin'] ,\
                       [2, 'Xmax'] ,\
                       [3, 'Ymin'] ,\
                       [4, 'Ymax'] ,\
                       [5, 'Zmin'] ,\
                       [6, 'Zmax'] ]
                   
        _belongs = [1,1,2,2,3,3,4,4,5,5,6,6]
        
        super(Box, self).__init__(points = _points, faces = _faces, physical = _physical, belongs = _belongs)

