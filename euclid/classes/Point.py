
class Point(list):
    def __init__(self,*largs, **kwargs):
        """
            class Point -> array of three floats (the coordinates), only cartesian
            also has the attributes .x , .y and .z --> same as self[0], self[1], self[2]
            can have an index
        """
        super(Point,self).__init__(*largs)
