import triangleLib as tl
import numpy as np
import sympy as sp
import quaternions as qtr




def diagBlock(Bloc,N):
    b = Bloc.shape[0]
    m = np.zeros((b*N,b*N),dtype=float)
    for i in range(N):
        idep = i*b
        for j in range(3):
            for k in range(3):
                m[idep+j,idep+k]=Bloc[j,k]
    return m



class Node(tl.Point):
    def __init__(self , x, y, z,ddl  = 'fixed',index=None):
        super(Node,self).__init__(x,y,z)
        self.ddl = ddl
    def setPosition(self, pos ):
        self[0]=pos[0]        
        self[1]=pos[1]        
        self[2]=pos[2]        

class Bar(tl.Segment):
    def __init__(self , pt1 , pt2 ,stif = 1.0,index=None):
        super(Bar,self).__init__(pt1,pt2) 
        self.stif = stif
        ux = tl.Vector(pt2[0]-pt1[0],pt2[1]-pt1[1],pt2[2]-pt1[2])
        ux.normalize()
        self.Ux = ux
        self.Uz = tl.Vector(0.,0.,1.)
        self.Uy = tl.cross(self.Uz,self.Ux)
        matrix = np.zeros((3,3),dtype = float)
        matrix[:,0] = [self.Ux[i] for i in range(3)]
        matrix[:,1] = [self.Uy[i] for i in range(3)]
        matrix[:,2] = [self.Uz[i] for i in range(3)]
        # local stiffness
        k=np.zeros((2,2),dtype=float)
        k=[[1,0,0,-1,0,0],\
           [0,0,0,0,0,0],\
           [0,0,0,0,0,0],\
           [-1,0,0,1,0,0],\
           [0,0,0,0,0,0],\
           [0,0,0,0,0,0]]
        print '----------------'
        print np.dot(np.array(k),self.stif)
        print '----------------'
        print np.linalg.inv(matrix)
        self.matrix = matrix
        mp = diagBlock(self.matrix,2)
        print '===========================' 
        print '============?????==========' 
        print np.dot(mp,k)
        print '===========================' 
        quat = qtr.Quaternion(1.,0.,0.,0.)
         

if __name__ == '__main__':
    L = 1.0
    nd1 = Node(0,0,0,ddl = 'fixed')
    nd2 = Node(0,L,0,ddl = 'fixed')
    nd3 = Node(L,L,0, ddl = ['x','y'])
    b1 = Bar(nd2,nd3)
    b2 = Bar(nd1,nd3)

    print np.dot(np.linalg.inv(b2.matrix),[0.1,0.1,0])


