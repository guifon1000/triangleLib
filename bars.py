import triangleLib as tl
import numpy as np
import sympy as sp
import quaternions as qtr
import matplotlib.pyplot as plt



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
    def __init__(self , pt1 =  None, pt2 = None,stif = 1.0,index=None):
        super(Bar,self).__init__(pt1,pt2) 
        self.stif = stif
        ux = tl.Vector(pt2[0]-pt1[0],pt2[1]-pt1[1],pt2[2]-pt1[2])
        ux.normalize()
        self.Ux = ux
        self.Uz = tl.Vector(0.,0.,1.)
        self.Uy = tl.cross(self.Uz,self.Ux)
        matrix = np.zeros((3,3),dtype = float)
        matrix[0,:] = [self.Ux[i] for i in range(3)]
        matrix[1,:] = [self.Uy[i] for i in range(3)]
        matrix[2,:] = [self.Uz[i] for i in range(3)]
        self.matrix = matrix
        # local stiffness
        k=np.zeros((2,2),dtype=float)
        k=[[1,0,0,-1,0,0],\
           [0,0,0,0,0,0],\
           [0,0,0,0,0,0],\
           [-1,0,0,1,0,0],\
           [0,0,0,0,0,0],\
           [0,0,0,0,0,0]]


	#print '\n\n\n # # # # # # NEW ELEMENT # # # # # # # #'
        #print '======= Rigidity Matrix (local) ========' 
        rm = np.dot(np.array(k),self.stif)
	#print rm
	#print '============ Passage Matrix ============'     
        mp = diagBlock(self.matrix,2)
	#print mp
	#print '======= Rigidity Matrix (global) ======='    
        mpt=np.transpose(mp)
	self.rml = np.dot(mpt,np.dot(rm,mp))
         

def createGlobalStiffnessMatrix(nodes,elements):
 
    gg = np.zeros((3*len(nodes),3*len(nodes)),dtype = float)
    for e in elements:
        n1 = nodes[e[0]]
        n2 = nodes[e[1]]
	x = np.zeros((3*len(nodes),3*len(nodes)),dtype = float)
	b=Bar(n1,n2,stif = kk)
	y = b.rml
	x[e[0]*3:e[0]*3+3,e[1]*3:e[1]*3+3]=y[0:3,3:6]
	x[e[0]*3:e[0]*3+3,e[0]*3:e[0]*3+3]=y[0:3,0:3]
	x[e[1]*3:e[1]*3+3,e[1]*3:e[1]*3+3]=y[3:6,3:6]
	x[e[1]*3:e[1]*3+3,e[0]*3:e[0]*3+3]=y[3:6,0:3]
        gg+=x
    return gg




def createMatrixDDL(nodes,gg):
    ntot = 0
    ddl = []
    for i,n in enumerate(nodes):
        if n.ddl is not 'fixed' and isinstance(n.ddl,list):
	    ntot+=len(n.ddl)
            li = []
	    li.append(i)
	    for d in n.ddl:
	        li.append(d)
            ddl.append(li)
    indices = []
    for li in ddl:
        ipos = li[0]*3
	if 'x' in li[1:]:
	    indices.append(ipos)
	if 'y' in li[1:]:
	    indices.append(ipos+1)
	if 'z' in li[1:]:
	    indices.append(ipos+2)
    mat = np.zeros((ntot,ntot),dtype = float)
    for i in range(ntot): 
        for j in range(ntot): 
            mat[i,j]=gg[indices[i],indices[j]]
    return indices,mat 



def plot2D(nodes,elements,positions,col):
    x=[]
    y=[]
    z=[]
    for i in range(len(positions)/3):
        x.append(positions[3*i])
        y.append(positions[3*i+1])
        z.append(positions[3*i+2])
    plt.scatter(x,y)
    for e in elements:
        xl = []
        yl = []
        xl.append(x[e[0]])
        xl.append(x[e[1]])
        yl.append(y[e[0]])
        yl.append(y[e[1]])
        plt.plot(xl,yl,col)
    plt.axis('equal')

if __name__ == '__main__':
    L = 10.
    nodes = []

    nodes.append(Node(0,L,0,ddl = 'fixed'))
    nodes.append(Node(L,0,0,ddl = 'fixed'))
    nodes.append(Node(0,-L,0,ddl = 'fixed'))
    nodes.append(Node(-L,0,0,ddl = 'fixed'))
    nodes.append(Node(-L/2,L/2.,0,ddl = ['x','y']))    
    nodes.append(Node(L/2,-L/2.,0,ddl = ['x','y']))    
    nodes.append(Node(0,0,0,ddl = ['y']))    
    
    
    
    positions = np.zeros((3*len(nodes)),dtype = float)
    for i,p in enumerate(nodes):
        positions[i*3:i*3+3]=[p.x,p.y,p.z]
    kk = 2.1e11*0.0001/L 
    print 'positions :' 
    print positions

    elements = [[0,4],[1,4],[2,5],[3,5],[4,6],[6,5]]


    f= np.array([0,0,0,0,-10000])
    print f
    gg=createGlobalStiffnessMatrix(nodes,elements)
    print '---------------------------------'
    #print gg
    indices,mat = createMatrixDDL(nodes,gg)
    dim = np.shape(mat)[0]
    print mat
    plt.matshow(mat)
    plt.colorbar()
    plt.show()
    plt.clf()
    print np.linalg.det(mat)
    plot2D(nodes,elements,positions,'k')
    dp =  np.dot(np.linalg.inv(mat),f)
    for i,j in enumerate(indices):
        positions[j]+=dp[i]
        
    print positions
    plot2D(nodes,elements,positions,'r')
    plt.show() 

    


    #print np.dot(np.linalg.inv(b2.matrix),[0.1,0.1,0])


