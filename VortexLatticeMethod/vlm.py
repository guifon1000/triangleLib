import numpy as np
import triangleLib as tl
import matplotlib.pyplot as plt
import fvlm

class fluidPoint(tl.Point):
    def __init__(self,x,y,z):
        super(fluidPoint, self).__init__(x,y,z)
        self.vel=[0.,0.,0.]
        self.magVel=0.

    def computeVelocity(self,others):
        self.vel=[0.,0.,0.]
        self.magVel=0.
        for other in others:
            uAdd=other.getVelocity(self.x,self.y,self.z)
            self.vel[0] += uAdd[0]
            self.vel[1] += uAdd[1]
            self.vel[2] += uAdd[2]
        self.magVel=np.sqrt(self.vel[0]**2+self.vel[1]**2+self.vel[2]**2)
    def invD(self,X,Y,Z):
        if ((X-self.p.x)**2+(Y-self.p.y)**2+(Z-self.p.z)**2)!=0.:
            invD = 1./np.sqrt((X-self.p.x)**2+(Y-self.p.y)**2+(Z-self.p.z)**2)
        else:invD=0.
        return invD
    def invD2(self,X,Y,Z):
        return self.invD(X,Y,Z)**2.




class freeStream(fluidPoint):
    def __init__(self,ux,uy,uz):
        super(fluidPoint, self).__init__(0.,0.,0.)
        self.v=tl.Vector(ux,uy,uz)
    def getVelocity(self,X,Y,Z):
        return self.v




class Source(fluidPoint):
    def __init__(self,strength,x,y,z):
        super(Source, self).__init__(x,y,z)
        self.strength=strength
    def getVelocity(self,X,Y,Z):
        """Returns the velocity field generated by a source/sink.
    
        Arguments
        ---------
        """
        invD3=(self.invD(X,Y,Z))**3.
        u = self.strength/(4*np.pi)*(X-self.p.x)*invD3
        v = self.strength/(4*np.pi)*(Y-self.p.y)*invD3
        w = self.strength/(4*np.pi)*(Z-self.p.z)*invD3
        return tl.Vector(u, v, w)


class lineVortex:
    def __init__(self,A,B,gamma):
        self.p1=A
        self.p2=B
        self.gamma=gamma
    def getVelocity(self,x,y,z):
        #print '-------------\nfortran :'
        v = fvlm.fvlm.vortxl(x,y,z,self.p1,self.p2,self.gamma)
        #print v 
        # Katz
        A=self.p1
        B=self.p2
        r0 = tl.Vector(B[0]-A[0],B[1]-A[1],B[2]-A[2])
        r1 = tl.Vector(x-A[0],y-A[1],z-A[2])
        r2 = tl.Vector(x-B[0],y-B[1],z-B[2])
        pv = tl.cross(r1,r2,norm=False)
        ps1 = tl.dot(r0,r1)
        ps2 = tl.dot(r0,r2)

        normVP2 = pv.norm0**2.
        eps=1.e-9 
        if (r1.norm0<eps or r2.norm0<eps or normVP2<eps):
            K=0.
        else:
            K = (self.gamma/(4.*np.pi*normVP2))*((ps1/r1.norm0)-(ps2/r2.norm0))
        #print 'python :'
        #print tl.Vector(K*pv.x,K*pv.y,K*pv.z)
        return tl.Vector(v[0],v[1],v[2])




class dipolePanel(tl.Panel):
    def __init__(self,p0,p1,p2,p3):
        super(dipolePanel,self).__init__(p0,p1,p2,p3)
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.pan = tl.Panel(p0,p1,p2,p3)
        self.cg = self.pan.cg
        self.n = self.pan.ng
        self.intensity = 0.
    def setIntensity(self,intensity):
        self.intensity = intensity
        self.l0 = lineVortex(self.p0,self.p1,self.intensity)
        self.l1 = lineVortex(self.p1,self.p2,self.intensity)
        self.l2 = lineVortex(self.p2,self.p3,self.intensity)
        self.l3 = lineVortex(self.p3,self.p0,self.intensity)
    def influence(self,other):
        l0 = lineVortex(self.p0,self.p1,1.)
        l1 = lineVortex(self.p1,self.p2,1.)
        l2 = lineVortex(self.p2,self.p3,1.)
        l3 = lineVortex(self.p3,self.p0,1.)
        p = other.cg
        vel = l0.getVelocity(p.x,p.y,p.z)+\
                l1.getVelocity(p.x,p.y,p.z)+\
                l2.getVelocity(p.x,p.y,p.z)+\
                l3.getVelocity(p.x,p.y,p.z)
        return tl.dot(vel,other.n)
    def getVelocity(self,x,y,z):
        v = fvlm.fvlm.voring(x,y,z,self.p0,self.p1,self.p2,self.p3,self.intensity)
        #print '--------- \nvoring fortran : '
        #print v
        #l0=self.l0
        #l1=self.l1
        #l2=self.l2
        #l3=self.l3
        #ux = l0.getVelocity(x,y,z)[0]+l1.getVelocity(x,y,z)[0]+l2.getVelocity(x,y,z)[0]+l3.getVelocity(x,y,z)[0]
        #uy = l0.getVelocity(x,y,z)[1]+l1.getVelocity(x,y,z)[1]+l2.getVelocity(x,y,z)[1]+l3.getVelocity(x,y,z)[1]
        #uz = l0.getVelocity(x,y,z)[2]+l1.getVelocity(x,y,z)[2]+l2.getVelocity(x,y,z)[2]+l3.getVelocity(x,y,z)[2]
        #print 'python :'
        #print (ux,uy,uz)
        return tl.Vector(v[0],v[1],v[2])

class horseShoeVortex:
    def __init__(self,p0,p1):
        self.p0=p0
        self.p1=p1





class Surface:
    def __init__(self,tab,sf):   # sf is a bigQuad
        self.tab = tab
        self.bigquad = sf
        self.cell_data = {}
        self.point_data = {}

    def dipoleMatrix(self,tab):
        mat=np.zeros((len(tab),len(tab)))
        for i,pi in enumerate(tab):  
            for j,pj in enumerate(tab):  
                mat[i,j] = pj.influence(pi)
        self.M=mat

    def addCellData(self,val,name,quad=True):
        tmp=[]
        for i in range((self.bigquad.Nu-1)*(self.bigquad.Nv-1)):
            tmp.append(val[i])
            tmp.append(val[i])
        self.cell_data[name]=np.array(tmp)
    def writeVTK(self,name,vec):
        meshio.write('test.vtu', self.pts, self.cells,cell_data={'xx' : vec})
        #tv.triangle_faces_to_VTK(name,self.bigquad.X,self.bigquad.Y,self.bigquad.Z,np.array(self.bigquad.faces),None,self.cell_data)
    def plot(self):
        print ''