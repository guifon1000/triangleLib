import numpy as np
import triangleLib as tl
import triangles_to_VTK as tv
import matplotlib.pyplot as plt
import pyevtk.hl as vhl

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




class Mesh:
    def __init__(self,Nx,Ny,Nz,Lx,Ly,Lz):
        X=np.linspace(-0.5*Lx,0.5*Lx,Nx)
        Y=np.linspace(-0.5*Ly,0.5*Lx,Ny)
        Z=np.linspace(-0.5*Lz,0.5*Lz,Nz)

        self.x = np.zeros((Nx, Ny, Nz)) 
        self.y = np.zeros((Nx, Ny, Nz)) 
        self.z = np.zeros((Nx, Ny, Nz))
        self.cells=[]
        for i in range(Nx): 
            for j in range(Ny):
                for k in range(Nz): 
                    self.x[i,j,k] = X[i] 
                    self.y[i,j,k] = Y[j]
                    self.z[i,j,k] = Z[k]
        ic=-1
        for i in range(Nx-1):
            for j in range(Ny-1):
                for k in range(Nz-1):
                    ic+= 1
                    xc = (1./8.)*(self.x[i,j,k]+\
                            self.x[i+1,j,k]+ \
                            self.x[i,j+1,k]+\
                            self.x[i,j,k+1]+\
                            self.x[i+1,j+1,k]+\
                            self.x[i,j+1,k+1]+\
                            self.x[i+1,j,k+1]+\
                            self.x[i+1,j+1,k+1])
                    yc = (1./8.)*(self.y[i,j,k]+\
                            self.y[i+1,j,k]+\
                            self.y[i,j+1,k]+\
                            self.y[i,j,k+1]+\
                            self.y[i+1,j+1,k]+\
                            self.y[i,j+1,k+1]+\
                            self.y[i+1,j,k+1]+\
                            self.y[i+1,j+1,k+1])
                    zc = (1./8.)*(self.z[i,j,k]+\
                            self.z[i+1,j,k]+\
                            self.z[i,j+1,k]+\
                            self.z[i,j,k+1]+\
                            self.z[i+1,j+1,k]+\
                            self.z[i,j+1,k+1]+\
                            self.z[i+1,j,k+1]+\
                            self.z[i+1,j+1,k+1])

                    self.cells.append(fluidPoint(xc,yc,zc)) 

        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
    def compute(self,others):
        for p in self.cells:
            p.computeVelocity(others)
            v=np
        self.d={'ux' : np.array([p.vel[0] for p in self.cells]).reshape((self.Nx-1),(self.Ny-1),(self.Nz-1)),\
                'uy' : np.array([p.vel[1] for p in self.cells]).reshape((self.Nx-1),(self.Ny-1),(self.Nz-1)),\
                'uz' : np.array([p.vel[2] for p in self.cells]).reshape((self.Nx-1),(self.Ny-1),(self.Nz-1))}
         
    def writeVTK(self):
        vhl.gridToVTK("./structured", self.x, self.y, self.z, cellData = self.d)
         


class freeStream(fluidPoint):
    def __init__(self,ux,uy,uz):
        super(fluidPoint, self).__init__(0.,0.,0.)
        self.v=tl.Vector(ux,uy,uz)
    def getVelocity(self,X,Y,Z):
        return self.v




class Source(fluidPoint):
    def __init__(self,strength,x,y,z):
        super(fluidPoint, self).__init__(x,y,z)
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
        # Katz
        A=self.p1
        B=self.p2
        r0 = tl.Vector(B.x-A.x,B.y-A.y,B.z-A.z)
        r1 = tl.Vector(x-A.x,y-A.y,z-A.z)
        r2 = tl.Vector(x-B.x,y-B.y,z-B.z)
        pv = tl.cross(r1,r2,norm=False)
        ps1 = tl.dot(r0,r1)
        ps2 = tl.dot(r0,r2)

        normVP2 = pv.norm0**2.
        eps=1.e-9 
        if (r1.norm0<eps or r2.norm0<eps or normVP2<eps):
            K=0.
        else:
            K = (self.gamma/(4.*np.pi*normVP2))*((ps1/r1.norm0)-(ps2/r2.norm0))
        return tl.Vector(K*pv.x,K*pv.y,K*pv.z)




class dipolePanel:
    def __init__(self,p0,p1,p2,p3):
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
        l0 = lineVortex(self.p0,self.p1,self.intensity)
        l1 = lineVortex(self.p1,self.p2,self.intensity)
        l2 = lineVortex(self.p2,self.p3,self.intensity)
        l3 = lineVortex(self.p3,self.p0,self.intensity)
        ux = l0.getVelocity(x,y,z)[0]+l1.getVelocity(x,y,z)[0]+l2.getVelocity(x,y,z)[0]+l3.getVelocity(x,y,z)[0]
        uy = l0.getVelocity(x,y,z)[1]+l1.getVelocity(x,y,z)[1]+l2.getVelocity(x,y,z)[1]+l3.getVelocity(x,y,z)[1]
        uz = l0.getVelocity(x,y,z)[2]+l1.getVelocity(x,y,z)[2]+l2.getVelocity(x,y,z)[2]+l3.getVelocity(x,y,z)[2]
        return tl.Vector(ux,uy,uz)

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
    def writeVTK(self,name):
        tv.triangle_faces_to_VTK(name,self.bigquad.X,self.bigquad.Y,self.bigquad.Z,np.array(self.bigquad.faces),None,self.cell_data)
        
    
