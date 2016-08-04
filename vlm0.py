import numpy as np
import triangleLib as tl
import matplotlib.pyplot as plt
import pyevtk.hl as vhl

class fluidPoint:
    def __init__(self,x,y,z):
        self.p=tl.Point(x,y,z)
        self.vel=[0.,0.,0.]
        self.magVel=0.
        self.streamFunction=0.
        self.potential=0.
    def computeVelocity(self,others):
        self.vel=[0.,0.,0.]
        self.magVel=0.
        for other in others:
            uAdd=other.getVelocity(self.p.x,self.p.y,self.p.z)
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
        for k in range(Nz): 
            for j in range(Ny):
                for i in range(Nx): 
                    self.x[i,j,k] = X[i] 
                    self.y[i,j,k] = Y[j]
                    self.z[i,j,k] = Z[k]
        ic=-1
        for k in range(Nz-1):
            for j in range(Ny-1):
                for i in range(Nx-1):
                    ic+=1
                    xc=(1./8.)*(self.x[i,j,k]+self.x[i+1,j,k]+self.x[i,j+1,k]+self.x[i,j,k+1]+self.x[i+1,j+1,k]+self.x[i,j+1,k+1]+self.x[i+1,j,k+1]+self.x[i+1,j+1,k+1])
                    yc=(1./8.)*(self.y[i,j,k]+self.y[i+1,j,k]+self.y[i,j+1,k]+self.y[i,j,k+1]+self.y[i+1,j+1,k]+self.y[i,j+1,k+1]+self.y[i+1,j,k+1]+self.y[i+1,j+1,k+1])
                    zc=(1./8.)*(self.z[i,j,k]+self.z[i+1,j,k]+self.z[i,j+1,k]+self.z[i,j,k+1]+self.z[i+1,j+1,k]+self.z[i,j+1,k+1]+self.z[i+1,j,k+1]+self.z[i+1,j+1,k+1])
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
        #self.d={'ux' : np.array([(p.vel[0],p.vel[1],p.vel[2]) for p in self.cells]).reshape(self.Nx-1,self.Ny-1,self.Nz-1).transpose()}
         
    def writeVTK(self):
        vhl.gridToVTK("./structured", self.x, self.y, self.z, cellData = self.d)
         
        #vhl.pointsToVTK("./points", np.array([p.p.x for p in self.cells]), np.array([p.p.y for p in self.cells]), np.array([p.p.z for p in self.cells]),data=None) 


class freeStream(fluidPoint):
    def __init__(self,ux,uy,uz):
        self.p=tl.Point(0.,0.,0.)
        self.v=[ux,uy,uz]
    def getVelocity(self,X,Y,Z):
        return self.v




class Source(fluidPoint):
    def __init__(self,strength,x,y,z):
        self.p=tl.Point(x,y,z)
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
        return [u, v, w]


class lineVortex:
    def __init__(self,A,B,gamma):
        self.p1=A
        self.p2=B
        self.gamma=gamma
    def getVelocity(self,x,y,z):
        M=fluidPoint(x,y,z)
        rn1=tl.Vector()
        rn1.fromPoints(M.p,self.p1)
        rn2=tl.Vector()
        rn2.fromPoints(M.p,self.p2)
        pv=tl.cross(rn1,rn2,norm=False)
        ps=tl.dot(rn1,rn2)
        
        ux = self.gamma/(4.*np.pi)*pv.x/pv.norm0**2.*(rn1.norm0+rn2.norm0)*(1.-ps)/(rn1.norm0*rn2.norm0)
        uy = self.gamma/(4.*np.pi)*pv.y/pv.norm0**2.*(rn1.norm0+rn2.norm0)*(1.-ps)/(rn1.norm0*rn2.norm0)
        uz = self.gamma/(4.*np.pi)*pv.z/pv.norm0**2.*(rn1.norm0+rn2.norm0)*(1.-ps)/(rn1.norm0*rn2.norm0)
        return [ux,uy,uz]


class Doublet(fluidPoint):
    def __init__(self,kappa,x,y,z,normal):
        self.kappa=kappa
        self.p=tl.Point(x,y,z)
        self.normal=normal

    def getVelocity0(self,X,Y,Z):
        dx=X-self.p.x
        dy=Y-self.p.y
        u=-self.kappa/(2*np.pi)*(dx*dx-dy*dy)/(dx**2+dy**2)**2
        v=-self.kappa/(np.pi)*dx*dy/(dx**2+dy**2)**2
        w=0.
        return [u,v,w]
    def getVelocity(self,X,Y,Z):
        
        P=tl.Point(X,Y,Z)
        Q=tl.Point(self.p.x+self.normal.x,self.p.y+self.normal.y,self.p.z+self.normal.z)
        u=tl.Vector()
        u.fromPoints(self.p,P)
        v=tl.Vector()
        v.fromPoints(self.p,Q)
        angle=tl.angle(v,u)
        print '#################################################################################################'
        print '#                                    3D DOUBLET VELOCITY                                        #'
        print '#                                                                                               #'
        print '      P        = '+str(P)
        print '      Q        = '+str(Q)
        print '      doublet  = '+str(self.p)
        print '      u        = '+str(u)
        print '      v        = '+str(v)
        print '      angle    = '+str(angle)+' rad    ;    '+str(angle*180./(np.pi))+' degrees'
        print '#                                                                                               #'
        print '#                                                                                               #'
        print '#################################################################################################'
        print ''
        return[0.,0.,0.]
    #def get_stream_function(self,X,Y):
        #dx=X-self.x
        #dy=Y-self.y
        #psi=-self.kappa/(2*np.pi)*dy/(dx**2+dy**2)**2
        #return psi



