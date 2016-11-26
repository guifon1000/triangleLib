import numpy as np
import bidLib as bid
import matplotlib.pyplot as plt





def loadXfoilFile(name):
    f = open(name,'r').readlines()
    points = []
    for i in range(1,len(f)):
        l = f[i].split()
        p = PointElement(float(l[0]),float(l[1]))
        points.append(p)
    panels = []
    for i in range(len(points)-1):
        p0 = points[i]
        p1 = points[i+1]
        panels.append(DipolePanel(p0,p1,10.))
    plt.clf()
    for p in panels:
        plt.plot([p.xa,p.xb],[p.ya,p.yb])
        pt = p.globXY(0.001,0.001)
        plt.scatter(pt[0],pt[1])
    plt.axis('equal')
    plt.show() 
    return panels


def orientation(el,vec,orient = 'x'):
    if orient == 'x':
        el.angleX = bid.angle([1.,0.],vec)
    elif orient =='y':
        el.angleX = 0.5*np.pi - bid.angle(vec,[0.,1.])
    pmat = np.zeros((2,2),dtype = float)
    pmat[0][0] = np.cos(el.angleX)
    pmat[0][1] = np.sin(el.angleX) 
    pmat[1][0] = -np.sin(el.angleX)
    pmat[1][1] = np.cos(el.angleX)
    el.pmat = pmat
    el.ipmat = np.linalg.inv(el.pmat)

    el.t = [el.pmat[0,0],el.pmat[0,1]  ]
    el.n = [el.pmat[1,0],el.pmat[1,1]  ]



def localPos(self,x,y):
    Nx = np.shape(x)[1]
    Ny = np.shape(x)[0]
    xl = np.zeros_like(x)
    yl = np.zeros_like(y)
    for i in range(Ny): 
        for j in range(Nx):
            prod = np.dot(self.pmat,np.array([x[i][j]-self.origin[0],y[i][j]-self.origin[1]])) 
            xl[i,j] = prod[0]
            yl[i,j] = prod[1]
    return xl,yl

def globalVel(self,ut,un):
    ux = np.zeros_like(ut)
    uy = np.zeros_like(un)
    for i in range(Ny):
        for j in range(Nx):
            prod = np.dot(self.ipmat,np.array([ut[i][j],un[i][j]]))
            ux[i,j] = prod[0]
            uy[i,j] = prod[1]
    return [ux,uy]

class PointElement(list):
    def __init__(self,x,y):
        super(PointElement,self).__init__([x,y])
        self.origin = [self[0],self[1]]


class LinearElement(list):
    def __init__(self,p0,p1):
        super(LinearElement,self).__init__([p0,p1])
        u=[p1[0]-p0[0],p1[1]-p0[1]]
        orientation(self,u)
        self.xa = p0[0]
        self.ya = p0[1]
        self.xb = p1[0]
        self.yb = p1[1]
        self.length = np.sqrt((self.xb-self.xa)**2.+(self.yb-self.ya)**2.)
        self.origin = [p0[0],p0[1]]
    def locXY(self,x,y):
        ori = np.array([self[0][0],self[0][1]])
        prod = np.dot(self.pmat,np.array([x-ori[0],y-ori[1]]))
        out = prod
        return out
    def globXY(self,x,y):
        ori = np.array([self[0][0],self[0][1]])
        prod = np.dot(self.ipmat,np.array([x,y]))
        out = [prod[0]+ori[0],prod[1]+ori[1]]
        return out

class Vortex(PointElement):
    def __init__(self,x,y,gamma=0.):
        super(Vortex,self).__init__(x,y)
        self.gamma = gamma
        orientation(self,[1.,0.])
    def velocities(self,x,y,mode='local'):
        fac = self.gamma/(2.*np.pi)
        dx = x - self[0]
        dy = y - self[1]
        ut = fac * dy /(dy**2.+dx**2.)
        un = -fac * dx /(dy**2.+dx**2.)
        if mode == 'global':
            pass
        return [ut, un]
            


class Dipole(PointElement):
    def __init__(self,x,y,normal,gamma=0.):
        super(Dipole,self).__init__(x,y)
        self.gamma = gamma
        self.normal = normal
        orientation(self,self.normal,orient='y')
    def velocities(self,x,y,mode='local'):
        """
        warning : x,y in GLOBAL  COORDINATES
        """
        xl,yl = localPos(self,x,y)
        fac = self.gamma/(np.pi)
        dx = xl
        dy = yl
        ut = fac * dx*dy /(dy**2.+dx**2.)**2.
        un = -0.5*fac * (dx**2.-dy**2.) /(dy**2.+dx**2.)**2.
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un]


class DipolePanel(LinearElement):
    def __init__(self,p0,p1,mu=0.):
        super(DipolePanel,self).__init__(p0,p1)
        self.mu = mu
    def velocities(self,x,y,mode='local'):
        """
        warning : x,y in GLOBAL  COORDINATES
        """
        xl,yl = localPos(self,x,y)
        fac = self.mu/(2.*np.pi)
        dx1 = xl
        dx2 = xl - self.length
        dy = yl
        ut = -fac * ( (dy/(dx1**2. + dy**2.)) - (dy/(dx2**2. + dy**2.)) )
        un = fac * ( (dx1/(dx1**2. + dy**2.)) - (dx2/(dx2**2. + dy**2.)) )
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un] 


b = []
panels = []
A = PointElement(-0.5,-0.5)
B = PointElement(0.5,0.5)

dp = DipolePanel(A,B,1.)


pan = Vortex(0.,0.,1.)
#panels.append(pan)
Nx = 250
Ny = 230
eps = 0.00001
borne = 1.2
Px = np.linspace(-borne+eps,borne+eps,num = Nx)
Py = np.linspace(-borne+eps,borne+eps,num = Ny)
X,Y = np.meshgrid(Px,Py)
vinf = np.array([1.,0.])
u , v = np.zeros_like(X),np.zeros_like(X)
print ' dipole ==============='
d = Dipole(0.,0.,[-0.1,1.],gamma = 1.0)
#panels.append(d)
print ' ======================'

d = Dipole(0.,0.,[-0.1,1.],gamma = 1.0)
plt.plot
dp = DipolePanel(A,B,1.414)

panels.append(dp)


panels = loadXfoilFile('profiles/NACAcamber0012.dat')
u[:] = vinf[0]
v[:] = vinf[1]
for p in panels:
    uloc = p.velocities(X,Y,mode='global')
    plt.plot([p.xa,p.xb],[p.ya,p.yb])
    u += uloc[0]
    v += uloc[1]
plt.axis('equal')
plt.streamplot(X, Y, u,v, density=4, linewidth=1, arrowsize=1, arrowstyle='->')
plt.show()


print b
print '--------'
