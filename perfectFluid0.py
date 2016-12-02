import numpy as np
import bidLib as bid
import matplotlib.pyplot as plt





def loadXfoilFile(name, method = 'Dipole'):
    f = open(name,'r').readlines()
    points = []
    for i in range(1,len(f)):
        l = f[i].split()
        p = PointElement(float(l[0]),3.*float(l[1]))
        points.append(p)
    panels = []
    for i in range(len(points)-1):
        p0 = points[i]
        p1 = points[i+1]
        if method=='Dipole':
            panels.append(DipolePanel(p0,p1,1.))
        if method=='Source':
            panels.append(SourcePanel(p0,p1,1.))
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
    try :
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

    except:
        return np.dot(self.pmat,np.array([x-self.origin[0],y-self.origin[1]]))
    return



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
    def plot(self,plt,axis = False):
        plt.scatter(self.origin[0],self.origin[1],marker='*')

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
        xm = 0.5*(self.xa+self.xb)
        ym = 0.5*(self.ya+self.yb)
        self.middle = [xm,ym]
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
    def plot(self,plt,axis = False):
        for p in self:
            plt.plot([self.xa,self.xb],[self.ya,self.yb],'--k')
            if axis:
                xloc = [1.,0.]
                yloc = [0., 1.0]
                xglob = np.dot(self.ipmat,xloc)
                yglob = np.dot(self.ipmat,yloc)
                plt.plot([self.origin[0],self.origin[0]+xglob[0]],\
                         [self.origin[1],self.origin[1]+xglob[1]],\
                         '*k')
                plt.plot([self.origin[0],self.origin[0]+yglob[0]],\
                         [self.origin[1],self.origin[1]+yglob[1]],\
                         '*k')




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
            return [ut, un]
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
        un = 0.5*fac * (dx**2.-dy**2.) /(dy**2.+dx**2.)**2.
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un]


class DipolePanel(LinearElement):
    def __init__(self,p0,p1,mu=1.):
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
        ut = fac * ( (dy/(dx1**2. + dy**2.)) - (dy/(dx2**2. + dy**2.)) )
        un = fac * ( (dx1/(dx1**2. + dy**2.)) - (dx2/(dx2**2. + dy**2.)) )
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un] 
    def autoInflu(self):
        self.autoInf = 2./(np.pi*self.length) 
        return 2./(np.pi*self.length)

class SourcePanel(LinearElement):
    def __init__(self,p0,p1,mu=0.):
        super(SourcePanel,self).__init__(p0,p1)
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
        ut = 0.5 * fac * np.log(   ( dx1**2. + dy**2.)/(dx2**2.+dy**2.))
        un = fac * (  np.arctan2(dy,dx2) - np.arctan2(dy,dx1) )
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un] 
    def autoInflu(self):
        self.autoInf = 0.5
        return 0.5

class PanelArray(list):
    def __init__(self,tab):
        super(PanelArray,self).__init__(tab)
    def influenceMatrix(self):
        print 'initializing the influence matrix'


class Particule(Vortex):
    def __init__(self,x,y,gamma=0.,t=0.):
        super(Particule,self).__init__(x,y,gamma)
        self.time = t
        self.vel = np.array([0.,0.])
    def influence(self,panel):
        return np.dot(self.velocities(panel.middle[0],panel.middle[1]),panel.pmat)
    
    
    
def particuleVelocity(part,x,y):
    for (i,p) in enumerate(part):
        return [p.influence(x,y)[0],p.influence(x,y)[1]]


Nx = 400
Ny = 400
eps = 0.00001
borne = 1.2
Px = np.linspace(-0.05,1.,num = Nx)
Py = np.linspace(-0.15,0.15,num = Ny)
X,Y = np.meshgrid(Px,Py)
vinf = np.array([1.,0.01])

mode = 'mano'
mode = 'wake'
singul = 'Dipole'
if mode == 'mano':
    panels = []
    start = PointElement(0.,-0.1)
    Nl = 1
    L=0.2
    p0 = start
    for i in range(1,Nl+1):
        p = PointElement(start[0],start[1]+float(i)*L/float(Nl))
        if singul == 'Dipole': pan = DipolePanel(start,p)
        if singul == 'Source': pan = SourcePanel(start,p)
        p0 = p
        panels.append(pan)
if mode =='xfoil':
    panels = loadXfoilFile('NACAcamber0012.dat',singul)
if mode =='wake':
    panels = [] 
    A = PointElement(0.,0.05)
    B = PointElement(0.4,-0.05)

    C = PointElement(0.3,0.05)
    D = PointElement(0.7,-0.1)
    if singul == 'Dipole': 
        pan = DipolePanel(A,B)
        panels.append(pan)
        pan = DipolePanel(B,C)
        #panels.append(pan)
        pan = DipolePanel(C,D)
        #panels.append(pan)
    if singul == 'Source': 
        pan = SourcePanel(A,B)
        panels.append(pan)
        pan = SourcePanel(B,C)
        panels.append(pan)



dt = 0.1

tredge = PointElement(panels[-1][-1][0],panels[-1][-1][1])
Pl = PointElement(panels[-1][-1][0]+vinf[0]*dt,panels[-1][-1][1]+vinf[1]*dt)


particules = []

lastPanel = DipolePanel(tredge,Pl)
panels.append(lastPanel)
       
N = len(panels)
M = np.zeros((N-1,N-1),dtype=float)
for i,pi in enumerate(panels[:-1]):
    
    for j,pj in enumerate(panels[:-1]):
        if i==j:
            M[i,j] = pi.autoInflu()
        else:
            M[i,j] = pj.velocities(pi.middle[0],pi.middle[1])[1]
        
u , v = np.zeros_like(X),np.zeros_like(X)
plt.clf()

particules.append(Vortex(Pl[0],Pl[1],1.))

for it in range(2):
    vinf2 = vinf
    Npart = len(particules)
    Mpart = np.zeros((Npart,Npart),dtype=float)

   # advection
    for ip,pi in enumerate(particules):
        pi.vel = [vinf[0],vinf[1]]
        for pa in panels:
            uu = pa.velocities(pi[0],pi[1],mode = 'global')
            pi.vel[0]+=uu[0]
            pi.vel[1]+=uu[1]
        for jp,pj in enumerate(particules) :
            if ip!=jp:
                uu = pj.velocities(pi[0],pi[1],mode = 'global')
                pi.vel[0]+=uu[0] 
                pi.vel[1]+=uu[1] 
        pi[0]+=pi.vel[0]*dt
        pi[1]+=pi.vel[1]*dt

    b = np.zeros(N-1,dtype = float)
    for i,p in enumerate(panels[:-1]):
        vv = vinf2
        for pp in particules:
            up = pp.velocities(p.middle[0],p.middle[1],mode = 'global')
            vv[0] += up[0]
            vv[1] += up[1]
        b[i] = -(vv[0]*p.n[0]+vv[1]*p.n[1])


    stren = np.linalg.solve(M,b)
    for i,p in enumerate(panels[:-1]):
        p.mu = stren[i]

    particules.append(Particule(Pl[0],Pl[1],1.))
    particules[-1].vel = vinf


    u[:] = vinf[0]
    v[:] = vinf[1]
    for i,p in enumerate(panels):
        if p!= panels[-1]:
            p.plot(plt,axis = False)
        else:
            plt.scatter(p.xb,p.yb) 
        uloc = p.velocities(X,Y,mode='global')
        u += uloc[0]
        v += uloc[1]
    for i,p in enumerate(particules):
        plt.scatter(p[0],p[1]) 
        uloc = p.velocities(X,Y,mode='global')
        u += uloc[0]
        v += uloc[1]
    plt.axis('equal')
    plt.streamplot(X, Y, u,v, density=4, linewidth=1, arrowsize=1, arrowstyle='->')
    plt.show()

    print "zob"





if singul == 'Dipole':
    b = [-(vinf[0]*p.n[0]+vinf[1]*p.n[1]) for p in panels[:-1]]
elif singul == 'Source':
    b = [-(vinf[0]*p.n[0]+vinf[1]*p.n[1]) for p in panels[:-1]]



stren = np.linalg.solve(M,b)
for i,p in enumerate(panels[:-1]):
    p.mu = stren[i]
plt.clf()


panels[-1].mu=stren[-1]

lastVec = panels[-1].t


#panels.append(Vortex(lastPanel[0]+fac*lastVec[0],lastPanel[1]+fac*lastVec[1],stren[-1]*-1.))


u[:] = vinf[0]
v[:] = vinf[1]
for i,p in enumerate(panels):
    if p!= panels[-1]:
         p.plot(plt,axis = False)
    else:
         plt.scatter(p.xb,p.yb) 
    uloc = p.velocities(X,Y,mode='global')
    u += uloc[0]
    v += uloc[1]
plt.axis('equal')
plt.streamplot(X, Y, u,v, density=4, linewidth=1, arrowsize=1, arrowstyle='->')
plt.show()


print '--------'