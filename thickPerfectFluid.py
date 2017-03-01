import numpy as np
import bidLib as bid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess as sub




def loadXfoilFile(name, method = 'Dipole',alpha = 0.,scale = 1.):
    pmat = np.zeros((2,2),dtype = float)
    pmat[0][0] = np.cos(alpha)*scale
    pmat[0][1] = np.sin(alpha)*scale
    pmat[1][0] = -np.sin(alpha)*scale
    pmat[1][1] = np.cos(alpha)*scale
    f = open(name,'r').readlines()
    points = []
    for i in range(1,len(f)):
        l = f[i].split()
        locP = np.array([float(l[0]),3.*float(l[1])])
        globP = np.dot(pmat,locP)
        p = PointElement(globP[0],globP[1])
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
            plt.plot([self.xa,self.xb],[self.ya,self.yb],'sk',markersize = 3.)
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
            
class VortexParticule(Vortex):
    def __init__(self,x,y,gamma=0.,vinf = [0.,0.]):
        super(VortexParticule,self).__init__(x,y,gamma)
	orientation(self,[1.,0.])
	self.x = x
	self.y = y
        self.speed = [vinf[0],vinf[1]]


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
        # - + : like in katz
        ut = -fac * ( (dy/(dx1**2. + dy**2.)) - (dy/(dx2**2. + dy**2.)) )
        un = fac * ( (dx1/(dx1**2. + dy**2.)) - (dx2/(dx2**2. + dy**2.)) )
        if mode == 'global':
            return globalVel(self,ut,un)
        else:
            return [ut, un] 
    def autoInflu(self):
        self.autoInf = -2./(np.pi*self.length) 
        return -2./(np.pi*self.length)

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
        N = len(panels)
        M = np.zeros((N,N),dtype=float)
        for i,pi in enumerate(panels):
            pi.mu = 1. 
            for j,pj in enumerate(panels):
                pj.mu = 1.0
                if i==j:
                    M[i,j] = -pi.autoInflu()
                else:
                    vv = np.array(pj.velocities(pi.middle[0],pi.middle[1]))
                    p1 = np.dot(vv,pj.ipmat)
                    p2 = np.dot(p1,pi.pmat)
                    M[i,j] = p2[1]
        return  M



    def normalVel(self,vinf,others = None):
        nv = np.zeros(len(self),dtype = float)
	for i,pi in enumerate(self):
	    mid = pi.middle
	    nv[i] = (vinf[0]*pi.n[0]+vinf[1]*pi.n[1])
            for j,pj in enumerate(self):
                if i!= j:
                    vv = np.array(pj.velocities(mid[0],mid[1]))
                    p1 = np.dot(vv,pj.ipmat)
                    p2 = np.dot(p1,pi.pmat)
                    nv[i]+=p2[1]
                else :
                    nv[i]+= 2.*pi.mu/(pi.length*np.pi)
	    for o in others:
		ve = o.velocities(mid[0],mid[1])
		nv[i]+=ve[1]
        return nv




    def rhs_freestream(self,vinf,others = None):
        if others != None :
	    rhs = np.zeros(len(self),dtype = float)
	    for i,p in enumerate(self):
	        mid = p.middle
		rhs[i] = -(vinf[0]*p.n[0]+vinf[1]*p.n[1])
		for o in others:
		    vo = np.dot(o.ipmat,o.velocities(mid[0],mid[1]))
		    ve = np.dot(p.pmat,vo)
		    rhs[i]+=-vo[1]
	    return rhs
	else : 
            if singul == 'Dipole':
                b = [-(vinf[0]*p.n[0]+vinf[1]*p.n[1]) for p in self]
            elif singul == 'Source':
                b = [-(vinf[0]*p.n[0]+vinf[1]*p.n[1]) for p in self]
        return b

def update_particle_speeds(part,panels = None):
    for i,oi in enumerate(part):
        oi.speed=[0.,0.]    
        for j,oj in enumerate(others):
	    if i!=j : 
	        oi.speed[0] += oj.velocities(oi[0],oi[1],mode = 'global')[0]    
	        oi.speed[1] += oj.velocities(oi[0],oi[1],mode = 'global')[1]
        if panels != None :
            for j,pj in enumerate(panels):
	        oi.speed[0] += np.dot(pj.ipmat,pj.velocities(oi[0],oi[1]))[0] 
	        oi.speed[1] += np.dot(pj.ipmat,pj.velocities(oi[0],oi[1]))[1] 
            oi.speed[0] += vinf[0]
            oi.speed[1] += vinf[1]    
 
def advection_particules(part,dt):
    for i,oi in enumerate(part):
        oi[0]+=oi.speed[0]*dt
        oi[1]+=oi.speed[1]*dt




def showAll(it):
    plt.clf()
    u , v = np.zeros_like(X),np.zeros_like(X)
    u[:] = vinf[0]
    v[:] = vinf[1]
    for p in ppp:
        p.plot(plt,axis = False)
        uloc = p.velocities(X,Y,mode='global')
        u += uloc[0]
        v += uloc[1]
    for p in others:
    #p.plot(plt,axis = False)
        plt.scatter(p[0],p[1],c = 'r')
        uloc = p.velocities(X,Y,mode='global')
        u += uloc[0]
        v += uloc[1]
    plt.streamplot(X, Y, u,v, density=4, linewidth=1, arrowsize=1, arrowstyle='->')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin,ymax)
    plt.savefig('./imgpf/img_'+str(it)+'.png')


xmin = -0.25
xmax =  3.5

Nx = 50
Ny = 50


ymin = - (0.5 * (xmax-xmin))
ymax =  (0.5 * (xmax-xmin))


eps = 0.00001
borne = 1.2
Px = np.linspace(xmin,xmax,num = Nx)
Py = np.linspace(ymin,ymax,num = Ny)
X,Y = np.meshgrid(Px,Py)
vinf = np.array([1.,0.])

geom = 'xfoil'
mode = 'unsteady'
singul = 'Dipole'
alpha_L = 01.02
dt = 0.075
Nl = 5
L = 0.4

sub.call(['rm','./imgpf/*'])

others=[]


if mode == 'stationary':
    Nit = 1
elif mode == 'unsteady':
    Nit = 2000



for it in range(Nit):
    print it
    pl = None
    if geom =='mano' :
        panels = []
        for i in range(1,Nl+1):
            p0 = PointElement(float(i-1)*L/float(Nl) ,\
                -0.15*np.cos(12.*float(it)*dt)*float(i-1)*float(i-1)*L/float(Nl)    )
            p1 = PointElement(float(i)*L/float(Nl)  ,\
                -0.15*np.cos(12.*float(it)*dt)*float(i)*float(i)*L/float(Nl)    )
            pan = DipolePanel(p0,p1)
            panels.append(pan)
        ppp = PanelArray(panels)
        pl = [ppp[-1][1][0]+vinf[0]*dt,ppp[-1][1][1]+vinf[1]*dt]
    elif geom =='xfoil' :
        incid = 0.3*np.cos(5.*float(it)*dt)
        panels = loadXfoilFile('NACAcamber0012.dat',singul,alpha = incid ,scale = 1.)
        trex = 0.5*(panels[0][0][0]+panels[-1][1][0])
        trey = 0.5*(panels[0][0][1]+panels[-1][1][1])
        trex+=vinf[0]*dt*np.cos(incid)
        trey+=vinf[1]*dt*np.sin(incid)
        pl = [trex,trey] 

    ppp = PanelArray(panels)
    M = ppp.influenceMatrix()
    b = ppp.rhs_freestream(vinf,others = others)
    stren = np.linalg.solve(M,b)
    for i,p in enumerate(ppp):
        p.mu = stren[i]
    #plt.clf()
    #plt.plot(ppp.normalVel(vinf,others))
    #plt.show()
    print ppp.normalVel(vinf,others)
    update_particle_speeds(others,panels)
    advection_particules(others,dt)
    others.append(VortexParticule(pl[0],pl[1],gamma = 0.5*(stren[-1]-stren[0])))
    showAll(it)
print '--------'
