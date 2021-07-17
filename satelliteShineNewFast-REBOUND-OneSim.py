import rebound
import numpy as np
import matplotlib.pylab as plt


""
def getXYZVVV(M, a, omega, e, Omega, inc, m0, m1):
    sim = rebound.Simulation()
    sim.add(m=m0)
    sim.add(M=M,a=a,omega=omega,e=e,Omega=Omega,m=m1,inc=inc)
    return sim.particles[1].xyz + sim.particles[1].vxyz


####################
# Written by Aaron Boley, June 2021
# Edited by Sam Lawler, June 2021
# ####################

twopi=np.pi*2
MEarth = 5.97e24
REarth = 6378.135e3
REkm = 6378.135

aukm=1.496e8
au=aukm*1e3

tilt=-23.4 * np.pi/180
#place sun on negative x axis, meaning we need to reverse the sign of the tilt
tilt*=-1

ICs=[]

ICs.append({'NPLANES':7178,'SATPP':1,'INC':30,'ALT':328})
if 1:
    ICs.append({'NPLANES':7178,'SATPP':1,'INC':40,'ALT':334})
    ICs.append({'NPLANES':7178,'SATPP':1,'INC':53,'ALT':345})
    ICs.append({'NPLANES':40,'SATPP':50,'INC':96.9,'ALT':360})
    ICs.append({'NPLANES':1998,'SATPP':1,'INC':75,'ALT':373})
    ICs.append({'NPLANES':4000,'SATPP':1,'INC':53,'ALT':499})
    ICs.append({'NPLANES':12,'SATPP':12,'INC':148,'ALT':604})
    ICs.append({'NPLANES':18,'SATPP':18,'INC':115.7,'ALT':614})

    ICs.append({'NPLANES':2547,'SATPP':1,'INC':53,'ALT':345.6})
    ICs.append({'NPLANES':2478,'SATPP':1,'INC':48,'ALT':340.8})
    ICs.append({'NPLANES':2493,'SATPP':1,'INC':42,'ALT':335.9})
    ICs.append({'NPLANES':32,'SATPP':50,'INC':53,'ALT':550})
    ICs.append({'NPLANES':72,'SATPP':22,'INC':53.2,'ALT':540})
    ICs.append({'NPLANES':36,'SATPP':20,'INC':70,'ALT':570})
    ICs.append({'NPLANES':6,'SATPP':58,'INC':97.6,'ALT':560})
    ICs.append({'NPLANES':4,'SATPP':43,'INC':97.6,'ALT':560.1})
    ICs.append({'NPLANES':18,'SATPP':40,'INC':87.9,'ALT':1200})
    ICs.append({'NPLANES':36,'SATPP':49,'INC':87.9,'ALT':1201})
    ICs.append({'NPLANES':32,'SATPP':72,'INC':40,'ALT':1202})
    ICs.append({'NPLANES':32,'SATPP':72,'INC':55,'ALT':1203})

    ICs.append({'NPLANES':16,'SATPP':30,'INC':85,'ALT':590})
    ICs.append({'NPLANES':40,'SATPP':50,'INC':50,'ALT':600})
    ICs.append({'NPLANES':60,'SATPP':60,'INC':55,'ALT':508})
    ICs.append({'NPLANES':48,'SATPP':36,'INC':30,'ALT':1145})
    ICs.append({'NPLANES':48,'SATPP':36,'INC':40,'ALT':1145})
    ICs.append({'NPLANES':48,'SATPP':36,'INC':50,'ALT':1145})
    ICs.append({'NPLANES':48,'SATPP':36,'INC':60,'ALT':1145})

    ICs.append({'NPLANES':34,'SATPP':34,'INC':51.9,'ALT':630})
    ICs.append({'NPLANES':36,'SATPP':36,'INC':42,'ALT':610})
    ICs.append({'NPLANES':28,'SATPP':28,'INC':33,'ALT':509})

xobs0 = REarth*1
yobs0 = 0.
zobs0 = 0.

def RotateZ(x,y,z,a):
    xp = x*np.cos(a)-y*np.sin(a)
    yp = x*np.sin(a)+y*np.cos(a)
    return xp,yp,z

def RotateY(x,y,z,a):
    xp=x*np.cos(a)+z*np.sin(a)
    zp=-x*np.sin(a)+z*np.cos(a)
    return xp,y,zp


sim = rebound.Simulation()
sim.G = 6.67430e-11
sim.add(m=MEarth)

for ic in ICs:
    nplanes=ic['NPLANES']
    nsat=ic['SATPP']
    ntot=nplanes*nsat
    sma = (ic['ALT']+REkm)*1000.

    dplane=twopi/nplanes
    mashuf = np.linspace(0,twopi,nsat,endpoint=False)
    for ilocal in range(nsat):
        for iplane in range(nplanes):
            Omega = iplane*dplane
            sim.add(M="uniform",a=sma,omega=0,e=0,Omega=Omega,inc=ic['INC']*np.pi/180.)


""
sim.getWidget(orbits=False,size=(1300,1300),pointsize=6)

""
sim.dt = sim.particles[1].P/15000.
sim.integrator = "whfast"
sim.N_active = 1

""
sim.t = 0
sim.integrate(10000.)

""

