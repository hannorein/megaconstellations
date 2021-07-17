import rebound
import numpy as np
import numpy.ma as ma
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
# Edited by Hanno Rein, July 2021
# ####################

twopi=np.pi*2
MEarth = 5.97e24
REarth = 6378.135e3
A=4. # m^2   -> effective cross section of satellite
albedo=0.2 # -> effective albedo 
EarthTilt=23.4 * np.pi/180

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

sim = rebound.Simulation()
sim.G = 6.67430e-11
sim.add(m=MEarth)

for ic in ICs:
    nplanes=ic['NPLANES']
    nsat=ic['SATPP']
    sma = ic['ALT']*1000.+REarth

    dplane=twopi/nplanes
    for ilocal in range(nsat):
        for iplane in range(nplanes):
            Omega = iplane*dplane
            sim.add(M="uniform",a=sma,omega=0,e=0,Omega=Omega,inc=ic['INC']*np.pi/180.)


""
sim.dt = sim.particles[1].P/15000.
sim.integrator = "whfast"
sim.N_active = 1


""
def rotY(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,0,-s],[0,1,0],[s,0,c]])
    return xyz @ M
def rotZ(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    return xyz @ M


""
timeOfYear = 0. # in rad
def getStereographic(latitude, tilt):
    sun = np.array([-1.4959787e+11,0,0]) # in m
    earth = np.array([0,0,0])
    obs_n = np.array([np.cos(latitude),0,np.sin(latitude)])
    obs = obs_n * REarth    
    
    xyz = np.zeros((sim.N,3),dtype="float64")
    sim.serialize_particle_data(xyz=xyz)
    xyz = xyz[1:] # remove earth    
    xyz = rotY(xyz, EarthTilt)
    xyz = rotZ(xyz, timeOfYear)
    
    lit = np.logical_or(xyz[:,0]<0,np.linalg.norm(xyz[:,1:3],axis=1)>REarth)
    xyz = xyz[lit] 
    
    xyz = xyz - obs
    xyz_d = np.linalg.norm(xyz,axis=1)
    xyz_n = xyz/xyz_d[:,np.newaxis]
    
    phase = np.arccos(np.clip(xyz_n[:,0], -1.0, 1.0)) # assume sun is in -x direction
    
    fac1 = 2/(3*np.pi**2)
    magV = -26.74 -2.5*np.log10(fac1 * A * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) ) ) + 5 * np.log10(xyz_d)

    
    xyz_n = rotY(xyz_n,latitude)
    
    above_horizon = xyz_n[:,0]>0
    xyz_n = xyz_n[above_horizon]
    magV = magV[above_horizon]
    
    return xyz_n[:,1:3]/(1.+xyz_n[:,0,np.newaxis]), magV

""
latitudes = [90,40,0]
hours = [0,10]
magVmin, magVmax = 4, 8
fig, axs = plt.subplots(len(latitudes),len(hours),squeeze=False, figsize=(4*len(hours),4*len(latitudes)))
cm = plt.cm.get_cmap('viridis')
for i, latitude in enumerate(latitudes):
    axs[i,0].set_ylabel("Latitude: %.0f"%latitude)
    for j, hour in enumerate(hours):
        axs[0,j].set_title("Hour: %.0f"%hour)
        ax = axs[i,j]
        ax.set_aspect("equal")
        plot_size = 1.14
        ax.set_xlim(-plot_size,plot_size)
        ax.set_ylim(-plot_size,plot_size)
        
        xyzf_stereographic, magV = getStereographic(latitude/180.*np.pi, tilt/180.*np.pi)
        im=ax.scatter(xyzf_stereographic[:,0],xyzf_stereographic[:,1],s=2, c=magV,cmap=cm,vmin=magVmin,vmax=magVmax)
        
        card_pos = 1.07
        ax.text(0,card_pos, "N",ha="center",va="center",weight="bold")
        ax.text(0,-card_pos, "S",ha="center",va="center",weight="bold")
        ax.text(card_pos,0, "W",ha="center",va="center",weight="bold")
        ax.text(-card_pos,0, "E",ha="center",va="center",weight="bold")
        circle1 = plt.Circle((0, 0), 1, fill=False)
        ax.add_patch(circle1)
        ax.text(-1,1,"N=%d"%len(xyzf_stereographic))
#fig.tight_layout()
fig.colorbar(im,ax=axs.ravel().tolist(),label="magV");


""
np.logical_or(xyzr[:,0]<0,np.linalg.norm(xyzr[:,1:3],axis=1)>REarth)

""
np.linalg.norm(xyzr[:,1:3],axis=1)>REarth

""
np.linalg.norm(xyzr[:,1:3],axis=1)

""
obs_n = np.array([np.cos(latitude),0,np.sin(latitude)])
obs = obs_n * REarth

""
xyz.shape
xyz_n = xyz/np.linalg.norm(xyz,axis=1)[:,np.newaxis]

""
angle = np.arccos(np.clip(np.dot(xyz_n, obs_n), -1.0, 1.0))

""
angle[3212]

""

