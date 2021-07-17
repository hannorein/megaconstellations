import rebound
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as plt
from matplotlib.collections import PatchCollection
from matplotlib.animation import FFMpegWriter


""
def getXYZVVV(M, a, omega, e, Omega, inc, m0, m1):
    sim = rebound.Simulation()
    sim.add(m=m0)
    sim.add(M=M,a=a,omega=omega,e=e,Omega=Omega,m=m1,inc=inc)
    return sim.particles[1].xyz + sim.particles[1].vxyz

""
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
def getStereographic(latitude, tilt, hour):
    sun = np.array([-1.4959787e+11,0,0]) # in m
    sun = rotY(sun, tilt)
    sun_n = sun/np.linalg.norm(sun)
    
    obs = np.array([REarth, 0, 0])
    obs = rotY(obs, -latitude)
    obs = rotZ(obs, hour)
    obs_n = obs/np.linalg.norm(obs)
    
    xyz = np.zeros((sim.N,3),dtype="float64")
    sim.serialize_particle_data(xyz=xyz)
    xyz = xyz[1:] # remove earth    

    lit = np.linalg.norm(np.cross(xyz,sun_n),axis=1)>REarth
    
    xyz = xyz[lit]
    
    xyz_n = xyz/np.linalg.norm(xyz,axis=1)[:,np.newaxis]
    xyz_r = xyz - obs
    xyz_rd = np.linalg.norm(xyz_r,axis=1)
    xyz_rn = xyz_r/xyz_rd[:,np.newaxis]
    
    phase = np.arccos(np.clip(np.dot(xyz_rn, sun_n), -1.0, 1.0)) # assume sun is in -x direction
    
    fac1 = 2/(3*np.pi**2)
    magV = -26.74 -2.5*np.log10(fac1 * A * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) ) ) + 5 * np.log10(xyz_rd)

        
    elevation = np.pi/2.-np.arccos(np.dot(xyz_rn,obs_n))
    
    xyz = rotZ(xyz, -hour)
    xyz = rotY(xyz, latitude)
    xyz_r = xyz - np.array([REarth, 0, 0])
    xyz_rd = np.linalg.norm(xyz_r,axis=1)
    xyz_rn = xyz_r/xyz_rd[:,np.newaxis]

    xyz_rn = xyz_rn[elevation>0.]
    magV = magV[elevation>0.]


    return xyz_rn[:,1:3]/(1.+xyz_rn[:,0,np.newaxis]), magV

""
latitudes = [-30,20,-50]
hours = np.linspace(-7,7,1000,endpoint=True)
sim.t = 0
for k, hour in enumerate(hours):
    if k>0:
        sim.dt = (hours[k-1]-hours[k])*60*60
        sim.step()
    timesOfYear = {"winter solstice":-np.pi/2.,"equinox":0.,"summer solstice":np.pi/2.}
    magVmin, magVmax =5, 8
    fig, axs = plt.subplots(len(latitudes),len(timesOfYear),squeeze=False, figsize=(2+4*len(timesOfYear),4*len(latitudes)))
    cm = plt.cm.get_cmap('viridis_r')
    for i, latitude in enumerate(latitudes):
        axs[i,0].set_ylabel("Latitude: %.0f"%latitude)
        for j, timeOfYear in enumerate(timesOfYear):
            axs[0,j].set_title(timeOfYear)
            ax = axs[i,j]
            ax.set_aspect("equal")
            plot_size = 1.14
            ax.set_xlim(-plot_size,plot_size)
            ax.set_ylim(-plot_size,plot_size)
            ax.get_xaxis().set_ticks([])
            ax.get_yaxis().set_ticks([])
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)

            circle = plt.Circle((0, 0), 1)
            coll = PatchCollection([circle], zorder=-10, color="black")
            ax.add_collection(coll)


            tilt = 23.4*np.sin(timesOfYear[timeOfYear])

            xyzf_stereographic, magV = getStereographic(latitude/180.*np.pi, tilt/180.*np.pi, hour/12.*np.pi)
            im=ax.scatter(xyzf_stereographic[:,0],xyzf_stereographic[:,1],s=4, c=magV,cmap=cm,vmin=magVmin,vmax=magVmax)

            card_pos = 1.07
            ax.text(0,card_pos, "N",ha="center",va="center",weight="bold")
            ax.text(0,-card_pos, "S",ha="center",va="center",weight="bold")
            ax.text(card_pos,0, "W",ha="center",va="center",weight="bold")
            ax.text(-card_pos,0, "E",ha="center",va="center",weight="bold")
            ax.text(-1,0.85,"N=%d"%len(xyzf_stereographic))
    fig.suptitle("Midnight %+.1fh"%hour, fontsize=16)
    fig.tight_layout()
    fig.colorbar(im,ax=axs.ravel().tolist(),label="magV");
    fig.savefig("plot_%05d.png"%k, facecolor="white")
    plt.close(fig)

""

