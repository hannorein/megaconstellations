import numpy as np
import rebound
MEarth = 5.97e24
REarth = 6378.135e3

constellations_all = {
    "Starlink": [ {'NPLANES':7178,'SATPP':1,'INC':30,'ALT':328},
        {'NPLANES':7178,'SATPP':1,'INC':40,'ALT':334},
        {'NPLANES':7178,'SATPP':1,'INC':53,'ALT':345},
        {'NPLANES':40,'SATPP':50,'INC':96.9,'ALT':360},
        {'NPLANES':1998,'SATPP':1,'INC':75,'ALT':373},
        {'NPLANES':4000,'SATPP':1,'INC':53,'ALT':499},
        {'NPLANES':12,'SATPP':12,'INC':148,'ALT':604},
        {'NPLANES':18,'SATPP':18,'INC':115.7,'ALT':614},
        {'NPLANES':2547,'SATPP':1,'INC':53,'ALT':345.6},
        {'NPLANES':2478,'SATPP':1,'INC':48,'ALT':340.8},
        {'NPLANES':2493,'SATPP':1,'INC':42,'ALT':335.9},
        {'NPLANES':32,'SATPP':50,'INC':53,'ALT':550},
        {'NPLANES':72,'SATPP':22,'INC':53.2,'ALT':540},
        {'NPLANES':36,'SATPP':20,'INC':70,'ALT':570},
        {'NPLANES':6,'SATPP':58,'INC':97.6,'ALT':560},
        {'NPLANES':4,'SATPP':43,'INC':97.6,'ALT':560.1},],
    "OneWeb": [ {'NPLANES':18,'SATPP':40,'INC':87.9,'ALT':1200},
        {'NPLANES':36,'SATPP':49,'INC':87.9,'ALT':1200},
        {'NPLANES':32,'SATPP':72,'INC':40,'ALT':1200},
        {'NPLANES':32,'SATPP':72,'INC':55,'ALT':1200},],
    "StarNet/GW": [ {'NPLANES':16,'SATPP':30,'INC':85,'ALT':590},
        {'NPLANES':40,'SATPP':50,'INC':50,'ALT':600},
        {'NPLANES':60,'SATPP':60,'INC':55,'ALT':508},
        {'NPLANES':48,'SATPP':36,'INC':30,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':40,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':50,'ALT':1145},
        {'NPLANES':48,'SATPP':36,'INC':60,'ALT':1145},],
    "Kuiper": [ {'NPLANES':34,'SATPP':34,'INC':51.9,'ALT':630},
        {'NPLANES':36,'SATPP':36,'INC':42,'ALT':610},
        {'NPLANES':28,'SATPP':28,'INC':33,'ALT':509},],
    }

def getAirmass(z):
    # z is the Zenith enagle
    X = 1./(np.cos(z) + 0.50572*(6.07995+90-z*180/np.pi)**(-1.6364))  # Kasten and Young (1989)
    #X = 1./np.cos(z) * (1-0.0012*np.tan(z)**2)  # Young and Irvine (1967)
    return X

def add_to_simulation(sim, ICs, debug=False):
    for IC in ICs:
        nplanes=IC['NPLANES']
        nsat=IC['SATPP']
        a = IC['ALT']*1000.+REarth

        Omegas = np.linspace(0.,2.*np.pi,nplanes)
        for i, Omega in enumerate(Omegas):
            # 5 percent jitter
            Ms = np.linspace(0.,2.*np.pi,nsat)+ 2.*np.pi/nsat*0.25*np.random.normal(size=nsat)
            for j, M in enumerate(Ms):
                sim.add(M=M, a=a, omega=0, e=0, Omega=Omega, inc=IC['INC']*np.pi/180.)
                if debug and sim.N>100:
                    return
def rotY(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,0,-s],[0,1,0],[s,0,c]])
    return xyz @ M
def rotZ(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    return xyz @ M

def length_of_night(month,latitude, p=0):
    # https://www.ikhebeenvraag.be/mediastorage/FSDocument/171/Forsythe+-+A+model+comparison+for+daylength+as+a+function+of+latitude+and+day+of+year+-+1995.pdf
    # p=18 for astronomical twilight
    day = month/12*365.25+79
    theta = 0.2163108+2.*np.arctan(0.9671396*np.tan(0.00860*(day-186)))
    phi = np.arcsin(0.39795*np.cos(theta))
    arccosarg = (np.sin(p*np.pi/180.)+np.sin(latitude/180.*np.pi)*np.sin(phi))/(np.cos(latitude/180.*np.pi)*np.cos(phi))
    if abs(arccosarg)>=1.:
        return 0.0
    return 24./np.pi * np.arccos(arccosarg)

def get_stereographic_data(sims, latitude=0., month=0., hour=0., albedo=0.2, area=4., airmassCoeff=0.2, randomCoeff=0.5, elevation_cut = 0):
    # latitude in degrees
    # month in months from spring euquinox
    # hours in hours since midnight
    latitude = latitude/180.*np.pi 
    tilt = 23.4*np.sin(month/6.*np.pi)/180.*np.pi
    hour = hour/12.*np.pi
    xy, mag = [], []     
    for name in sims:
        sim = sims[name]
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

        phase = np.arccos(np.clip(np.dot(xyz_rn, -sun_n), -1.0, 1.0)) # assume sun is in -x direction

        fac1 = 2/(3*np.pi**2)
        pfac = 3.1
        m_sun = -26.47 # g' band 
        magV = m_sun -2.5*np.log10(fac1 * area * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) ) ) + 5 * np.log10(xyz_rd)
        #magV = m_sun -2.5*np.log10(2/(3*np.pi**(pfac+1)) * area * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) )**pfac ) + 5 * np.log10(xyz_rd)


        elevation = (np.pi/2.-np.arccos(np.dot(xyz_rn,obs_n)))/np.pi*180.

        xyz = rotZ(xyz, -hour)
        xyz = rotY(xyz, latitude)
        xyz_r = xyz - np.array([REarth, 0, 0])
        xyz_rd = np.linalg.norm(xyz_r,axis=1)
        xyz_rn = xyz_r/xyz_rd[:,np.newaxis]

        #elevation_cut = 45
        xyz_rn = xyz_rn[elevation>elevation_cut]
        magV = magV[elevation>elevation_cut]

        airmass = getAirmass((90.-elevation[elevation>elevation_cut])*np.pi/180.)
        magV += airmassCoeff*airmass
        if randomCoeff>0.:
            magV += randomCoeff*np.random.normal(0.,1.,size=len(magV))

        xy.append(xyz_rn[:,1:3]/(1.+xyz_rn[:,0,np.newaxis]))
        mag.append(magV)
    if len(xy)>0:
        return np.concatenate(xy), np.concatenate(mag) 
    else:
        return None, None


def get_simulations(constellations=None, use_cache=True):
    if constellations is None:
        constellations = constellations_all
    sims = {}
    for c in constellations.keys():
        sim = None
        if use_cache:
            filename = "mega_"+"".join(x for x in c if x.isalnum())+".bin"
            try:
                sim = rebound.Simulation(filename)
            except:
                # need to create simulation
                pass
        if sim is None:
            sim = rebound.Simulation()
            sim.G = 6.67430e-11
            sim.add(m=MEarth)
            sim.N_active = 1
            add_to_simulation(sim, constellations[c])
            if use_cache:
                sim.save(filename)
        sims[c] = sim
    return sims
