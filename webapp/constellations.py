import numpy as np
MEarth = 5.97e24
REarth = 6378.135e3

constellations = {
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
        {'NPLANES':36,'SATPP':49,'INC':87.9,'ALT':1201},
        {'NPLANES':32,'SATPP':72,'INC':40,'ALT':1202},
        {'NPLANES':32,'SATPP':72,'INC':55,'ALT':1203},],
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


def add_to_sim(sim, ICs, debug=True):
     for IC in ICs:
        nplanes=IC['NPLANES']
        nsat=IC['SATPP']
        a = IC['ALT']*1000.+REarth

        Omegas = np.linspace(0.,2.*np.pi,nplanes)
        for i, Omega in enumerate(Omegas):
            # 5 percent jitter
            Ms = np.linspace(0.,2.*np.pi,nsat)+ 2.*np.pi/nsat*0.05*np.random.normal(size=nsat)
            for j, M in enumerate(Ms):
                sim.add(M=M, a=a, omega=0, e=0, Omega=Omega, inc=IC['INC']*np.pi/180.)
                if debug and sim.N>100:
                    return
