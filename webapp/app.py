from flask import Flask, render_template, request
import rebound
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from matplotlib.collections import PatchCollection
import matplotlib.transforms as transforms
import constellations
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

twopi=np.pi*2
MEarth = 5.97e24
REarth = 6378.135e3
A=4. # m^2   -> effective cross section of satellite
albedo=0.2 # -> effective albedo
EarthTilt=23.4 * np.pi/180
cm = plt.cm.get_cmap('plasma_r')

sims = {}
for constellation in constellations.names:
    print(constellation)
    sim = rebound.Simulation()
    sim.G = 6.67430e-11
    sim.add(m=MEarth)
    constellations.add_to_sim(sim,constellations.satellites[constellation], debug=False)
    sims[constellation] = sim

app = Flask(__name__)


@app.route("/", methods=['post','get'])
def hello_world():
    latitude = 50.
    timeofyear = 3
    timeofday = 0.
    enabled_constellations = [k for k in constellations.names]
    return render_template('index.html', latitude=latitude, timeofyear=timeofyear, timeofday=timeofday, constellations=constellations.names, enabled_constellations=enabled_constellations)


def rotY(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,0,-s],[0,1,0],[s,0,c]])
    return xyz @ M
def rotZ(xyz,alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    M = np.array([[c,-s,0],[s,c,0],[0,0,1]])
    return xyz @ M
def lengthOfNight(timeOfYear,latitude, p=0):
    # https://www.ikhebeenvraag.be/mediastorage/FSDocument/171/Forsythe+-+A+model+comparison+for+daylength+as+a+function+of+latitude+and+day+of+year+-+1995.pdf
    # p=18 for astronomical twilight
    theta = 2.*np.arctan(0.9671396*np.tan(-timeOfYear/2.+np.pi/4.))
    phi = np.arcsin(0.39795*np.cos(theta))
    return 24./np.pi * np.arccos((np.sin(p*np.pi/180.)+np.sin(latitude*np.pi/180.)*np.sin(phi))/(np.cos(latitude*np.pi/180.)*np.cos(phi)))

def getStereographic(sim, latitude, tilt, hour):
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
    magV = -26.74 -2.5*np.log10(fac1 * A * albedo * ( (np.pi-phase)*np.cos(phase) + np.sin(phase) ) ) + 5 * np.log10(xyz_rd)


    elevation = (np.pi/2.-np.arccos(np.dot(xyz_rn,obs_n)))/np.pi*180.

    xyz = rotZ(xyz, -hour)
    xyz = rotY(xyz, latitude)
    xyz_r = xyz - np.array([REarth, 0, 0])
    xyz_rd = np.linalg.norm(xyz_r,axis=1)
    xyz_rn = xyz_r/xyz_rd[:,np.newaxis]

    elevation_cut = 0
    xyz_rn = xyz_rn[elevation>elevation_cut]
    #magV = phase[elevation>elevation_cut]/np.pi*180.
    magV = magV[elevation>elevation_cut]


    return xyz_rn[:,1:3]/(1.+xyz_rn[:,0,np.newaxis]), magV

@app.route('/plot.png')
def plot_png():
    try:
        latitude = float(request.args.get('latitude'))
    except:
        latitude = 50.
    try:
        timeofyear = float(request.args.get('timeofyear')) 
    except:
        timeofyear = 3.
    try:
        timeofday = float(request.args.get('timeofday'))
    except:
        timeofday = 0.
    try:
        enabled_constellations = request.args.getlist('constellation')
    except:
        enabled_constellations = []
    if len(enabled_constellations)==0: 
        enabled_constellations = [k for k in constellations.names]

    magVmin, magVmax =5, 8
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.set_aspect("equal")
    plot_size = 1.1
    ax.set_xlim(-plot_size,plot_size)
    ax.set_ylim(-plot_size,plot_size)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Background
    circle = plt.Circle((0, 0), 1)
    coll = PatchCollection([circle], zorder=-10, color="black")
    ax.add_collection(coll)
    
    # 30 degree 
    elevation = 30
    circle = plt.Circle((0, 0), np.cos(elevation/180.*np.pi)/(1.+np.sin(elevation/180.*np.pi)))
    coll = PatchCollection([circle], zorder=3,edgecolor="dimgray",facecolor="none")
    ax.add_collection(coll)
    elevation = 0.
    circle = plt.Circle((0, 0), np.cos(elevation/180.*np.pi)/(1.+np.sin(elevation/180.*np.pi)))
    coll = PatchCollection([circle], zorder=3,edgecolor="k",facecolor="none",lw=3)
    ax.add_collection(coll)


    tilt = 23.4*np.sin(timeofyear/6.*np.pi)

    N = 0
    for k in enabled_constellations:
        xyzf_stereographic, magV = getStereographic(sims[k], latitude/180.*np.pi, tilt/180.*np.pi, timeofday/12.*np.pi)
        im=ax.scatter(xyzf_stereographic[:,0],xyzf_stereographic[:,1],s=4, c=magV,cmap=cm,vmin=magVmin,vmax=magVmax)
        N += len(xyzf_stereographic);

    #im=ax.scatter(xyzf_stereographic[:,0],xyzf_stereographic[:,1],s=4,color="r")
    card_pos = 1.11
    ax.text(0,card_pos, "N",ha="center",va="center")
    ax.text(0,-card_pos, "S",ha="center",va="center")
    ax.text(card_pos,0, "W",ha="center",va="center")
    ax.text(-card_pos,0, "E",ha="center",va="center")
    ax.text(0.7,-0.95,"N=%d"%N,fontsize=12)

    fig.tight_layout()

    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')


