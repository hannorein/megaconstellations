import os
os.environ[ 'MPLCONFIGDIR' ] = '/home/rein/megaconstellations/webapp/'
os.chdir("/home/rein/megaconstellations/webapp/")

from flask import Flask, render_template, request
import rebound
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from matplotlib.collections import PatchCollection
import matplotlib.transforms as transforms
import mega
import io
import random
from flask import Response
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

cm = plt.cm.get_cmap('plasma_r')
sims = mega.get_simulations()

app = Flask(__name__,static_folder='static')
app.config['TEMPLATES_AUTO_RELOAD'] = True

latitudes = {"North pole":90., "Canada":50., "Hawaii":20., "Equator": 0., "Chile":-30., "South pole":-90.};

@app.route("/", methods=['post','get'])
def index():
    latitude = 50.
    month = 3
    hour = 0.
    albedo = 0.2
    area = 4.0
    enabled_constellations = [k for k in mega.constellations]
    return render_template('index.html', latitudes=latitudes, latitude=latitude, month=month, area=area, albedo=albedo, hour=hour, constellations=sims, enabled_constellations=enabled_constellations)


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

@app.route('/plot.png')
def plot_png():
    try:
        latitude = float(request.args.get('latitude'))
    except:
        latitude = 50.
    try:
        month = float(request.args.get('month')) 
    except:
        month = 3.
    try:
        hour = float(request.args.get('hour'))
    except:
        hour = 0.
    try:
        enabled_constellations = request.args.getlist('constellation')
    except:
        enabled_constellations = []
    try:
        albedo = float(request.args.get('albedo'))
    except:
        albedo = 0.2
    try:
        area = float(request.args.get('area'))
    except:
        area = 4.

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


    N = 0
    esims = {}
    for k in enabled_constellations:
        esims[k] = sims[k]

    xyzf_stereographic, magV = mega.get_stereographic_data(esims, latitude=latitude, month=month, hour=hour, albedo=albedo, area=area)
    if xyzf_stereographic is not None:
        N = len(xyzf_stereographic)
        im=ax.scatter(xyzf_stereographic[:,0],xyzf_stereographic[:,1],s=4, c=magV,cmap=cm,vmin=magVmin,vmax=magVmax)

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


