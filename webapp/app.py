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

sims = {}
for constellation in constellations.names:
    print(constellation)
    sim = rebound.Simulation()
    sim.G = 6.67430e-11
    sim.add(m=MEarth)
    constellations.add_to_sim(sim,constellations.satellites[constellation], debug=True)
    sims[constellation] = sim

app = Flask(__name__)

@app.route("/", methods=['post','get'])
def hello_world():
    # Default:
    location = 50.
    timeofyear = 0.
    print(np.pi)
    timeofday = 0.
    enabled_constellations = [k for k in constellations.names]

    if request.method == 'POST':
        location = float(request.form.get('location'))  # access the data inside
        timeofyear = float(request.form.get('timeofyear'))  # access the data inside
        timeofday = float(request.form.get('timeofday'))  # access the data inside
        enabled_constellations = request.form.getlist('constellation')
    return render_template('index.html', location=location, timeofyear=timeofyear, timeofday=timeofday, constellations=constellations.names, enabled_constellations=enabled_constellations)


def create_figure():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    xs = range(100)
    ys = [random.randint(1, 50) for x in xs]
    axis.plot(xs, ys)
    return fig

@app.route('/plot.png')
def plot_png():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')
