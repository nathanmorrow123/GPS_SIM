#Nathan Morrow AFIT Fall 2023
#All GPS data is downloaded from NAVCE
#https://www.navcen.uscg.gov/gps-nanus-almanacs-opsadvisories-sof

from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from collections import namedtuple
from functools import partial
import itertools
import numpy as np
import time
import math
import yuma_decode
RADIUS_EARTH=6.371e6 #meters
ROTATION_EARTH=7.2921150e-5 #rad/s
GRAV_CONSTANT=6.67430e-11 #N*m^2/kg^2
EARTH_MASS=5.97219e24 #kg
ANTENNA_LOCATION = [39.765659548785884, -84.19151722210734] #Latitude/Longitude in degrees
#Assuming antenna is normal to earth surface

def getSatPos(args, t):
	"""
    Calculate the position of a satellite in Cartesian coordinates.

    Parameters:
    args (list): List containing satellite parameters.

    Returns:
    list: Satellite ID and its x, y, z coordinates.
    """
	[satID,health,e,toa,inc,ascRate,sqrt_a,ra_week,w,E,Af0,Af1] = args
	t = t + toa # Sync Animation to TOA
	a = sqrt_a**2
	b = math.sqrt(a**2*(1-e**2)) #https://jtauber.github.io/orbits/019.html
	n = math.sqrt(GRAV_CONSTANT * EARTH_MASS / a**3)  # Mean motion
	M0 = E - e * math.sin(E)  # Mean anomaly at reference time
	M = M0 + n * (t - toa)  # Mean anomaly at time t

    # Solve Kepler's equation for Eccentric Anomaly (E)
	E_new = M
	for _ in range(10):  # Iterate to solve for E
		E_new = M + e * math.sin(E_new)

    # True anomaly
	v = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E_new / 2), math.sqrt(1 - e) * math.cos(E_new / 2))
	
	#position in the orbit plane
	p = a*(math.cos(E)-e) #https://en.wikipedia.org/wiki/Kepler%27s_equation
	q = b*math.sin(E)
	#rotate by argument of periapsis
	x = math.cos(w) * p - math.sin(w) * q
	y = math.sin(w) * p + math.cos(w) * q
	#rotate by inclination
	z = math.sin(inc) * y
	y = math.cos(inc) * y
	#rotate by longitude of ascending node
	xtemp = x
	x = math.cos(ra_week) * xtemp - math.sin(ra_week) * y
	y = math.sin(ra_week) * xtemp + math.cos(ra_week) * y
	return [satID,x,y,z]
def meshTolines(a,b,c):
	"""
    Convert mesh grid data to line segments for 3D plotting.

    Parameters:
    a (ndarray): X coordinates of the mesh grid.
    b (ndarray): Y coordinates of the mesh grid.
    c (ndarray): Z coordinates of the mesh grid.

    Returns:
    tuple: Arrays of start and end points for line segments.
    """
	m,n,o,p,q,r = ([],[],[],[],[],[])
	for long , elem in enumerate(a):
		for lat, elem in enumerate(a[long]):
			if(long<len(a)-1 and lat <len(a[long])-1):
				m.append(a[long,lat])
				p.append(a[long+1,lat])
				m.append(a[long,lat])
				p.append(a[long,lat+1])
	for long , elem in enumerate(b):
		for lat, elem in enumerate(b[long]):
			if(long<len(b)-1 and lat <len(b[long])-1):
				n.append(b[long,lat])
				q.append(b[long+1,lat])
				n.append(b[long,lat])
				q.append(b[long,lat+1])
	for long , elem in enumerate(c):
		for lat, elem in enumerate(c[long]):
			if(long<len(c)-1 and lat <len(c[long])-1):
				o.append(c[long,lat])
				r.append(c[long+1,lat])
				o.append(c[long,lat])
				r.append(c[long,lat+1])
	return (np.stack(m),np.stack(n),np.stack(o),np.stack(p),np.stack(q),np.stack(r))
def plotEarth(ax,tof,mps):
	"""
    Plot the Earth rotation as line segments for the given time of flight.

    Parameters:
    ax (Axes3D): Matplotlib 3D axes object.
    tof (float): Time of flight in seconds.
    mps (float): Number of simulated location points per second.

    Returns:
    ndarray: Array of line segments for each time step.
    """
	times = range(int(tof*mps))
	earthMovement = []
	earth_segs = []
	for t in times:
		c, s = np.cos((t/mps)*ROTATION_EARTH), np.sin((t/mps)*ROTATION_EARTH)
		R = np.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		x = RADIUS_EARTH*np.cos(u)*np.sin(v)
		y = RADIUS_EARTH*np.sin(u)*np.sin(v)
		z = RADIUS_EARTH*np.cos(v)
		preImage = [x,y,z]
		postImage = [x*R[0,0]+y*R[0,1],x*R[1,0]+y*R[1,1],z]
		earth_segs.append(meshTolines(postImage[0],postImage[1],postImage[2]))
	return np.stack(earth_segs)
def plotAntenna(ax,tof,mps,ant):
	"""
    Plot the antenna movement as it rotates with the Earth.

    Parameters:
    ax (Axes3D): Matplotlib 3D axes object.
    tof (float): Time of flight in seconds.
    mps (float): Number of simulated location points per second.
    ant (list): Antenna position in spherical coordinates (radius, theta, phi).

    Returns:
    ndarray: Array of antenna positions for each time step.
    """
	times = range(int(tof*mps))
	antennaMovement = []
	antCart = sphericalTocartesian(ant)
	for t in times:
		ant[2]=-1*ROTATION_EARTH*t/mps
		antCart = sphericalTocartesian(ant)
		antennaMovement.append([antCart[0],antCart[1],antCart[2]])
	return np.stack(antennaMovement)
def plotSats(ax,satCoords):
	"""
    Plot the satellite positions on the 3D plot.

    Parameters:
    ax (Axes3D): Matplotlib 3D axes object.
    satCoords (ndarray): Array of satellite coordinates.

    Returns:
    None
    """
	for satCoord in satCoords:
		ax.scatter(satCoord[1],satCoord[2],satCoord[3],c='green',marker='.',s=20)
		ax.text(satCoord[1],satCoord[2],satCoord[3],int(satCoord[0]))

def sphericalTocartesian(sCoords): #sCoords (rad,theta,phi)
	"""
    Convert spherical coordinates to Cartesian coordinates.

    Parameters:
    sCoords (list): Spherical coordinates [radius, theta, phi].

    Returns:
    list: Cartesian coordinates [x, y, z].
	"""
	cartX=sCoords[0]*math.sin(sCoords[1])*math.cos(sCoords[2])
	cartY=sCoords[0]*math.sin(sCoords[1])*math.sin(sCoords[2])
	cartZ=sCoords[0]*math.cos(sCoords[1])
	return [cartX,cartY,cartZ]
def init_fig(fig,ax,artists):
	"""
    Initialize the 3D plot with labels and settings.

    Parameters:
    fig (Figure): Matplotlib figure object.
    ax (Axes3D): Matplotlib 3D axes object.
    artists (Artists): Named tuple containing the wireframe and quiver artists.

    Returns:
    Artists: The updated artists.
    """
	ax.set_xlabel("X axis")
	ax.set_ylabel("Y Axis")
	ax.set_zlabel("Z Axis")
	#ax.grid(False)
	ax.xaxis.set_pane_color((.1, .1, .1, .1))
	ax.yaxis.set_pane_color((.1, .1, .1, .1))
	ax.zaxis.set_pane_color((.1, .1, .1, .1))
	ax.set_xlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
	ax.set_ylim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
	ax.set_zlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
	plotSats(ax,satCoords[1:-1]) # Skip the first one it's bad data
	return artists

def update_artists(frames,artists):
	"""
    Update the artists for each animation frame.

    Parameters:
    frames (tuple): Frame data containing segments for wireframe and quiver plots.
    artists (Artists): Named tuple containing the wireframe and quiver artists.

    Returns:
    Artists: The updated artists.
    """
	x,y,z,u,v,w,m,n,o,p,q,r = frames
	temp = np.array([x,y,z,u,v,w]).reshape(6,-1)
	qSegs = [[[x,y,z],[u,v,w]]for x,y,z,u,v,w in zip(*temp.tolist())]
	temp = np.array([m,n,o,p,q,r]).reshape(6,-1)
	wSegs = [[[m,n,o],[p,q,r]]for m,n,o,p,q,r in zip(*temp.tolist())]
	artists.wireframe.set_segments(wSegs)
	artists.quiver.set_segments(qSegs)
	return artists

##Void Run Method for simulation
#Create Antenna (radius,theta,phi) Starting Antenna Position
tof =  86400 # in seconds
mps =  1/360 # Number of Simulated Location points per second
[constData,activeSats]=yuma_decode.gatherData()
satCoords = []
for satData in constData:
	satCoords.append(getSatPos(satData))
satCoords = np.stack(satCoords)
ant = [RADIUS_EARTH,np.radians(ANTENNA_LOCATION[0]),np.radians(ANTENNA_LOCATION[1])]
plt.style.use('dark_background')
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
earthMotion = plotEarth(ax,tof,mps)
antMotion = plotAntenna(ax,tof,mps,ant)
sigStr = 1.5
Nfrm = 240
fps = 24
def compute_segs(t):
	m,n,o,p,q,r = (earthMotion[t,0],earthMotion[t,1],earthMotion[t,2],earthMotion[t,3],earthMotion[t,4],earthMotion[t,5])
	x,y,z,u,v,w = (antMotion[t,0],antMotion[t,1],antMotion[t,2],sigStr*antMotion[t,0],sigStr*antMotion[t,1],sigStr*antMotion[t,2])
	return x,y,z,u,v,w,m,n,o,p,q,r

Artists = namedtuple("Artists", ("wireframe", "quiver"))

artists = Artists(
	ax.plot_wireframe(np.array([[]]),np.array([[]]),np.array([[]]),color="blue"),
	ax.quiver([],[],[],[],[],[],color="green"),
)


def frame_iter(from_second, until_second):
	for t in range(from_second, until_second):
		x,y,z,u,v,w,m,n,o,p,q,r = compute_segs(t)
		yield(x,y,z,u,v,w,m,n,o,p,q,r)
init = partial(init_fig, fig=fig, ax=ax, artists=artists)
step = partial(frame_iter, from_second=0, until_second=int(tof*mps))
update = partial(update_artists, artists=artists)
ani = animation.FuncAnimation(
	fig=fig, 
	func=update, 
	frames = step,
	interval=1000/fps,
	blit=True,
	init_func=init,
	save_count=len(list(step())),
    repeat_delay=0,
)

plt.show()
ani.save(
  filename='gps_sim.gif',
  fps=60,
  dpi=300,
)
