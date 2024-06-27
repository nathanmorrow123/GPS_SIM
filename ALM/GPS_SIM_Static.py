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

def compute_mean_anomaly(E, e, t, toa, GRAV_CONSTANT, EARTH_MASS):
    a = E ** 2
    n = math.sqrt(GRAV_CONSTANT * EARTH_MASS / a**3)  # Mean motion
    M0 = E - e * math.sin(E)  # Mean anomaly at reference time
    return M0 + n * (t - toa)  # Mean anomaly at time t

def solve_keplers_equation(M, e, iterations=10):
    E_new = M
    for _ in range(iterations):  # Iterate to solve for E
        E_new = M + e * math.sin(E_new)
    return E_new

def compute_orbital_position(a, e, E_new):
    r = a * (1 - e * math.cos(E_new))
    v = 2 * math.atan2(math.sqrt(1 + e) * math.sin(E_new / 2), math.sqrt(1 - e) * math.cos(E_new / 2))
    x_orb = r * math.cos(v)
    y_orb = r * math.sin(v)
    return x_orb, y_orb

def rotate_position(x_orb, y_orb, w, inc, ra_week):
    x_prime = math.cos(w) * x_orb - math.sin(w) * y_orb
    y_prime = math.sin(w) * x_orb + math.cos(w) * y_orb
    z_prime = math.sin(inc) * y_prime
    y_final = math.cos(inc) * y_prime
    x_final = math.cos(ra_week) * x_prime - math.sin(ra_week) * y_final
    y_final = math.sin(ra_week) * x_prime + math.cos(ra_week) * y_final
    return x_final, y_final, z_prime

def compute_sat_pos(args, t):
    [satID, health, e, toa, inc, ascRate, sqrt_a, ra_week, w, E, Af0, Af1] = args
    t = t + toa  # Sync Animation to TOA
    a = sqrt_a ** 2
    M = compute_mean_anomaly(E, e, t, toa, GRAV_CONSTANT, EARTH_MASS)
    E_new = solve_keplers_equation(M, e)
    x_orb, y_orb = compute_orbital_position(a, e, E_new)
    x_final, y_final, z_prime = rotate_position(x_orb, y_orb, w, inc, ra_week)
    return [satID, x_final, y_final, z_prime]

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

def computeAllEarthPos(times):
	"""
    Plot the Earth rotation as line segments for the given time of flight.

    Parameters:
    ax (Axes3D): Matplotlib 3D axes object.
    tof (float): Time of flight in seconds.
    mps (float): Number of simulated location points per second.

    Returns:
    ndarray: Array of line segments for each time step.
    """
	earth_segs = []
	for t in times:
		c, s = np.cos((t)*ROTATION_EARTH), np.sin((t)*ROTATION_EARTH)
		R = np.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		x = RADIUS_EARTH*np.cos(u)*np.sin(v)
		y = RADIUS_EARTH*np.sin(u)*np.sin(v)
		z = RADIUS_EARTH*np.cos(v)
		preImage = [x,y,z]
		postImage = [x*R[0,0]+y*R[0,1],x*R[1,0]+y*R[1,1],z]
		earth_segs.append(meshTolines(postImage[0],postImage[1],postImage[2]))
	return np.stack(earth_segs)

def computeAllAntPos(times,ant):
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
	antennaMovement = []
	antCart = sphericalTocartesian(ant)
	for t in times:
		ant[2]=-1*ROTATION_EARTH*t
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
		#ax.text(satCoord[1],satCoord[2],satCoord[3],int(satCoord[0]))

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
def init_fig(ax,artists,satCoords):
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
	plotSats(ax,satCoords) # Skip the first one it's bad data
	return artists

def update_artists(frame,artists):
	"""
    Update the artists for each animation frame.

    Parameters:
    frames (tuple): Frame data containing segments for wireframe and quiver plots.
    artists (Artists): Named tuple containing the wireframe and quiver artists.

    Returns:
    Artists: The updated artists.
    """
	earth_segs, sat_coord = frame
	x, y, z, u, v, w, m, n, o, p, q, r = earth_segs
	temp = np.array([x,y,z,u,v,w]).reshape(6,-1)
	qSegs = [[[x,y,z],[u,v,w]]for x,y,z,u,v,w in zip(*temp.tolist())]
	temp = np.array([m,n,o,p,q,r]).reshape(6,-1)
	wSegs = [[[m,n,o],[p,q,r]]for m,n,o,p,q,r in zip(*temp.tolist())]
	artists.wireframe.set_segments(wSegs)
	artists.quiver.set_segments(qSegs)
	artists.sat_pos.set_offsets([sat_coord[:, 1], sat_coord[:, 2], sat_coord[:, 3]])
	return artists

def get_segs(i):
	sigStr = 1.5
	m,n,o,p,q,r = (earthMotion[i,0],earthMotion[i,1],earthMotion[i,2],earthMotion[i,3],earthMotion[i,4],earthMotion[i,5])
	x,y,z,u,v,w = (antMotion[i,0],antMotion[i,1],antMotion[i,2],sigStr*antMotion[i,0],sigStr*antMotion[i,1],sigStr*antMotion[i,2])
	return x,y,z,u,v,w,m,n,o,p,q,r

def get_satPos(i):
	satCoords = satCoordsOverTime[i]
	return satCoords

def compute_all_sat_pos(const_data, times):
    sat_coords_over_time = []
    for t in times:
        sat_coords = [compute_sat_pos(sat_data, t) for sat_data in const_data]
        sat_coords_over_time.append(sat_coords)
    return np.array(sat_coords_over_time)

def frame_iter(times):
	# Helper funciton in animation
	for i in range(len(times)):
		x,y,z,u,v,w,m,n,o,p,q,r = get_segs(i)
		earth_segs = [x,y,z,u,v,w,m,n,o,p,q,r]
		satPos = get_satPos(i)
		yield(earth_segs,satPos)

##Void Run Method for simulation
#Create Antenna (radius,theta,phi) Starting Antenna Position
tof =  86400 # Total sim time (Seconds)
dt =  360 # Time step (seconds)
tof = 10
dt = 1
times = np.linspace(0,tof,int(tof/dt+1),endpoint=True)
[constData,activeSats]=yuma_decode.gatherData()
satCoordsOverTime = compute_all_sat_pos(constData, times)
ant = [RADIUS_EARTH,np.radians(ANTENNA_LOCATION[0]),np.radians(ANTENNA_LOCATION[1])]
earthMotion = computeAllEarthPos(times)
antMotion = computeAllAntPos(times,ant)

Nfrm = 240
fps = 24
Artists = namedtuple("Artists", ("wireframe", "quiver",'sat_pos'))
plt.style.use('dark_background')
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
artists = Artists(
	ax.plot_wireframe(np.array([[]]),np.array([[]]),np.array([[]]),color="blue"),
	ax.quiver([],[],[],[],[],[],color="green"),
	ax.scatter([],[],[],c='green',marker='.',s=20),
)

init = partial(init_fig, ax=ax, artists=artists, satCoords=satCoordsOverTime[0])
step = partial(frame_iter, times)
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
  dpi=100  
)
