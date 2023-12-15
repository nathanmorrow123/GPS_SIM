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

def getSatPos(args):
	[satID,health,e,toa,inc,ascRate,sqrt_a,ra_week,w,E,Af0,Af1] = args
	a = sqrt_a**2
	b = math.sqrt(a**2*(1-e**2)) #https://jtauber.github.io/orbits/019.html
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
def plotEarth(ax,tof,mps):
	times = range(int(tof*mps))
	earthMovement = []
	for t in times:
		c, s = np.cos((t/mps)*ROTATION_EARTH), np.sin((t/mps)*ROTATION_EARTH)
		R = np.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		x = RADIUS_EARTH*np.cos(u)*np.sin(v)
		y = RADIUS_EARTH*np.sin(u)*np.sin(v)
		z = RADIUS_EARTH*np.cos(v)
		preImage = [x,y,z]
		postImage = [x*R[0,0]+y*R[0,1],x*R[1,0]+y*R[1,1],z]
		earthMovement.append(postImage)
	return np.stack(earthMovement)
def plotAntenna(ax,tof,mps,ant):
	times = range(int(tof*mps))
	antennaMovement = []
	antCart = sphericalTocartesian(ant)
	for t in times:
		ant[2]=-1*ROTATION_EARTH*t/mps
		antCart = sphericalTocartesian(ant)
		antennaMovement.append([antCart[0],antCart[1],antCart[2]])
	return np.stack(antennaMovement)
def plotSats(ax,satCoords):
	for satCoord in satCoords:
		ax.scatter(satCoord[1],satCoord[2],satCoord[3],c='green',marker='.',s=20)
		ax.text(satCoord[1],satCoord[2],satCoord[3],int(satCoord[0]))
def sphericalTocartesian(sCoords): #sCoords (rad,theta,phi)
	cartX=sCoords[0]*math.sin(sCoords[1])*math.cos(sCoords[2])
	cartY=sCoords[0]*math.sin(sCoords[1])*math.sin(sCoords[2])
	cartZ=sCoords[0]*math.cos(sCoords[1])
	return [cartX,cartY,cartZ]

##Void Run Method for simulation
#Create Antenna (radius,theta,phi) Starting Antenna Position
tof =  86400 # in seconds
mps =  1/360 # Number of Simulatred Location points per second
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

	a,b,c = (earthMotion[t,0],earthMotion[t,1],earthMotion[t,2])
	x,y,z,u,v,w = (antMotion[t,0],antMotion[t,1],antMotion[t,2],sigStr*antMotion[t,0],sigStr*antMotion[t,1],sigStr*antMotion[t,2])
	return x,y,z,u,v,w,a,b,c

Artists = namedtuple("Artists", ("wireframe", "quiver"))

artists = Artists(
	ax.plot_wireframe(np.array([[]]),np.array([[]]),np.array([[]]),color="blue"),
	ax.quiver([],[],[],[],[],[],color="green"),
)

def init_fig(fig,ax,artists):
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
	plotSats(ax,satCoords)
	return artists
def meshTolines(a,b,c):
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
def update_artists(frames,artists):
	x,y,z,u,v,w,a,b,c = frames
	m,n,o,p,q,r = meshTolines(a,b,c)
	temp = np.array([x,y,z,u,v,w]).reshape(6,-1)
	qSegs = [[[x,y,z],[u,v,w]]for x,y,z,u,v,w in zip(*temp.tolist())]
	temp = np.array([m,n,o,p,q,r]).reshape(6,-1)
	wSegs = [[[m,n,o],[p,q,r]]for m,n,o,p,q,r in zip(*temp.tolist())]
	artists.wireframe.set_segments(wSegs)
	artists.quiver.set_segments(qSegs)
	return artists
def frame_iter(from_second, until_second):
	for t in range(from_second, until_second):
		x,y,z,u,v,w,a,b,c = compute_segs(t)
		yield(x,y,z,u,v,w,a,b,c)
init = partial(init_fig, fig=fig, ax=ax, artists=artists)
step = partial(frame_iter, from_second=0, until_second=int(tof*mps))
update = partial(update_artists, artists=artists)
ani = animation.FuncAnimation(
	fig=fig, 
	func=update, 
	frames = step,
	interval=1000/fps,
	blit=False,
	init_func=init,
	save_count=len(list(step())),
    repeat_delay=0,
)

plt.show()
ani.save(
  filename='/tmp/gps_sim.mp4',
  fps=24,
  extra_args=['-vcodec', 'libx264'],
  dpi=300,
)
