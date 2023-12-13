#Nathan Morrow AFIT Fall 2023
#All GPS data is downloaded from NAVCEN
#https://www.navcen.uscg.gov/gps-nanus-almanacs-opsadvisories-sof

from matplotlib import pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
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
def setSatellitePath(tof,mps,alt,theta,phi):
	times = range(int(tof*mps))#number of points to measure in a 10 second timespan
	times = np.divide(times,mps)
	r = alt+RADIUS_EARTH
	orbPeriod = 2*math.pi*math.sqrt(r**3/(GRAV_CONSTANT*EARTH_MASS))#Keeplers Equation
	satPosMat = []
	c1, s1 = np.cos(theta), np.sin(theta)
	c2, s2 = np.cos(phi), np.sin(phi)
	R1 = np.matrix([[1, 0, 0], [0, c2, -s2], [0, s2, c2]])
	R2 = np.matrix([[c1, 0, s1], [0, 1, 0], [-s1, 0, c1]])
	for t in times:
		v= np.matrix([[r * math.cos(-t/orbPeriod)],[r * math.sin(-t/orbPeriod)],[0]])
		satPosMat.append((R2*(R1*v)).T)
	satPath = np.stack(satPosMat)
	print(satPath)
	return satPath
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
mps =  1/24 # Number of Simulatred Location points per second
[constData,activeSats]=yuma_decode.gatherData()
satCoords = []
for satData in constData:
	satCoords.append(getSatPos(satData))
satCoords = np.stack(satCoords)
ant = [RADIUS_EARTH,np.radians(ANTENNA_LOCATION[0]),np.radians(ANTENNA_LOCATION[1])]
plt.style.use('dark_background')
fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
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
earthMotion = plotEarth(ax,tof,mps)
antMotion = plotAntenna(ax,tof,mps,ant)
sigStr = 1.2
Nfrm = 3600
fps = 60
wframe = None
quiver = None
def init():
	t=0
	plotAntenna(ax,tof,mps,ant)
	return (ax.quiver(antMotion[t,0],antMotion[t,1],antMotion[t,2],sigStr*antMotion[t,0],sigStr*antMotion[t,1],sigStr*antMotion[t,2],color = 'green'),ax.plot_wireframe(earthMotion[t,0], earthMotion[t,1], earthMotion[t,2], color="b"))#normal)
def update(t):
	global wframe
	global quiver
	# If a line collection is already remove it before drawing.
	if wframe:
		ax.collections.remove(wframe)
	if quiver:
		ax.collections.remove(quiver)
	quiver = ax.quiver(antMotion[t,0],antMotion[t,1],antMotion[t,2],sigStr*antMotion[t,0],sigStr*antMotion[t,1],sigStr*antMotion[t,2],color = 'green')#normal
	wframe = ax.plot_wireframe(earthMotion[t,0], earthMotion[t,1], earthMotion[t,2], color="b")
ani = animation.FuncAnimation(fig, update, interval=10000/fps,blit=True,init_func=init)

plt.show()
fn = 'plot_wireframe_funcanimation'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)
ani.save(fn+'.gif',writer='imagemagick',fps=fps)
