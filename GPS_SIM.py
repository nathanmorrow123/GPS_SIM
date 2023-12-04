#Nathan Morrow AFIT Fall 2023
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time
import math
RADIUS_EARTH=6.371e6 #meters
ROTATION_EARTH=7.2921150e-5 #rad/s
GRAV_CONSTANT=6.67430e-11 #N*m^2/kg^2
EARTH_MASS=5.97219e24 #kg
ANTENNA_LOCATION = [39.765659548785884, -84.19151722210734] #Latitude/Longitude in degrees

def run3dTrackingSim(times,satPath1,tof,mps):
	#Create Antenna (radius,theta,phi) Starting Antenna Position
	ant = [RADIUS_EARTH,np.radians(ANTENNA_LOCATION[0]),np.radians(ANTENNA_LOCATION[1])]
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	plt.ion()
	plt.show()
	for t in times:
		start = time.time()
		satX = satPath1[t,0]
		satY = satPath1[t,1]
		satZ = satPath1[t,2]
		#SignalStrength array (xSignal,ySignal,zSignal)
		sigStr = -(100)/updateAntennaStrength(ant,satX,satY,satZ)
		ant = updateAntennaPosition(ant,t,mps)
		ax.clear()
		ax.set_xlabel("X axis")
		ax.set_ylabel("Y Axis")
		ax.set_zlabel("Z Axis")
		ax.set_xlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		ax.set_ylim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		ax.set_zlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		plotSatellite(ax,satX,satY,satZ,satPath1)
		plotAntennaToSat(ax,ant,satX,satY,satZ)
		plotAntenna(ax,ant,sigStr)
		plotEarth(ax,t,mps)
		plt.pause(1e-3)
		end = time.time()
		while (end-start<1e-2):
			time.sleep(1e-3)
			end = time.time()
		plt.show()

def updateAntennaPosition(ant,t,mps):	
	ant[2]=-1*ROTATION_EARTH*t/mps
	return ant

def updateAntennaStrength(ant,satX,satY,satZ):
	#Formula for estimating RSSI is C*10 * Log10(Distance) where C is outdoor Coefficient Normally 2
	radius = math.sqrt(satX**2 + satY**2 + satZ**2)
	antStr = -25 * math.log((math.sqrt((satX-radius*math.sin(ant[1])*math.cos(ant[2]))**2+(satY-radius*math.sin(ant[1])*math.sin(ant[2]))**2+(satZ-radius*math.cos(ant[1]))**2)),10) 
	print(antStr)# Power to satellite in dB
	return antStr

#Set orbit of L1 GPS satillite using ECEI frame 
def setSatellitePath(tof,mps,alt,rot):
	times = range(int(tof*mps))#number of points to measure in a 10 second timespan
	times = np.divide(times,mps)
	r = alt+RADIUS_EARTH
	orbPeriod = 2*math.pi*math.sqrt(r**3/(GRAV_CONSTANT*EARTH_MASS))#Keeplers Equation
	satPosMat = []
	c, s = np.cos(rot), np.sin(rot)
	R = np.matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
	for t in times:
		v= np.matrix([[r * math.cos(-t/orbPeriod)],[r * math.sin(-t/orbPeriod)],[0]])
		satPosMat.append((R*v).T)
	satPath = np.stack(satPosMat)
	print(satPath)
	return satPath
def plotEarth(ax,t,mps):
	c, s = np.cos((t/mps)*ROTATION_EARTH), np.sin((t/mps)*ROTATION_EARTH)
	R = np.matrix([[c, s, 0], [-s, c, 0], [0, 0, 1]])
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = RADIUS_EARTH*np.cos(u)*np.sin(v)
	y = RADIUS_EARTH*np.sin(u)*np.sin(v)
	z = RADIUS_EARTH*np.cos(v)
	preImage = [x,y,z]
	postImage = [x*R[0,0]+y*R[0,1],x*R[1,0]+y*R[1,1],z]
	ax.plot_wireframe(postImage[0], postImage[1], postImage[2], color="b")
def plotAntenna(ax,ant,sigStr):
	antCart = sphericalTocartesian(ant)
	ax.quiver(antCart[0],antCart[1],antCart[2],sigStr*antCart[0],sigStr*antCart[1],sigStr*antCart[2],color = 'purple')#xAxis Rx Antenna 
def plotAntennaToSat(ax,ant,satX,satY,satZ):
	antCart = sphericalTocartesian(ant)
	ax.quiver(antCart[0],antCart[1],antCart[2],satX, satY, satZ, color = 'black', length = 0.5)#Tx Antenna 
def plotSatellite(ax,satX,satY,satZ,satPath1):
	ax.scatter(satX, satY, satZ, c='red', marker='.', s=20)
	ax.scatter(satPath1[:,0],satPath1[:,1],satPath1[:,2], c='red', marker='.', s=10)
def sphericalTocartesian(sCoords): #sCoords (rad,theta,phi)
	cartX=sCoords[0]*math.sin(sCoords[1])*math.cos(sCoords[2])
	cartY=sCoords[0]*math.sin(sCoords[1])*math.sin(sCoords[2])
	cartZ=sCoords[0]*math.cos(sCoords[1])
	return [cartX,cartY,cartZ]
##Void Run Method for simulation
tof =  10*86400 # in seconds
mps =  1/240 # Number of Simulatred Location points per second
satPath1 = setSatellitePath(tof,mps,2e7,np.radians(45))
satPath2 = setSatellitePath(tof,mps,3e7,np.radians(20))	
while True:
	try:
		times = range(int(tof*mps))
		run3dTrackingSim(times,satPath1,tof,mps)
	except OSError:
		print(OSError)
	pass
window.close()
