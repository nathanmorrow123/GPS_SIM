#Nathan Morrow AFIT Fall 2023
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time
import math
RADIUS_EARTH=6.371e6 #meters
ROTATION_EARTH=7.2921150e-5 #rad/s

def run3dTrackingSim(times,satCoords):
	X0=0
	Y0=0
	Z0=0
	#Create Antenna Spherical Coordinates (RSSI,theta,phi) Starting Antenna Positions || to Corresponding Axis
	xAnt = [10,math.pi/2,0]
	yAnt = [10,math.pi/2,math.pi/2]
	zAnt = [10,0,0]
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')
	plt.ion()
	plt.show()
	ax.set_xlabel("X axis")
	ax.set_ylabel("Y Axis")
	ax.set_zlabel("Z Axis")
	for t in times:
		satCoords = satPath1[:,t]
		satX = satCoords[0]
		satY = satCoords[1]
		satZ = satCoords[2]
		#SignalStrength array (xSignal,ySignal,zSignal)
		sigStr = updateAntennaStrength(xAnt,yAnt,zAnt,satX,satY,satZ)
		xAnt[0] = -1e9/(sigStr[0])
		yAnt[0] = -1e9/(sigStr[1])
		zAnt[0] = -1e9/(sigStr[2])
		ax.clear()
		ax.set_xlim([-2*RADIUS_EARTH,2*RADIUS_EARTH])
		ax.set_ylim([-2*RADIUS_EARTH,2*RADIUS_EARTH])
		ax.set_zlim([-2*RADIUS_EARTH,2*RADIUS_EARTH])
		plotSatellite(ax,satX,satY,satZ)
		plotAntennaToSat(ax,X0,Y0,Z0,satCoords)
		plotAntenna(ax,X0,Y0,Z0,xAnt,yAnt,zAnt)
		plotEarth(ax)
		plt.pause(.01)
		plt.show()
    
def updateAntennaPosition(xAnt,yAnt,zAnt):	
	thetaErr = 0
	phiErr = 0
	thetaPhiErr = [thetaErr,phiErr]
	return thetaPhiErr

def updateAntennaStrength(xAnt,yAnt,zAnt,satX,satY,satZ):
	
	#Formula for estimating RSSI is C*10 * Log10(Distance) where C is outdoor Coefficient Normally 2
    radius = math.sqrt(satX**2 + satY**2 + satZ**2)
    xStr = -25 * math.log((math.sqrt((satX-radius*math.sin(xAnt[1])*math.cos(xAnt[2]))**2+(satY-radius*math.sin(xAnt[1])*math.sin(xAnt[2]))**2+(satZ-radius*math.cos(xAnt[1]))**2)),10) 
    yStr = -25 * math.log((math.sqrt((satX-radius*math.sin(yAnt[1])*math.cos(yAnt[2]))**2+(satY-radius*math.sin(yAnt[1])*math.sin(yAnt[2]))**2+(satZ-radius*math.cos(yAnt[1]))**2)),10) 
    zStr = -25 * math.log((math.sqrt((satX-radius*math.sin(zAnt[1])*math.cos(zAnt[2]))**2+(satY-radius*math.sin(zAnt[1])*math.sin(zAnt[2]))**2+(satZ-radius*math.cos(zAnt[1]))**2)),10) 
    sigStr = [xStr,yStr,zStr]
    print(sigStr)
    return sigStr
    
#Set orbit of L1 GPS satillite using ECEI frame 
def setSatellitePath(tof,mps,alt,rot):
	times = range(tof*mps)#number of points to measure in a 10 second timespan
	times = np.divide(times,mps)
	r = alt+RADIUS_EARTH
	satPosMat = []
	c, s = np.cos(rot), np.sin(rot)
	R = np.matrix([[c, s, 0], [-s, c,0], [0,0,1]])
	for t in times:
		v= np.matrix[[r * math.cos(t)],[r * math.sin(t)],[0]]
		satPosMat=satPosMat+R*v
	print(satPosMat)
	return satPosMat
def plotEarth(ax):
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = RADIUS_EARTH*np.cos(u)*np.sin(v)
	y = RADIUS_EARTH*np.sin(u)*np.sin(v)
	z = RADIUS_EARTH*np.cos(v)
	ax.plot_wireframe(x, y, z, color="b")
def plotAntenna(ax,X0,Y0,Z0,xAnt,yAnt,zAnt):
	ax.quiver(X0,Y0,Z0,xAnt[0]*math.sin(xAnt[1])*math.cos(xAnt[2]),xAnt[0]*math.sin(xAnt[1])*math.sin(xAnt[2]),xAnt[0]*math.cos(xAnt[1]),color = 'red')#xAxis Rx Antenna 
	ax.quiver(X0,Y0,Z0,yAnt[0]*math.sin(yAnt[1])*math.cos(yAnt[2]),yAnt[0]*math.sin(yAnt[1])*math.sin(yAnt[2]),yAnt[0]*math.cos(yAnt[1]),color = 'green')#yAxis Rx Antenna 
	ax.quiver(X0,Y0,Z0,zAnt[0]*math.sin(zAnt[1])*math.cos(zAnt[2]),zAnt[0]*math.sin(zAnt[1])*math.sin(zAnt[2]),zAnt[0]*math.cos(zAnt[1]),color = 'blue')#zAxis Rx Antenna 
def plotAntennaToSat(ax,X0,Y0,Z0,satCoords):
	ax.quiver(X0,Y0,Z0,satCoords[0], satCoords[1], satCoords[2], color = 'black', length = 0.5)#Tx Antenna 
def plotSatellite(ax,satX,satY,satZ):
	ax.scatter(satX, satY, satZ, c='red', marker='.', s=20)
	ax.scatter(satPath1[0],satPath1[1], satPath1[2], c='red', marker='.', s=5)
##Void Run Method for simulation
tof = 20 # Time of Flight = 10 Seconds
mps = 10 # Number of Simulatred Location points per second
satPath1 = setSatellitePath(tof,mps,2e7,np.radians(20))
satPath2 = setSatellitePath(tof,mps,3e7,np.radians(20))	
while True:
    try:
        times = range(tof*mps)
        run3dTrackingSim(times,satPath1)
    except OSError:
        print(OSError)
        pass
    window.close()

