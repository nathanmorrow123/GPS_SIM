#Nathan Morrow AFIT Fall 2023
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time
import math
RADIUS_EARTH=6.371e6 #meters
ROTATION_EARTH=7.2921150e-5 #rad/s

def run3dTrackingSim(times,satPath1,tof,mps):
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
		satX = satPath1[t,0]
		satY = satPath1[t,1]
		satZ = satPath1[t,2]
		#SignalStrength array (xSignal,ySignal,zSignal)
		sigStr = updateAntennaStrength(xAnt,yAnt,zAnt,satX,satY,satZ)
		xAnt[0] = -1e9/(sigStr[0])
		yAnt[0] = -1e9/(sigStr[1])
		zAnt[0] = -1e9/(sigStr[2])
		ax.clear()
		ax.set_xlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		ax.set_ylim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		ax.set_zlim([-3*RADIUS_EARTH,3*RADIUS_EARTH])
		plotSatellite(ax,satX,satY,satZ,satPath1)
		plotAntennaToSat(ax,X0,Y0,Z0,satX,satY,satZ)
		plotAntenna(ax,X0,Y0,Z0,xAnt,yAnt,zAnt)
		plotEarth(ax,t,mps)
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
	times = range(int(tof*mps))#number of points to measure in a 10 second timespan
	times = np.divide(times,mps)
	r = alt+RADIUS_EARTH
	satPosMat = []
	c, s = np.cos(rot), np.sin(rot)
	R = np.matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
	for t in times:
		v= np.matrix([[r * math.cos(t)],[r * math.sin(t)],[0]])
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
def plotAntenna(ax,X0,Y0,Z0,xAnt,yAnt,zAnt):
	ax.quiver(X0,Y0,Z0,xAnt[0]*math.sin(xAnt[1])*math.cos(xAnt[2]),xAnt[0]*math.sin(xAnt[1])*math.sin(xAnt[2]),xAnt[0]*math.cos(xAnt[1]),color = 'red')#xAxis Rx Antenna 
	ax.quiver(X0,Y0,Z0,yAnt[0]*math.sin(yAnt[1])*math.cos(yAnt[2]),yAnt[0]*math.sin(yAnt[1])*math.sin(yAnt[2]),yAnt[0]*math.cos(yAnt[1]),color = 'green')#yAxis Rx Antenna 
	ax.quiver(X0,Y0,Z0,zAnt[0]*math.sin(zAnt[1])*math.cos(zAnt[2]),zAnt[0]*math.sin(zAnt[1])*math.sin(zAnt[2]),zAnt[0]*math.cos(zAnt[1]),color = 'blue')#zAxis Rx Antenna 
def plotAntennaToSat(ax,X0,Y0,Z0,satX,satY,satZ):
	ax.quiver(X0,Y0,Z0,satX, satY, satZ, color = 'black', length = 0.5)#Tx Antenna 
def plotSatellite(ax,satX,satY,satZ,satPath1):
	ax.scatter(satX, satY, satZ, c='red', marker='.', s=20)
	ax.scatter(satPath1[:,0],satPath1[:,1],satPath1[:,2], c='red', marker='.', s=5)
##Void Run Method for simulation
tof =  864000 # in seconds
mps =  1/3600 # Number of Simulatred Location points per second
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

