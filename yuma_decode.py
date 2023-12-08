# Program to read the specific length 
# of characters in a yuma using read() function
import os
import re
import numpy as np
import pandas as pd

def gatherData():
    
    os.system("sudo ./get_data.sh")
    yuma = open("current_yuma.alm", "r")
    opsadvisory = open("current_opsadvisory.txt", "r")
    content1 = yuma.readlines()
    content2 = opsadvisory.readlines()
    rows = 32
    cols = 12
    constData=np.zeros((rows,cols))
    satID = 0

    #Example Yuma output
    #******** Week 244 almanac for PRN-02 ********
    #ID:                         02
    #Health:                     000
    #Eccentricity:               0.1624822617E-001
    #Time of Applicability(s):  61440.0000
    #Orbital Inclination(rad):   0.9674049839
    #Rate of Right Ascen(r/s):  -0.8046049436E-008
    #SQRT(A)  (m 1/2):           5153.857422
    #Right Ascen at Week(rad):  -0.1336742916E+001
    #Argument of Perigee(rad):  -1.292559701
    #Mean Anom(rad):            -0.1510177121E+001
    #Af0(s):                    -0.5178451538E-003
    #Af1(s/s):                   0.3637978807E-011
    #week:                        244
    
    sciReg = '[+\-]?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?'
    for line in content1:
        if "ID" in line:
            constData[satID][0] = int(re.findall(r'\d+',line.split()[1])[0])
        if "Health" in line:
            constData[satID][1] = int(re.findall(r'\d+',line.split()[1])[0])
        if "Eccentricity" in line:
            constData[satID][2] = re.findall(sciReg,line.split()[1])[0]
        if "Time of Applicability(s)" in line:
            constData[satID][3] = re.findall(sciReg,line.split()[3])[0]
        if "Orbital Inclination(rad)" in line:
            constData[satID][4] = re.findall(sciReg,line.split()[2])[0]
        if "Rate of Right Ascen(r/s)" in line:
            constData[satID][5] = re.findall(sciReg,line.split()[4])[0]
        if "SQRT(A)  (m 1/2)" in line:
            constData[satID][6] = re.findall(sciReg,line.split()[3])[0]
        if "Right Ascen at Week(rad)" in line:
            constData[satID][7] = re.findall(sciReg,line.split()[4])[0]
        if "Argument of Perigee(rad)" in line:
            constData[satID][8] = re.findall(sciReg,line.split()[3])[0]
        if "Mean Anom(rad)" in line:
            constData[satID][9] = re.findall(sciReg,line.split()[2])[0]
        if "Af0(s)" in line:
            constData[satID][10] = re.findall(sciReg,line.split()[1])[0]
        if "Af1(s/s)" in line:
            constData[satID][11] = re.findall(sciReg,line.split()[1])[0]
            satID += 1
    print(pd.DataFrame(constData))
    activeSats = np.zeros(32)
    i = 0
    #Reading in the Opsadvisory for active satellites
    for line in content2:
        if "BLOCK II : PRNS" in line:
            for str in re.findall(r'\d+',line):
                activeSats[i] = int(str)
                i+=1
        elif "BLOCK III: PRNS" in line:
            for str in re.findall(r'\d+',line):
                activeSats[i] = int(str)
                i+=1
    print(activeSats)
    yuma.close()
    opsadvisory.close()
    return [constData,activeSats]

if __name__ == '__main__': 
    gatherData()

