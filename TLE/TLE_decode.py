# Program to read the specific length 
# of characters in a yuma using read() function
import os
import re
import numpy as np
import pandas as pd

def gatherData():
    
    os.system("sudo ./get_data.sh")
    tle = open("gnss_tle.txt.alm", "r")
    content1 = tle.readlines()
    rows = 32
    cols = 12
    constData=np.zeros((rows,cols))
    satID = 0
    
    #Example TLE output from celestrak
    '''
    NAVSTAR 43 (USA 132)    
    1 24876U 97035A   24046.62828655 -.00000021  00000+0  00000+0 0  9991
    2 24876  55.6310 131.9233 0076080  52.7907 307.9354  2.00565331194855
    '''

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
    tle.close()
    return [constData,activeSats]

if __name__ == '__main__': 
    gatherData()

