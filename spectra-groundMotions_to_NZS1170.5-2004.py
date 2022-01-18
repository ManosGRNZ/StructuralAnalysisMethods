# Export the spectra graphs for the K1 (and optionally K2) factors identified

import math
import multiprocessing
from numpy.lib.function_base import average, disp
import requests
import matplotlib.pyplot as plt
from numpy import arange
from time import time
import pickle
from os import path

# Ground motion Records to run
numGMsRun = 14  # usually up to 11 - NZS1170.5 requires minimum 3
processGMList = list( range(numGMsRun) )

startTime = time()

Z = 0.40 * 9.81  # Hazard factor for Wellington, New Zealand

def sign(x):
    if x == 0:
        return 0
    else:
        return x / abs(x)

# Increment and range of natural periods for spectrum (secs)
T1 = 1.5
minT = 0.4 * T1  # should not be zero if mass is given and stiffness is calculated - avoid div by zero
maxT = 1.3 * T1
numOfT = 100
deltaT = (maxT - minT) / numOfT

# Set filenames here for the two directions of the same ground motion
path_Linux = "/home/user/YourPathGoesHere/"
path_Windows = "C:\\YourPathGoesHere\\"

if path.exists(path_Linux):
    pathMain = path_Linux
    pathOut = pathMain + "output/"
if path.exists(path_Windows):
    pathMain = path_Windows
    pathOut = pathMain + "output\\"

#list_of_GMs = ["RSN171", "RSN722", "RSN767", "RSN1084", "RSN1176", "RSN1481", "RSN1602", "RSN6960", "RSN8606", "RSN4032588", "RSN4032590" ]
#scaleFactorsGM = [0.93, 2.67, 0.90, 0.71, 1.72, 2.80, 0.77, 1.70, 1.24, 2.70, 2.00 ]

scaleFactorsGM = pickle.load( open(pathMain+"k1.dat", "rb") )
list_of_GMs_1 = pickle.load( open(pathMain+"GM1-List.dat", "rb") )
list_of_GMs_2 = pickle.load( open(pathMain+"GM2-List.dat", "rb") )

# Remove the last character of each filename (which is "1" or "2") to have a single name for each ground motion
list_of_GMs_names = list_of_GMs_1.copy()
for i in range( len(list_of_GMs_names) ) :
    list_of_GMs_names[i] = list_of_GMs_names[i][:-1]
    
damping = 0.05  # ratio of critical damping of the SDOF
alpha = 1  # damping exponent


k = 10.0 ** 5  # default stiffness in N/m - arbitrarily set, as the other values will result
# m = 100 # default mass in kg - arbitrarily set, as the other values will result

#Tspec = pickle.load( open(pathMain+"periodsSDOF.dat", "rb"))
Tspec = [minT]
for T in arange(minT, maxT, deltaT):
    Tspec.append(T+deltaT)
filename_Periods = pathOut + "periods.data"
pickle.dump(Tspec, open(filename_Periods, 'wb'))


def runGM(gmNum, avList, avListX, avListY):
    tempStartTime = time()

    gmID = gmNum
    accGMfileFP = open( pathMain+list_of_GMs_1[gmID] , "r")
    accGMfileFN = open( pathMain+list_of_GMs_2[gmID] , "r")

    # Get dt from one file and ignore the first line of the second file (its the same)
    str = accGMfileFP.readline()
    dt = float(str.split()[-1]) #ignore the first number, the second is the time intervalnext(accGMfileFN)
    str = accGMfileFN.readline()
    
    # multiplier for the ground motion numbers (e.g. 9.81 if the units are g)
    multiplier = 9.81 * scaleFactorsGM[gmID]

    # acceleration records are imported as lists
    accX = []
    accY = []

    for ax in accGMfileFP:
        accX.append(float(ax.split()[-1]) * multiplier)

    for ay in accGMfileFN:
        accY.append(float(ay.split()[-1]) * multiplier)

    accGMfileFP.close()
    accGMfileFN.close()

    # print(accX ,accY)

    # Integration of the two records into ground velocities and then ground displacements
    velX = [0.0]
    velY = [0.0]
    dispX = [0.0]
    dispY = [0.0]

    for a in accX:
        velX.append(velX[-1] + a * dt)
        dispX.append(dispX[-1] + velX[-1] * dt)

    for a in accY:
        velY.append(velY[-1] + a * dt)
        dispY.append(dispY[-1] + velY[-1] * dt)

    recNum = len(velX)

    #print(f"{pathMain+list_of_GMs_1[gmID]} \n velX: {len(velX)} \t velY: {len(velY)} \t dispX: {len(dispX)} \t dispY: {len(dispY)} \n")
    if len(velX) != len(velY):
        print(f"CHECK THIS --> {pathMain+list_of_GMs_1[gmID]} \n")

    # Lists of periods and accelerations for the response spectrum
    global Tspec

    #Initialise spectrums for T=0
    Axspec = [float(max(accX))]
    Ayspec = [float(max(accY))]
    Aspec = [( (Axspec[0])**2 + (Ayspec[0])**2 )**0.5 ]

    # iterate through the range of SDOFs to form response spectrum
    for T in arange(minT, maxT, deltaT):
        # mass calculation in kg for the default k and this T
        m = (T / (2 * math.pi)) ** 2 * k
        # k = m / ( (T / (2 * math.pi)) ** 2 ) # stiffness calculation for the default m and this T
        c = 2 * ((k * m) ** 0.5) * damping  # damping coefficient for each SDOF
        # print("T = {:.2f} \t k = {:.2f} \t c = {:.2f}".format(T,k,c))

        # Initialise
        Xa = Xv = Xd = 0
        Ya = Yv = Yd = 0
        Xdg = Ydg = 0

        maxA = maxAx = maxAy = 0

        # Apply equation of ground motion equilibrium
        # Spring acting for rectangular response (no biaxial interaction)
        for n in range(0, recNum, 1):
            Xdg = dispX[n]
            Ydg = dispY[n]

            # Calculation of spring force and then relative acceleration, velocity and displecement
            Fkx = (Xdg - Xd) * k
            Xa = (Fkx - c * (abs(Xv) ** alpha) * sign(Xv)) / m
            Xv += Xa * dt
            Xd += Xv * dt

            Fky = (Ydg - Yd) * k
            Ya = (Fky - c * (abs(Yv) ** alpha) * sign(Yv)) / m
            Yv += Ya * dt
            Yd += Yv * dt

            # maximum acceleration regardless of angle
            maxA = max(maxA, (Xa ** 2 + Ya ** 2) ** 0.5)
            maxAx = max(maxAx, abs(Xa))
            maxAy = max(maxAy, abs(Ya))

        Aspec.append(maxA)
        Axspec.append(maxAx)
        Ayspec.append(maxAy)

        if len(avList) < len(Aspec):
            avList.append(Aspec[-1]/numGMsRun)
            avListX.append(Axspec[-1]/numGMsRun)
            avListY.append(Ayspec[-1]/numGMsRun)
        else:
            i = len(Aspec) - 1
            avList[i] += Aspec[-1]/numGMsRun
            avListX[i] += Axspec[-1]/numGMsRun
            avListY[i] += Ayspec[-1]/numGMsRun

    tempTime = time() - tempStartTime

    print("Time to calculate {}: {:.2f}".format(list_of_GMs_names[gmNum], tempTime))
    
    fig, axs = plt.subplots(4, facecolor=(.2, .2, .2))
    #axs[0].plot(Tspec, Aspec, label='Resultant', linewidth=2, color='yellow')
    axs[0].plot(Tspec, Axspec, label='Spectral-X', linewidth=1, color='red')
    axs[0].plot(Tspec, Ayspec, label='Spectral-Y', linewidth=1, color='green')
    #print(f"len(Tspec) = {len(Tspec)} \t len(spectral_C) = {len(spectral_Shape_factor_class_C(Tspec))}")
    axs[0].plot(Tspec, spectral_Shape_factor_class_C(Tspec), label='Design Spectrum', linewidth=1, color='blue')
    axs[1].plot(accX, label='acceleration-X', linewidth=0.5, color='red')
    axs[1].plot(accY, label='acceleration-Y', linewidth=0.5, color='green')
    
    axs[2].plot(velX, label='velocity-X', linewidth=0.5, color='red')
    axs[2].plot(velY, label='velocity-Y', linewidth=0.5, color='green')
    axs[3].plot(dispX, label='displacement-X', linewidth=0.5, color='red')
    axs[3].plot(dispY, label='displacement-Y', linewidth=0.5, color='green')
    
    axs[0].set_facecolor('black')
    axs[1].set_facecolor('black')
    axs[2].set_facecolor('black')
    axs[3].set_facecolor('black')
    fig.suptitle(list_of_GMs_names[gmID], fontsize=12)

    #plt.show()
    plt.savefig(pathOut + list_of_GMs_names[gmID] + '.png', dpi=300)
    
    filename_FP = pathOut + list_of_GMs_names[gmID] + "_Spectrum_FP.data"
    filename_FN = pathOut + list_of_GMs_names[gmID] + "_Spectrum_FN.data"
    filename_Resultant = pathOut + list_of_GMs_names[gmID] + "_Spectrum_Resultant.data"
    pickle.dump(Axspec, open(filename_FP, 'wb'))
    pickle.dump(Ayspec, open(filename_FN, 'wb'))
    pickle.dump(Aspec, open(filename_Resultant, 'wb'))

    return tempTime


def spectral_Shape_factor_class_C(periods):
    Ch = []
    for Ti in periods:
        if Ti == 0:
            Ch.append(Z * 1.33)
        elif Ti < 0.1:
            Ch.append(Z * 1.33 + 1.60 * (Ti / 0.1))
        elif Ti < 0.3:
            Ch.append(Z * 2.93)
        elif Ti <= 1.5:
            Ch.append(Z * 2 * (0.5 / Ti)**0.75)
        elif Ti <= 3:
            Ch.append(Z * 1.32 / Ti)
        else:
            Ch.append(Z * 3.96 / (Ti**2))
    return Ch


def spectral_Shape_factor_class_D(periods):
    Ch = []
    for Ti in periods:
        if Ti == 0:
            Ch.append(Z * 1.12)
        elif Ti < 0.1:
            Ch.append(Z * 1.12 + 1.88 * (Ti / 0.1))
        elif Ti < 0.56:
            Ch.append(Z * 3.0)
        elif Ti <= 1.5:
            Ch.append(Z * 2.4 * (0.75 / Ti)**0.75)
        elif Ti <= 3:
            Ch.append(Z * 2.14 / Ti)
        else:
            Ch.append(Z * 6.42 / (Ti**2))
    return Ch


def spectral_Shape_factor_interpolation_C_D(periods, Tsite):
    Ch = []
    fi = 1.05 + 0.5 * (Tsite - 0.25)
    for Ti in periods:
        if Ti ==0 :
            Ch.append(Z * 1.33)
        elif Ti < 0.1:
            Ch.append(Z * 1.33 + 1.60 / (Ti / 0.1))
        elif Ti < 0.3:
            Ch.append(min(3 * Z, fi * Z * 2.93))
        elif Ti <= 1.5:
            Ch.append(min(3 * Z, fi * Z * 2 * (0.5 / Ti)**0.75))
        elif Ti <= 3:
            Ch.append(fi * Z * 1.32 / Ti)
        else:
            Ch.append(fi * Z * 3.96 / (Ti**2))
    return Ch


if __name__ == '__main__':
    manager = multiprocessing.Manager()
    averageSpectrum = manager.list()
    averageSpectrum_X = manager.list()
    averageSpectrum_Y = manager.list()

    print("multiprocessing starts")
    processes = []
    for process in processGMList:
        print('append: ',process)
        processes.append(multiprocessing.Process(target=runGM, args=[process, averageSpectrum, averageSpectrum_X, averageSpectrum_Y]))

    for process in processes:
        print('start: ',process)
        process.start()

    for process in processes:
        print('join: ',process)
        process.join()

    fig, axs = plt.subplots(3, facecolor=(.2, .2, .2))
    #axs[0].plot(Tspec, averageSpectrum, label='Resultant', linewidth=2, color='yellow')
    axs[0].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[1].plot(Tspec, averageSpectrum_X, label='Spectral-FP', linewidth=1, color='red')
    axs[1].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[2].plot(Tspec, averageSpectrum_Y, label='Spectral-FN', linewidth=1, color='green')
    axs[2].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[0].set_facecolor('black')
    axs[1].set_facecolor('black')
    axs[2].set_facecolor('black')
    fig.suptitle("Accumulated Average", fontsize=12)
    #plt.show()
    plt.savefig(pathOut + "average-spectra" + '.png')

print("Total time: {:.2f}".format(time() - startTime))
