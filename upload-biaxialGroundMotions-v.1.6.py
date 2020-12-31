# Acceleration response spectrum from biaxial ground motions
# Warning-1: This methodology does not aim to follow the strict theory of dynamics of structures
# Warning-2: This methodology is for demonstration and research purposes only
# Warning-3: DO NOT USE FOR ANY ACTUAL PROJECT
# Tested with Python 3.9

import math
import multiprocessing
from numpy.lib.function_base import average
import requests
import matplotlib.pyplot as plt
from numpy import arange
from time import time

# Ground motion Records to run
#processGMList = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
processGMList = [ 0, 1, 7 ]
numGMsRun = len(processGMList)

startTime = time()

Z = 0.40 * 9.81  # Hazard factor for Wellington

def sign(x):
    if x == 0:
        return 0
    else:
        return x / abs(x)


# Set filenames here for the two directions of the same ground motion
path = "/home/" # Set file path for ground motion files here - set the ground motions and scale factors below
list_of_GMs = ["RSN171", "RSN722", "RSN767", "RSN1084", "RSN1176", "RSN1481", "RSN1602", "RSN6960", "RSN8606", "RSN4032588", "RSN4032590" ]
scaleFactorsGM = [0.93, 2.17, 0.95, 0.91, 2.12, 1.80, 0.56, 1.33, 1.18, 2.10, 2.41 ]

damping = 0.05  # ratio of critical damping of the SDOF
alpha = 1  # damping exponent

# Increment and range of natural periods for spectrum (secs)
minT = 0.1  # should not be zero if mass is given and stiffness is calculated - avoid div by zero
maxT = 4.0 + minT
numOfT = 100
deltaT = (maxT - minT) / numOfT

k = 10.0 ** 5  # default stiffness in N/m - arbitrarily set, as the other values will result
# m = 100 # default mass in kg - arbitrarily set, as the other values will result

Tspec = [0.0]
for T in arange(minT, maxT, deltaT):
    Tspec.append(T)



def runGM(gmNum, avList, avListX, avListY):
    tempStartTime = time()

    gmID = gmNum
    accGMfileFP = open(path + list_of_GMs[gmID] + "_FP.at2", "r")
    accGMfileFN = open(path + list_of_GMs[gmID] + "_FN.at2", "r")

    # Ignore the first line for FN file and two lines of FP file
    next(accGMfileFP)
    next(accGMfileFN)
    next(accGMfileFN)

    # get time interval which corresponds to the common dt of the GM records
    strWithDT = accGMfileFP.readline()
    print(strWithDT)
    strWithDT = strWithDT[:-4]
    while strWithDT[1] != '.':
        strWithDT = strWithDT[1:]
    dt = float(strWithDT)

    # multiplier for the ground motion numbers (e.g. 9.81 if the units are g)
    multiplier = 9.81 * scaleFactorsGM[gmID]

    # acceleration records are imported as lists
    accX = []
    accY = []

    for ax in accGMfileFP:
        accX.append(float(ax) * multiplier)

    for ay in accGMfileFN:
        accY.append(float(ay) * multiplier)

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

    print("Time to calculate {}: {:.2f}".format(list_of_GMs[gmNum], tempTime))

    fig, axs = plt.subplots(2, facecolor=(.2, .2, .2))
    axs[0].plot(Tspec, Aspec, label='Resultant', linewidth=2, color='yellow')
    axs[0].plot(Tspec, Axspec, label='Spectral-X', linewidth=1, color='red')
    axs[0].plot(Tspec, Ayspec, label='Spectral-Y', linewidth=1, color='green')
    axs[0].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[1].plot(accX, label='acceleration-X', linewidth=0.5, color='red')
    axs[1].plot(accY, label='acceleration-Y', linewidth=0.5, color='green')
    axs[0].set_facecolor('black')
    axs[1].set_facecolor('black')
    fig.suptitle(list_of_GMs[gmID], fontsize=12)
    plt.show()
    #plt.savefig(path + 'output\\' + list_of_GMs[gmID] + '.png')

    return tempTime


def spectral_Shape_factor_class_C(periods):
    Ch = []
    for Ti in periods:
        if Ti == 0:
            Ch.append(Z * 1.33)
        if Ti < 0.1:
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
        if Ti < 0.1:
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
    axs[0].plot(Tspec, averageSpectrum, label='Resultant', linewidth=2, color='yellow')
    axs[0].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[1].plot(Tspec, averageSpectrum_X, label='Spectral-FP', linewidth=1, color='red')
    axs[1].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[2].plot(Tspec, averageSpectrum_Y, label='Spectral-FN', linewidth=1, color='green')
    axs[2].plot(Tspec, spectral_Shape_factor_interpolation_C_D(Tspec, 0.8), label='Design Spectrum', linewidth=1, color='blue')
    axs[0].set_facecolor('black')
    axs[1].set_facecolor('black')
    axs[2].set_facecolor('black')
    fig.suptitle("Accumulated Average", fontsize=12)
    plt.show()
    #plt.savefig(path + 'output\\' + "average-spectra" + '.png')

print("Total time: {:.2f}".format(time() - startTime))
