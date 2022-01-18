# Scale ground motions to soil-C according to the methodology of NZS1170.5:2004

import scaleGMtoNZS
from time import time
from os import path
import pickle


startTime = time()

# Set paths here
path_Linux = "/home/user/YourPathGoesHere/"
path_Windows = "C:\\YourPathGoesHere\\"

if path.exists(path_Linux):
    pathMain = path_Linux
elif path.exists(path_Windows):
    pathMain = path_Windows


# Set of periods to be checked
T1 = 2.79 # period of main translational mode in sec
minT = 0.4 * T1
maxT = 1.3 * T1 
periodIncrement = 1.02 # see NZS1170.5:2004, cl.5.5.2,(c),(v) "...each period within 10% of the preceeding one..."
print(f"T = {T1:1.2f} \t Tmin = {minT:1.2f} \t Tmax = {maxT:1.2f}")

# Define a list of the full set of periods to be checked
periodsSDOF = [minT]
while periodsSDOF[-1]<maxT:
    periodsSDOF.append( periodsSDOF[-1] * periodIncrement )

list_of_GM1 = []
list_of_GM2 = []

# Create list with all the filenames of the ground motions
for i in range(1,15):
    stGM = "773-WLGGE285842_TH"
    if i<10:
        stGM += "0"
        ext = ".txt"
    else:
        ext = "_bc.txt"
    
    stGM += str(i) + "_H"

    list_of_GM1.append( stGM + "1" + ext)
    list_of_GM2.append( stGM + "2" + ext)

k1_list1 = []
k1_list2 = []

targetSpec = scaleGMtoNZS.spectrum_class_CDint(periodsSDOF, Z=0.85*0.4*9.81, R=1, Tsite=0.7)
"""print("Periods:")
for n in periodsSDOF:
    print(f"{n:1.3f}")
print("Ch factors:")
for n in targetSpec:
    print(f"{n:1.3f}")"""

print(" Calculating K1 for X")
for f in list_of_GM1:
    #gmRec = scaleGMtoNZS.importGM(pathMain+f)
    resp_Spec = scaleGMtoNZS.respSpectrumGM(pathMain+f, periodsSDOF)
    #print("This is the unscaled spectrum of the ground motion: \n", resp_Spec)
    k1_list1.append( scaleGMtoNZS.factorK1(resp_Spec, targetSpec, periodsSDOF, T1) )
    pickle.dump( resp_Spec, open(pathMain+f+".shortSpectrum.dat", "wb") )

print("\n Calculating K1 for Y")
for f in list_of_GM2:
    #gmRec = scaleGMtoNZS.importGM(pathMain+f)
    resp_Spec = scaleGMtoNZS.respSpectrumGM(pathMain+f, periodsSDOF)
    k1_list2.append( scaleGMtoNZS.factorK1(resp_Spec, targetSpec, periodsSDOF, T1) )
    pickle.dump( resp_Spec, open(pathMain+f+".shortSpectrum.dat", "wb") )

#choose the minimum k1 for each ground motion
k1_list = []   
for n in range( len(k1_list1) ):
    k1_list.append( min(k1_list1[n] , k1_list2[n]) )

i=0
for n in k1_list:
    i += 1
    print(f"GM-{i} \t k1 = {n:.2f}")

print("Total time: {:.2f}".format(time() - startTime)) 

pickle.dump(k1_list, open(pathMain+"k1.dat", "wb"))
pickle.dump(list_of_GM1, open(pathMain+"GM1-List.dat", "wb"))
pickle.dump(list_of_GM2, open(pathMain+"GM2-List.dat", "wb"))
pickle.dump(periodsSDOF, open(pathMain+"periodsSDOF.dat", "wb"))


