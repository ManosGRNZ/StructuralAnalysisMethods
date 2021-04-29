import pickle
import matplotlib.pyplot as plt
from os import path
import tkinter 
import tksheet

path_Linux_DL = "/home/manos/Programs/python/structural/List-GM/output/"
path_Windows = "D:\\Programs\\python\\structural\\List-GM\\output\\"
path_Linux_LT = "/home/mb/Programs/python/structural/List-GM/output/"

pathMain = path_Linux_DL
if path.exists(path_Windows):
    pathMain = path_Windows
if path.exists(path_Linux_LT):
    pathMain = path_Linux_LT

list_of_GMs = ["RSN171", "RSN722", "RSN767", "RSN1084", "RSN1176", "RSN1481", "RSN1602", "RSN6960", "RSN8606", "RSN4032588", "RSN4032590" ]

#scaleFactorsGM = [1.00] * len(list_of_GMs)  # default value of 1 - to be changed
scaleFactorsGM = [0.93, 2.67, 0.90, 0.71, 1.72, 2.80, 0.77, 1.70, 1.24, 2.70, 2.00 ]

numGMsRun = len(list_of_GMs)
processGMList = list(range(numGMsRun))

Tspec = pickle.load(open(pathMain + "periods.data", 'rb'))

def spectral_Shape_factor_interpolation_C_D(periods=Tspec, Tsite=0.8, Z=0.4, R=1):
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


spectra = []

for gm in list_of_GMs:
    spectra.append(pickle.load(open(pathMain + gm + "_Spectrum_FP.data", 'rb')))


def plot_spectra(scaleFactorsStr):
    #scaleFactorsInput = [float(i) for i in scaleFactorsStr]
    scaleFactorsInput = [] 
    for f in scaleFactorsStr[0]:
        scaleFactorsInput.append(float(f))
    average = []
    for n in range(len(spectra[0])):
        s = 0
        i=0
        for a in spectra:
            s += a[n]*scaleFactorsInput[i]
            i += 1
        average.append(s/(numGMsRun*9.81))

    plt.plot(Tspec, average, label='Average', linewidth=2, color='red')
    plt.plot(Tspec, spectral_Shape_factor_interpolation_C_D(), label='Design Spectrum', linewidth=1, color='blue')
    #plt.set_facecolor('black')
    #plt.set_facecolor('black')
    plt.show()

def present(event):
    plot_spectra([ list(sheet.get_row_data(1) ) ])

#Create sheet GUI
top = tkinter.Tk()
top.geometry('900x150')
sheet = tksheet.Sheet(top, height=130, width=880)
sheet.grid()
scaleFactorsTemp = scaleFactorsGM
sheet.set_sheet_data([list_of_GMs,scaleFactorsTemp])
#scaleFactorsTemp = list(plot_spectra([ list( sheet.get_row_data(1) ) ] ) ) 
#sheet.enable_bindings( "end_edit_cell", list(plot_spectra([ float(n) for n in list(sheet.get_row_data(1) ) ] ) ) )
sheet.enable_bindings()
#sheet.extra_bindings("end_edit_cell", plot_spectra([ list(sheet.get_row_data(1) ) ] ) )
#sheet.mainloop()
buttonWindow = tkinter.Tk()
widget = tkinter.Button(buttonWindow, text='calculate')
widget.pack()
widget.bind('<Button-1>',  present)
widget.mainloop() 
