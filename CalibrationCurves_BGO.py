# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:08:40 2024

@author: adria
"""

#IMPORTS
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo


############################~~~Definitions~~~##################################

'''
Function that skips text and stores two columns of data into a dictionary,
Given a file name and the line of text before the data
Returns a dictionary
'''
def process_file(filename, first_line, final_line):
    
    #Initialise Dictionary and counter
    Flux = []
    count = 0
    
    #Open and set the file to be Read
    with open(filename, 'r') as file:
        
        #Begin loop to cycle through the rows
        for line in file.readlines():
            if(line.startswith(final_line)):
                count = 2
               
            #If statement to determine when to begin storing the data
            if(count == 1 ):#Begin storing data
                
     #Tell the program to seperate the values and record the columns seperately
     #strip() removes any whitespaces and split() seperates the values
     #map the values as floats and append them to the list in the dictionary
                
                flux = float(line.strip())
                Flux.append((flux))
            
            #Searching for the final text line to update the counter
            elif line.startswith(first_line):
                count = 1
               
    return Flux


'''
A Function that shortens 2 arrays to data regions of interest
Given a boundary of data (ROI, 2 values stored in an array)
Returns an array of all the data between the boundaries
'''
def ROI_Maker(x,y,ROI):
    
    
    ROI_x = []
    ROI_y = []
    for i in range(len(x)):
        
        if(x[i] >= ROI[0] and x[i] <= ROI[1]):
            
            ROI_x.append(x[i])
            ROI_y.append(y[i])
        elif(x[i] >= ROI[1]):
            break
            
    return np.asarray(ROI_x), np.asarray(ROI_y)

'''
Function that returns a Gaussian line
Used for fitting the Gaussian to the photopeak
Given x values, Amplitude, centroid Wavelength, sigma, line function
'''
def Gaussinator(x,A,mu,std):
    
    Gauss = (A/(np.sqrt(2*np.pi)*std))*np.exp(-((x - mu)**2)/(2*std**2)) 
            
    return Gauss

    

'''
Function that returns an array of estimates for Gaussian parameteres
Returns the average value of the data 
Formula from the manual plus the line
Given; X and Y values, and Region of Interests
'''
def Paraminator(x, y,ROI):
    
    x, y = ROI_Maker(x, y, ROI)
    
    Sum=0
    #Finding the wavelength at which the maximum flux occurs
    for i in range(len(y)):
        Sum += y[i]
 
        if(y[i] == max(y)):
            pos = i
    
    #avg = Sum/len(y)
    
    #Finding the Amplitude
    A = max(y)
    
    #Setting the Wavelength of maximum Flux        
    mu = x[pos]
    
    #Calculating the Standard Deviation from the ROI flux data
    std = np.std(y)
   
    #std = np.std(y)
    #std = 0.8
    #Return values for graphing and printing, values rounded after calculation
    return A,mu , std#,avg

def resolution_model(E, a, b, c):
     #return a*(E*E) + (b*E) + c
     return 1/(a / np.sqrt(E) + b + c * np.sqrt(E))
 
    
def efficiency_model(E,a,b,c):
    return a + (b*np.log(E) + (c*(np.log(E))**2))

###############################################################################
###############################~~~MAIN~~~######################################


#########FILE READER##############
#Setting filenames and line boundaries of data in the text file
filename_AM241 = 'BGO AM241 0deg 300sec.spe'
filename_Ba133 = 'BGO Ba133 0deg 300sec.spe'
filename_Cs137 = 'BGO Cs137 0deg 300sec.spe'
filename_Co60 = 'BGO Co60 0deg 600sec.spe'
filename_Background = 'BGO Background.spe'
filename_Am241_Off90 = 'BGO AM241 90deg 300sec.spe'


first_line = '0 1023'
final_line = '$ROI:'
#Creating Arrays of data
Flux_Am241 = np.asarray(process_file(filename_AM241, first_line, final_line))
Flux_Ba133 = np.asarray(process_file(filename_Ba133, first_line, final_line))
Flux_Cs137 = np.asarray(process_file(filename_Cs137, first_line, final_line))
Flux_Co60 = np.asarray(process_file(filename_Co60, first_line, final_line))
Flux_Bg = np.asarray(process_file(filename_Background, first_line, final_line))
Flux_Am241_Off90 = np.asarray(process_file(filename_Am241_Off90, first_line, final_line))


#Scaling Background to Co60
ScaledFlux_Bg = 2*Flux_Bg
#Subtracting the Background from Co60
ReducedFlux_Co60 = np.zeros_like(Flux_Co60)

for i in range(len(Flux_Co60)):

    #To avoid negative flux values
    if(Flux_Co60[i] < ScaledFlux_Bg[i]):
        ReducedFlux_Co60[i] = 0 
    else:
        ReducedFlux_Co60[i] = Flux_Co60[i] - ScaledFlux_Bg[i]




#Creating Channels
Channels= np.zeros(len(Flux_Am241))

for i in range(len(Flux_Am241) - 1):
    
    Channels[i] = i



######Calibration Curve#########
#Expected Energies of PhotoPeaks, in order of expected occurence
#Some peaks are assumed to have merged and the two values are averaged
ExpectedEvPeaks = np.array([59.54, 80.31,661.657, 1173.228,1332.49])



##########GAUSSIAN 230 400
#Setting Region of Interest
ROI_Cs = (290,390)
ROI_Am = [15,45]
ROI_Ba = [(25,55),(150,210)]
ROI_Co = [(550,630),(640,740)]

Flux_Co1 = []
Flux_Co2 = []

#Due to Gaussian fitting problems, the arrays were adjusted to aid the functions
#Data that was affecting the peaks was removed 
for i in range(len(Flux_Co60)):
    if(i < 630):    
        Flux_Co1.append(ReducedFlux_Co60[i])
        Flux_Co2.append(0)
    else:    
        Flux_Co1.append(0)
        Flux_Co2.append(ReducedFlux_Co60[i])
Flux_Co1 = np.array(Flux_Co1)
Flux_Co2 = np.array(Flux_Co2)

#CREATING BOUNDS
_b0_Cs = ([0,300,0],[np.inf,np.inf,np.inf])
_b0_Co = ([0,595,0],[np.inf,np.inf,np.inf])
_b0_Co2 = ([0,630,0],[np.inf,np.inf,np.inf])
_b0_Am = ([0,20,0],[np.inf,40,np.inf])
_b0_Ba = ([0,35,0],[12487,60,np.inf])
_b0_Ba2 = ([0,160,0],[38000,200,np.inf])

#Using the paraminator for gaussian estimations
_p0_Cs = Paraminator(Channels,Flux_Cs137,ROI_Cs)
_p0_Ba = Paraminator(Channels,Flux_Ba133,ROI_Ba[0])
_p0_Ba2 = Paraminator(Channels,Flux_Ba133,ROI_Ba[1])
_p0_Am = Paraminator(Channels,Flux_Am241,ROI_Am)
_p0_Co = Paraminator(Channels, Flux_Co1,ROI_Co[0])
_p0_Co2 = Paraminator(Channels,Flux_Co2,ROI_Co[1])
_p0_Am90 = Paraminator(Channels,Flux_Am241_Off90,ROI_Am)


#Finding true std,mu and A
popt_Am90, pcov_Am90 = spo.curve_fit(Gaussinator,Channels,Flux_Am241_Off90, p0=_p0_Am,bounds = _b0_Am)
popt_Cs137, pcov_Cs137 = spo.curve_fit(Gaussinator,Channels,Flux_Cs137, p0=_p0_Cs,bounds = _b0_Cs)
popt_Ba, pcov_Ba = spo.curve_fit(Gaussinator,Channels,Flux_Ba133, p0=_p0_Ba,bounds = _b0_Ba)
popt_Ba2, pcov_Ba2 = spo.curve_fit(Gaussinator,Channels,Flux_Ba133, p0=_p0_Ba2,bounds = _b0_Ba2)
popt_Am, pcov_Am = spo.curve_fit(Gaussinator,Channels,Flux_Am241, p0=_p0_Am,bounds = _b0_Am)
popt_Co, pcov_Co = spo.curve_fit(Gaussinator,Channels,Flux_Co1, p0=_p0_Co,bounds = _b0_Co)
popt_Co2, pcov_Co2 = spo.curve_fit(Gaussinator,Channels,Flux_Co2, p0=_p0_Co2,bounds = _b0_Co2)

#Creating gaussian Fit
yfit_Cs137 = Gaussinator(Channels, *popt_Cs137)
yfit_Ba = Gaussinator(Channels, *popt_Ba)
yfit_Ba2 = Gaussinator(Channels, *popt_Ba2)
yfit_Am = Gaussinator(Channels, *popt_Am)
yfit_Co = Gaussinator(Channels, *popt_Co)
yfit_Co2 = Gaussinator(Channels, *popt_Co2)


###############################################################################
#############################~~CALCULATIONS~~~#################################


#Scaling the Channels to Energies
#
Channel_Scale_Cs = ExpectedEvPeaks [2]/popt_Cs137[1]
Channel_Scale_Am = ExpectedEvPeaks[0]/popt_Am[1]
Channel_Scale_Ba = ExpectedEvPeaks[1]/popt_Ba[1]
#Channel_Scale_Ba2 = Exp_Peaks_Ba133[5]/popt_Ba2[1]
Channel_Scale_Co = ExpectedEvPeaks[3]/popt_Co[1]
Channel_Scale_Co2 = ExpectedEvPeaks[4]/popt_Co2[1]

#Averaging Value
Channel_Scale = (Channel_Scale_Am + Channel_Scale_Ba + Channel_Scale_Co
        + Channel_Scale_Cs + Channel_Scale_Co2)/5

#Creating Energy axis
eV_Channels = np.zeros_like(Channels)
for i in range(len(Channels)):
    
    eV_Channels[i] = Channel_Scale*Channels[i]






###############################################################################
###############################~~PLOTTING~~~###################################


#Plotting Spectra and Gaussians
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Flux [counts] Vs Channels')
ax.set(xlabel = 'Channels', ylabel = 'Flux [counts]')

#PLotting Gaussians
ax.plot(Channels,yfit_Cs137)
ax.plot(Channels,yfit_Ba)
ax.plot(Channels,yfit_Ba2)
ax.plot(Channels,yfit_Am)
ax.plot(Channels,yfit_Co)
ax.plot(Channels,yfit_Co2)


ax.plot(Channels,Flux_Am241, color = 'm', marker = 'd',
           label = 'Am241 counts')

ax.plot(Channels,Flux_Ba133, color = 'red', marker = 'd',
           label = 'Ba133 counts')

ax.plot(Channels,Flux_Cs137, color = 'blue', marker = 'd',
           label = 'Cs137 counts')

ax.plot(Channels,Flux_Co60, color = 'green', marker = 'd',
           label = 'Co60 counts')
ax.plot(Channels,Flux_Bg, color = 'black', marker = 'd',
           label = 'Background counts')

ax.legend(loc = 'best', fontsize = 10)
plt.show()


#Plotting Energies
fig, ax = plt.subplots(figsize = (15,10))


ax.set_title('Flux [counts] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Flux [counts]')

ax.plot(eV_Channels,Flux_Am241, color = 'm', marker = 'd',
           label = 'Am241 counts')

ax.plot(eV_Channels,Flux_Ba133, color = 'red', marker = 'd',
           label = 'Ba133 counts')

ax.plot(eV_Channels,Flux_Cs137, color = 'blue', marker = 'd',
           label = 'Cs137 counts')

ax.plot(eV_Channels,Flux_Co60, color = 'green', marker = 'd',
           label = 'Co60 counts')
ax.plot(eV_Channels,Flux_Bg, color = 'black', marker = 'd',
           label = 'Background counts')

ax.legend(loc = 'best', fontsize = 10)
plt.show()


#Plotting Calibration
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Energy [keV] Vs Channels')
ax.set(xlabel = 'Channels', ylabel = 'Energy [keV] ')


ChannelPeaks = np.array([popt_Am[1],popt_Ba[1],popt_Cs137[1], popt_Co[1], popt_Co2[1]])


a,b,c = np.polyfit(ChannelPeaks, ExpectedEvPeaks, 2)

ax.scatter(ChannelPeaks, ExpectedEvPeaks)
plt.plot(ChannelPeaks, a*(ChannelPeaks*ChannelPeaks) + b*ChannelPeaks + c ,
         label = f'For a quadratic curve: a = {a:0.4f}, b = {b:0.4f} & c = {c:0.4f}')
ax.legend(loc = 'upper left', fontsize = 12)

#Plotting Cobalt with no Background
#Check if Bg has values greater than Co60

'''
for i in range(len(Flux_Co60)):
    
    if(Flux_Bg[i] > Flux_Co60[i]):
        print('Background has value greater than Co_60 at:')
        print(i)
'''

#Plotting Cobalt
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Flux [counts] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Flux [counts]')  
ax.plot(Channels, Flux_Bg, color = 'blue', label = 'Background Flux')
ax.plot(Channels,yfit_Co)
ax.plot(Channels,yfit_Co2)
ax.plot(Channels, ReducedFlux_Co60, color = 'black', label = 'Scaled Co-60 Counts')

ax.legend(loc = 'best', fontsize = 10)


####RESOLUTION###################################
#Resolution
'''
 
 R= δE/E
 where,
 E is the photopeak energy
 δE is the FWHM of the photopeak.

'''

#Finding change in E
δE_Am = popt_Am[2]*(2*np.sqrt(2*(np.log(2))))
δE_Cs = popt_Cs137[2]*(2*np.sqrt(2*(np.log(2))))
δE_Ba = popt_Ba[2]*(2*np.sqrt(2*(np.log(2))))
δE_Ba2 = popt_Ba2[2]*(2*np.sqrt(2*(np.log(2))))
δE_Co = popt_Co[2]*(2*np.sqrt(2*(np.log(2))))
δE_Co2 = popt_Co2[2]*(2*np.sqrt(2*(np.log(2))))


#Finding E
E_Am = popt_Am[1]*Channel_Scale
E_Cs = popt_Cs137[1]*Channel_Scale
E_Ba = popt_Ba[1]*Channel_Scale
E_Ba2 = popt_Ba2[1]*Channel_Scale
E_Co = popt_Co[1]*Channel_Scale
E_Co2 = popt_Co2[1]*Channel_Scale


#Calculating Resolution
R_Am = δE_Am/E_Am
R_Cs = δE_Cs/E_Cs
R_Ba = δE_Ba/E_Ba
R_Ba2 = δE_Ba2/E_Ba2
R_Co = δE_Co/E_Co
R_Co2 = δE_Co2/E_Co2


#Creating Arrays for curve fitting and plotting
R = np.array([R_Am,R_Cs,R_Co,R_Co2])
E = np.array([E_Am,E_Cs,E_Co,E_Co2])

#Curve Fitting
initial_guess = [2, 0.9, 0.09]
popt_R, pcov_R = spo.curve_fit(resolution_model, E, R, p0=initial_guess)
#Saving values to print
a,b,c = popt_R[0],popt_R[1],popt_R[2]


#Plotting Resolution Curve
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Resolution [%] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Resolution [%]')
ax.set(xscale = 'log', yscale = 'log')


E_fit = np.linspace(min(E), max(E), 100)  # finer range for smooth curve
R_fit = resolution_model(E_fit, *popt_R)
plt.plot(E_fit, R_fit, label=f'For a quadratic curve: a = {a:0.3f}, b = {b:0.3f} & c = {c:0.3f}', color="blue")
ax.scatter(E,R)
ax.legend(loc = 'best', fontsize = 20)


#Efficiencies
#Amplitudes
A_Cs = (popt_Cs137[0]/(np.sqrt(2*np.pi)*popt_Cs137[2]))
#Due to poor Gaussian fit, Ba-133 had to have amplitude hardcoded
A_Ba = 1750
A_Am90 = (popt_Am90[0]/(np.sqrt(2*np.pi)*popt_Am90[2]))
A_Am = (popt_Am[0]/(np.sqrt(2*np.pi)*popt_Am[2]))
A_Co = (popt_Co[0]/(np.sqrt(2*np.pi)*popt_Co[2]))
A_Co2 = (popt_Co2[0]/(np.sqrt(2*np.pi)*popt_Co2[2]))


#Absolute Efficiency
#Activity
N_Am = 409892
N_Ba = 19753
N_Cs = 160584
N_Co = 1033

Eff_Abs_Am = (A_Am/300)/N_Am
Eff_Abs_Ba = (A_Ba/300)/N_Ba
Eff_Abs_Cs = (A_Cs/300)/N_Cs
Eff_Abs_Co = (A_Co/600)/N_Co
Eff_Abs_Co2 = (A_Co2/600)/N_Co
Eff_Abs_Am90 = (A_Am90/300)/N_Am


#Intrinsic Efficiencies
SolidAngle = (((2.54E-2)**2)*np.pi)/(4*np.pi*(14.4E-2))
SolidAngle90 = ((2*2.54E-2)*5.08E-2)/(4*np.pi*(14.4E-2))

PhotonRate_Am = N_Am * 300 * (SolidAngle/(4*np.pi))
PhotonRate_Ba = N_Ba * 300 * (SolidAngle/(4*np.pi))
PhotonRate_Cs = N_Cs * 300 * (SolidAngle/(4*np.pi))
PhotonRate_Co = N_Co * 300 * (SolidAngle/(4*np.pi))
PhotonRate_Am90 = N_Am * 300 * (SolidAngle90/(4*np.pi))



Eff_Int_Am90 = (A_Am90/300)/PhotonRate_Am90
Eff_Int_Am = (A_Am/300)/PhotonRate_Am
Eff_Int_Ba = (A_Ba/300)/PhotonRate_Ba
Eff_Int_Cs = (A_Cs/300)/PhotonRate_Cs
Eff_Int_Co = (A_Co/300)/PhotonRate_Co
Eff_Int_Co2 = (A_Co2/300)/PhotonRate_Co
Eff = np.array([Eff_Int_Am,Eff_Int_Ba,Eff_Int_Cs,Eff_Int_Co,Eff_Int_Co2])


#New E array
E = np.array([E_Am,E_Ba,E_Cs,E_Co,E_Co2])


#Curve fitting for Efficiency
popt_Eff, pcov_Eff = spo.curve_fit(efficiency_model, E, Eff)
a,b,c = popt_Eff[0],popt_Eff[1],popt_Eff[2]


#Plotting Efficiency
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Efficiency [%] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Efficiency [%]')
ax.set(xscale = 'log', yscale = 'log')


E_fit = np.linspace(min(E), max(E), 100)  # finer range for smooth curve
Eff_fit = efficiency_model(E_fit, *popt_Eff)
plt.plot(E_fit, Eff_fit, label=f'For a quadratic curve: a = {a:0.3E}, b = {b:0.3E} & c = {c:0.3E}', color="blue")
ax.scatter(E,Eff)
ax.legend(loc = 'best', fontsize = 20)





