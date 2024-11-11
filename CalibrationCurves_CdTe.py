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
     
     return 1/(a / np.sqrt(E) + b + c * np.sqrt(E))

def efficiency_model(E,a,b,c):
    return a + (b*np.log(E) + (c*(np.log(E))**2))


###############################################################################
###############################~~~MAIN~~~######################################


#########FILE READER##############
#Setting filenames and line boundaries of data in the text file
filename_AM241 = 'CdTe AM241 0deg 600sec.mca'
filename_Ba133 = 'CdTe Ba133 0deg 600sec.mca'
filename_Cs137 = 'CdTe Cs137 0deg 600sec.mca'
filename_Co60 = 'CdTe Co60 0deg 1200sec.mca'
filename_Background = 'CdTe Background.mca'
filename_AM241_sample = 'SampleAM241CdTeDATA.mca'
filename_CdTeAm241_Off30 = 'CdTe AM241 30deg 600sec.mca'
filename_CdTeAm241_Off60 = 'CdTe AM241 60deg 600sec.mca'
filename_CdTeAm241_Off90 = 'CdTe AM241 90deg 600sec.mca'

first_line = '<<DATA>>'
final_line = '<<END>>'
#Creating Arrays of data
Flux_Am241 = np.asarray(process_file(filename_AM241, first_line, final_line))
Flux_Ba133 = np.asarray(process_file(filename_Ba133, first_line, final_line))
Flux_Cs137 = np.asarray(process_file(filename_Cs137, first_line, final_line))
Flux_Co60 = np.asarray(process_file(filename_Co60, first_line, final_line))
Flux_Bg = np.asarray(process_file(filename_Background, first_line, final_line))
Flux_Am241_sample = np.asarray(process_file(filename_AM241_sample, first_line, final_line))
Flux_CdTeAM241_Off30 = np.asarray(process_file(filename_CdTeAm241_Off30, first_line, final_line))
Flux_CdTeAM241_Off60 = np.asarray(process_file(filename_CdTeAm241_Off60, first_line, final_line))
Flux_CdTeAM241_Off90 = np.asarray(process_file(filename_CdTeAm241_Off90, first_line, final_line))

#Creating Channels
Channels= np.zeros(len(Flux_Am241))

for i in range(len(Flux_Am241) - 1):
    
    Channels[i] = i



######Calibration Curve#########
#Expected Energies of PhotoPeaks, in order of expected occurence
#Some peaks are assumed to have managed and the two values are averaged
Exp_Peaks_Am241 = [26.34,33.20,59.54]



##########GAUSSIAN 230 400
#Setting Region of Interest
ROI_Am = [(265,290),(340,370),(1165,1200)]




#CREATING BOUNDS
_b0_Am = ([0,260,0],[np.inf,290,np.inf])
_b0_Am2 = ([0,330,0],[np.inf,370,np.inf])
_b0_Am3 = ([0,1165,0],[np.inf,1200,np.inf])


#Using the paraminator for gaussian estimations
_p0_Am3 = Paraminator(Channels,Flux_Am241,ROI_Am[2])
_p0_Am = Paraminator(Channels,Flux_Am241,ROI_Am[0])
_p0_Am2 = Paraminator(Channels,Flux_Am241,ROI_Am[1])
_p0_Am_90 = Paraminator(Channels,Flux_CdTeAM241_Off90,ROI_Am[0])
_p0_Am2_90 = Paraminator(Channels,Flux_CdTeAM241_Off90,ROI_Am[1])
_p0_Am3_90 = Paraminator(Channels,Flux_CdTeAM241_Off90,ROI_Am[2])



#Curve Fitting
popt_Am3, pcov_Am3 = spo.curve_fit(Gaussinator,Channels,Flux_Am241, p0=_p0_Am3,bounds = _b0_Am3)
popt_Am, pcov_Am = spo.curve_fit(Gaussinator,Channels,Flux_Am241, p0=_p0_Am,bounds = _b0_Am)
popt_Am_90, pcov_Am_90 = spo.curve_fit(Gaussinator,Channels,Flux_CdTeAM241_Off90, p0=_p0_Am_90,bounds = _b0_Am)
popt_Am2_90, pcov_Am_90 = spo.curve_fit(Gaussinator,Channels,Flux_CdTeAM241_Off90, p0=_p0_Am2_90,bounds = _b0_Am2)
popt_Am3_90, pcov_Am_90 = spo.curve_fit(Gaussinator,Channels,Flux_CdTeAM241_Off90, p0=_p0_Am3_90,bounds = _b0_Am3)
popt_Am2, pcov_Am2 = spo.curve_fit(Gaussinator,Channels,Flux_Am241, p0=_p0_Am2,bounds = _b0_Am2)


#Creating Gaussian Fits
yfit_Am3 = Gaussinator(Channels, *popt_Am3)
yfit_Am = Gaussinator(Channels, *popt_Am)
yfit_Am2 = Gaussinator(Channels, *popt_Am2)


###############################################################################
#############################~~CALCULATIONS~~~#################################


#Scaling the Channels to Energies

Channel_Scale_Am = Exp_Peaks_Am241[0]/popt_Am[1]
Channel_Scale_Am2 = Exp_Peaks_Am241[1]/popt_Am2[1]
Channel_Scale_Am3 = Exp_Peaks_Am241[2]/popt_Am3[1]

#Averaging Values
Channel_Scale = (Channel_Scale_Am + Channel_Scale_Am3 + Channel_Scale_Am2)/3

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

#PLotting Gaussian Sadness
ax.plot(Channels,yfit_Am)
ax.plot(Channels,yfit_Am2)
ax.plot(Channels,yfit_Am3)


#Plotting Spectra against Channels
ax.plot(Channels,Flux_Am241, color = 'm', marker = 'd',
           label = 'Am241 counts')

#ax.plot(Channels,Flux_Am241_sample, color = 'black', marker = 'd',
#           label = 'Am241 counts')

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
ax.set_title('Channels Vs Energy [keV]')
ax.set(xlabel = 'Channels', ylabel = 'Energy [keV]')
ChannelPeaks = np.array([popt_Am[1],popt_Am2[1],popt_Am3[1]])


a,b ,c= np.polyfit(ChannelPeaks, Exp_Peaks_Am241 , 2)

ax.scatter(ChannelPeaks, Exp_Peaks_Am241)
plt.plot(ChannelPeaks, a*(ChannelPeaks**2) + b*ChannelPeaks + c,
         label = f'For a quadratic curve: a = {a:0.2E}, b = {b:0.4f} & c = {c:0.4f}')

ax.legend(loc = 'best', fontsize = 15)


#Resolution
'''
 
 R= δE/E
 where,
 E is the photopeak energy
 δE is the FWHM of the photopeak.

'''
#Finding Change in E
δE_Am = popt_Am[2]*(2*np.sqrt(2*(np.log(2))))
δE_Am2 = popt_Am2[2]*(2*np.sqrt(2*(np.log(2))))
δE_Am3 = popt_Am3[2]*(2*np.sqrt(2*(np.log(2))))

#FInding E
E_Am = popt_Am[1]*Channel_Scale
E_Am2 = popt_Am2[1]*Channel_Scale
E_Am3 = popt_Am3[1]*Channel_Scale

#Calculating Resolution
R_Am = δE_Am/E_Am
R_Am2 = δE_Am2/E_Am2
R_Am3 = δE_Am3/E_Am3


#Creating arrays for fitting and plotting
R = np.array([R_Am,R_Am2,R_Am3])
E = np.array([E_Am,E_Am2,E_Am3])

#Fitting Resolution Model
initial_guess = [2, 0.9, 0.09]
popt_R, pcov_R = spo.curve_fit(resolution_model, E, R, p0=initial_guess)
a,b,c = popt_R[0],popt_R[1],popt_R[2]


#PRinting Resolution
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Resolution [%] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Resolution [%]')
ax.set(xscale = 'log', yscale = 'log')


E_fit = np.linspace(min(E), max(E), 100)  # finer range for smooth curve
R_fit = resolution_model(E_fit, *popt_R)
plt.plot(E_fit, R_fit, label=f'For a quadratic curve: a = {a:0.3f}, b = {b:0.3f} & c = {c:0.3f}', color="blue")
ax.scatter(E,R)
ax.legend(loc = 'best', fontsize = 20)



####Efficiencies

#Amplitudes
A_Am = (popt_Am[0]/(np.sqrt(2*np.pi)*popt_Am[2]))
A_Am_90 = (popt_Am_90[0]/(np.sqrt(2*np.pi)*popt_Am_90[2]))
A_Am2 = (popt_Am2[0]/(np.sqrt(2*np.pi)*popt_Am2[2]))
A_Am2_90 = (popt_Am2_90[0]/(np.sqrt(2*np.pi)*popt_Am2_90[2]))
A_Am3 = (popt_Am3[0]/(np.sqrt(2*np.pi)*popt_Am3[2]))
A_Am3_90 = (popt_Am3_90[0]/(np.sqrt(2*np.pi)*popt_Am3_90[2]))



#Absolute Efficiency
#Activity
N_Am = 409892


Eff_Abs_Am = (A_Am/600)/N_Am
Eff_Abs_Am_90 = (A_Am_90/600)/N_Am
Eff_Abs_Am2 = (A_Am2/600)/N_Am
Eff_Abs_Am2_90 = (A_Am2_90/600)/N_Am
Eff_Abs_Am3 = (A_Am3/600)/N_Am
Eff_Abs_Am3_90 = (A_Am3_90/600)/N_Am



#Intrinsic Efficiencies

SolidAngle = (((5E-3)**2))/(4*np.pi*(14.4E-2))
SolidAngle90 = (((5E-3)*1E-3))/(4*np.pi*(14.4E-2))

PhotonRate_Am = N_Am * 300 * (SolidAngle/(4*np.pi))
PhotonRate_Am90 =  N_Am * 300 * (SolidAngle90/(4*np.pi))




Eff_Int_Am = (A_Am/300)/PhotonRate_Am
Eff_Int_Am2 = (A_Am2/300)/PhotonRate_Am
Eff_Int_Am3 = (A_Am3/300)/PhotonRate_Am
Eff_Int_Am_90 = (A_Am_90/300)/PhotonRate_Am90
Eff_Int_Am2_90 = (A_Am2_90/300)/PhotonRate_Am90
Eff_Int_Am3_90 = (A_Am3_90/300)/PhotonRate_Am90


#Creating arrays for fitting and plotting
Eff = np.array([Eff_Int_Am,Eff_Int_Am2,Eff_Int_Am3])
E = np.array([E_Am,E_Am2,E_Am3])

popt_Eff, pcov_Eff = spo.curve_fit(efficiency_model, E, Eff)
a,b,c = popt_Eff[0],popt_Eff[1],popt_Eff[2]


#Plotting
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title('Efficiency [%] Vs Energy [keV]')
ax.set(xlabel = 'Energy [keV]', ylabel = 'Efficiency [%]')
ax.set(xscale = 'log', yscale = 'log')


E_fit = np.linspace(min(E), max(E), 100)  # finer range for smooth curve
Eff_fit = efficiency_model(E_fit, *popt_Eff)
plt.plot(E_fit, Eff_fit, label=f'For a quadratic curve: a = {a:0.3E}, b = {b:0.3E} & c = {c:0.3E}', color="blue")
ax.scatter(E,Eff)
ax.legend(loc = 'best', fontsize = 17)