
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




###############################################################################
###############################~~~MAIN~~~######################################


#########FILE READER##############
#Setting filenames and line boundaries of data in the text file
filename_BGOAm241_Off0 = 'BGO AM241 0deg 300sec.spe'
filename_BGOAm241_Off30 = 'BGO AM241 30deg 300sec.spe'
filename_BGOAm241_Off60 = 'BGO AM241 60deg 300sec.spe'
filename_BGOAm241_Off90 = 'BGO AM241 90deg 300sec.spe'

filename_NaIAm241_Off0 = 'NaI AM241 0deg 300sec.spe'
filename_NaIAm241_Off30 = 'NaI AM241 30deg 300sec.spe'
filename_NaIAm241_Off60 = 'NaI AM241 60deg 300sec.spe'
filename_NaIAm241_Off90 = 'NaI AM241 90deg 300sec.spe'

filename_CdTeAm241_Off0 = 'CdTe AM241 0deg 600sec.mca'
filename_CdTeAm241_Off30 = 'CdTe AM241 30deg 600sec.mca'
filename_CdTeAm241_Off60 = 'CdTe AM241 60deg 600sec.mca'
filename_CdTeAm241_Off90 = 'CdTe AM241 90deg 600sec.mca'

first_line = '0 1023'

final_line = '$ROI:'
first_line_CdTe = '<<DATA>>'
final_line_CdTe = '<<END>>'

#Creating Arrays of data
Flux_BGOAM241_Off0 = np.asarray(process_file(filename_BGOAm241_Off0, first_line, final_line))
Flux_BGOAM241_Off30 = np.asarray(process_file(filename_BGOAm241_Off30, first_line, final_line))
Flux_BGOAM241_Off60 = np.asarray(process_file(filename_BGOAm241_Off60, first_line, final_line))
Flux_BGOAM241_Off90 = np.asarray(process_file(filename_BGOAm241_Off90, first_line, final_line))

Flux_NaIAM241_Off0 = np.asarray(process_file(filename_NaIAm241_Off0, first_line, final_line))
Flux_NaIAM241_Off30 = np.asarray(process_file(filename_NaIAm241_Off30, first_line, final_line))
Flux_NaIAM241_Off60 = np.asarray(process_file(filename_NaIAm241_Off60, first_line, final_line))
Flux_NaIAM241_Off90 = np.asarray(process_file(filename_NaIAm241_Off90, first_line, final_line))

Flux_CdTeAM241_Off0 = np.asarray(process_file(filename_CdTeAm241_Off0, first_line_CdTe, final_line_CdTe))
Flux_CdTeAM241_Off30 = np.asarray(process_file(filename_CdTeAm241_Off30, first_line_CdTe, final_line_CdTe))
Flux_CdTeAM241_Off60 = np.asarray(process_file(filename_CdTeAm241_Off60, first_line_CdTe, final_line_CdTe))
Flux_CdTeAM241_Off90 = np.asarray(process_file(filename_CdTeAm241_Off90, first_line_CdTe, final_line_CdTe))



#Creating Channels
Channels= np.zeros(len(Flux_BGOAM241_Off0))
Channels2 = np.zeros(len(Flux_CdTeAM241_Off0))

for i in range(len(Flux_BGOAM241_Off0) - 1):
    
    Channels[i] = i
    Channels2[i] = i


###############################################################################
###############################~~PLOTTING~~~###################################


#Plotting off-Axis Spectra

##Plotting BGO
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title ( 'BGO Flux [counts] Vs Channel')


ax.plot(Channels,Flux_BGOAM241_Off0, color = 'm', marker = 'd',
           label = 'Am241 0 degrees')
ax.plot(Channels,Flux_BGOAM241_Off30, color = 'r', marker = 'd',
           label = 'Am241 30 degrees')
ax.plot(Channels,Flux_BGOAM241_Off60, color = 'blue', marker = 'd',
           label = 'Am241 60 degrees')
ax.plot(Channels,Flux_BGOAM241_Off90, color = 'black', marker = 'd',
           label = 'Am241 90 degrees')

ax.set(xlabel = 'Channel', ylabel = 'Flux [counts]')
ax.legend(loc = 'best', fontsize = 10)


##Plotting NaI
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title ( 'NaI Flux [counts] Vs Channel')


ax.plot(Channels,Flux_NaIAM241_Off0, color = 'm', marker = 'd',
           label = 'Am241 0 degrees')
ax.plot(Channels,Flux_NaIAM241_Off30, color = 'r', marker = 'd',
           label = 'Am241 30 degrees')
ax.plot(Channels,Flux_NaIAM241_Off60, color = 'blue', marker = 'd',
           label = 'Am241 60 degrees')
ax.plot(Channels,Flux_NaIAM241_Off90, color = 'black', marker = 'd',
           label = 'Am241 90 degrees')

ax.set(xlabel = 'Channel', ylabel = 'Flux [counts]')
ax.legend(loc = 'best', fontsize = 10)


#Plotting CdTe
fig, ax = plt.subplots(figsize = (15,10))
ax.set_title ( 'CdTe Flux [counts] Vs Channel')


ax.plot(Channels2,Flux_CdTeAM241_Off0, color = 'm', marker = 'd',
           label = 'Am241 0 degrees')
ax.plot(Channels2,Flux_CdTeAM241_Off30, color = 'r', marker = 'd',
           label = 'Am241 30 degrees')
ax.plot(Channels2,Flux_CdTeAM241_Off60, color = 'blue', marker = 'd',
           label = 'Am241 60 degrees')
ax.plot(Channels2,Flux_CdTeAM241_Off90, color = 'black', marker = 'd',
           label = 'Am241 90 degrees')

ax.set(xlabel = 'Channel', ylabel = 'FLux [counts]')
ax.legend(loc = 'best', fontsize = 10)

