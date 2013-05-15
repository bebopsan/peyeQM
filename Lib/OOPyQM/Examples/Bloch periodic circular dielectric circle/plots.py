# -*- coding: utf-8 -*-
"""
Created on Sat May 11 19:00:07 2013
Makes plots for solutions of reduced Brilluoin cell perimeter, in the 
context of a photonic crystals simulation.
@author: santiago
"""
import pickle 
import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import sqrt, around, pi 

filename = 'unit_cell_1-02_two_reg'
a = 1 # Lenght of the k grid

f = open(filename +'.pkl', 'r')
[freq, k_div] = pickle.load(f)  # frequency

fig = plt.figure()
ax = fig.add_subplot(111)
facecolors = ( 'green', 'red','blue', 'black')  #Define 4 colors for gaps

plt.title(r"Dispersion relations for Square lattice where r/a =0.2")

for i in range(0,freq.shape[1],2):
    ax.plot(freq[:,i])
    if i < freq.shape[1]-2 and min(freq[:,i+2]) - max(freq[:,i]) > 0.0001:
        # Draw a rectangle to show where is the bandgap        
        rec = plt.axhspan(max(freq[:,i]), min(freq[:,i+2]), alpha = 0.5)
        # Calculate gap-midgap ratio        
        dw = min(freq[:,i+2]) - max(freq[:,i]) # gap_width
        w = max(freq[:,i]) + dw/2              # frequency at the middle
        gmr = dw/w * 100  # gap-midgap ratio
        rec.set_label(r"$\frac{\Delta\omega}{\omega_m} =%.2f$%%"%(gmr) )
        rec.set_facecolor(facecolors[i/2])

# Plot lines defining the three important ponts of the perimeter
l = plt.axvline(x=k_div)
l = plt.axvline(x=k_div*2+1)
l = plt.axvline(x= freq.shape[0])
# Set figure handles
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles[::-1], labels[::-1], loc=1,bbox_to_anchor=(1, 0.8))
fig.subplots_adjust(left = 0.1, right=0.76)
# Define limits and letters for the corner points
plt.ylim(0, max(freq[:,-1]))
plt.xlim(0, i)
ax.set_xticks((0,k_div,k_div*2+1, freq.shape[0]))
ax.set_xticklabels(('$\Gamma$','$\chi$','$M$','$\Gamma$'),fontsize = 20)
ax.set_ylabel(r'Frequency $ \frac{\omega a}{2\pi c} $', fontsize = 20)
#Sa
plt.savefig("img/"+ filename +".svg")
plt.show()
f.close

