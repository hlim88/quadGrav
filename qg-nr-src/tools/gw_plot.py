###########################################################################
#
# May.21.2018
# Python script to plot GW waveform. Read the result data from Dendro 
# simulation and save a plot as pdf file. 
#
# Usage : pythong gw_plot.py <your data> <max_depth> <grid_size>
# Black solid line : Real part of (2,2) mode Psi4
# Red dot-dash line : Imaginary part of (2,2) mode Psi4
#
###########################################################################

#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

data = sys.argv[1]
max_depth = float(sys.argv[2])
grid_size = float(sys.argv[3])

# Calculate dt
cfl = 0.1
dx = grid_size/2**max_depth
dt = cfl*dx

xbh1, ybh1 = np.loadtxt(data, usecols=(0,2), unpack = True)
xbh2, ybh2 = np.loadtxt(data, usecols=(0,3), unpack = True)

fig = plt.figure(figsize=(10,6), dpi=100)
ax = fig.add_subplot(111)

#ax.axis([0,450,-0.001,0.001])
#ax.set_xlabel('x/M')
ax.set_ylabel('r$\Psi_{2,2}$')
ax.plot(xbh1*dt, ybh1, linestyle ="-", color ="black", label='Re[$\Psi_{2,2}$]')
ax.plot(xbh2*dt, ybh2, linestyle ="-.", color ="red", label='Im[$\Psi_{2,2}$]')
ax.legend(loc='upper left')

print("Plot GW extraction from spin-weighted spherical harmonics integration")
print("Your dt= %s" % dt)
plt.savefig('waveform.pdf')

