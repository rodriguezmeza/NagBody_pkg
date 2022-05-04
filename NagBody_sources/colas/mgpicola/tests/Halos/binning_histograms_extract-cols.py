
# source activate nbodykit-env

#========= Importing libs =========

from nbodykit.lab import *

import os
import matplotlib.pyplot as plt

import numpy as np
import math
import statistics

#========= Reading model files =========

# Setting variables and parameters:
#sim_name = "F6_Box1_27"

file_name = "Rockstar_M200c_F6_B1_B1024_NP1024_S27.dat"
prefix_name = "../../../../halo_catalogues/"
sim_name = "F6"
box_name = "box1"
snap_name = "27"

sim_name_save = sim_name+"_"+box_name+"_"+snap_name

#box_name = "../../../../halo_catalogues/F6/box1/27/Rockstar_M200c_F6_B1_B1024_NP1024_S27.dat"

file_path_name = prefix_name+sim_name+"/"+box_name+"/"+snap_name+"/"+file_name

save_name = "./"
LBox = 1024.0
Nparticles = 1024*1024*1024
print(Nparticles)

# Reading files
inputFile = np.loadtxt(file_path_name)

massvec=inputFile[:,2]
npvec=inputFile[:,7]
xvec=inputFile[:,8]
yvec=inputFile[:,9]
zvec=inputFile[:,10]
nval=massvec.shape[0]
npvecsum = np.sum(npvec)
print(nval)
print(npvecsum)
print(npvecsum/Nparticles)

for i in range(0,nval):
    xvec[i] = xvec[i]/LBox
    yvec[i] = yvec[i]/LBox
    zvec[i] = zvec[i]/LBox

massmax = np.amax(massvec)
massmin = np.amin(massvec)
#print('mass minimum:',massmin)
#print('mass maximum:',massmax)

print('mass minimum:','%e' % massmin)
print('mass maximum:','%e' % massmax)

npmin = np.min(npvec)
npmax = np.max(npvec)
print('np minimum:','%e' % npmin)
print('np maximum:','%e' % npmax)

yf = 15.5
yi = 10.5
Ny = 30
dy = (yf-yi)/int(Ny)
yarray = np.arange(yi,yf,dy)
xarray = np.arange(yi,yf,dy)

print(dy)
print(yarray)

for i in range(0,Ny):
    xarray[i] = np.power(10.0,yarray[i])

print(xarray)


#massbins = [1.0e11,1.0e12,1.0e13,1.0e14,1.0e15]
#massbins = [0.5e11,1.0e11,0.5e12,1.0e12,0.5e13,1.0e13,0.5e14,1.0e14,0.5e15,1.0e15]
massbins = xarray

masshisto = np.histogram(massvec,bins=massbins)
print(masshisto)

fig1 = plt.figure()

#plt.hist(massvec, bins='auto')
plt.hist(massvec, bins=xarray)
#plt.hist(masshisto)


plt.xlabel(r"$M$ [$M_{\odot}$]")
plt.ylabel(r"$N$")
#plt.xlim(0.005, 0.6)
#plt.xlim(0.01, 10.0)
#plt.ylim(1,100000)
plt.xscale('log')
#plt.yscale('log')

fig_txt_name = save_name+'halos_mass_hist_'+sim_name_save+'.pdf'

plt.savefig(fig_txt_name)

# Saving to file m x y z
# (1) (2) (3) (4)
# x y z mass
arraya=np.array([xvec,yvec,zvec,massvec])
arrayb=np.transpose(arraya)
print arrayb

array_txt_name = save_name+'halos_mxyz_'+sim_name_save+'.dat'

nvalstr = str(nval)
np.savetxt(array_txt_name,arrayb,delimiter='\t',header=nvalstr,newline='\r\n')
