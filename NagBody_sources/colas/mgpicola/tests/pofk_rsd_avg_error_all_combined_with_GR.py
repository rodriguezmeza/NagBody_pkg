
# source activate nbodykit-env

#========= Importing libs =========
from nbodykit.lab import *

import os
import matplotlib.pyplot as plt

import numpy as np
import math
import statistics

#========= Reading model files =========
data_path_dir = 'F4/'

sim_combined_name = "F4_All_rel_GR_z1p000"

sim_name = "z1p000"
ext = ".dat"
model = "F4"

# Reading files
box1_name = data_path_dir+"box1001/pofk_rsd_mu_2d_l_"+model+"_Box1001_"+sim_name+ext
box2_name = data_path_dir+"box2001/pofk_rsd_mu_2d_l_"+model+"_Box2001_"+sim_name+ext
box3_name = data_path_dir+"box3001/pofk_rsd_mu_2d_l_"+model+"_Box3001_"+sim_name+ext
box4_name = data_path_dir+"box4001/pofk_rsd_mu_2d_l_"+model+"_Box4001_"+sim_name+ext
box5_name = data_path_dir+"box5001/pofk_rsd_mu_2d_l_"+model+"_Box5001_"+sim_name+ext

inputFile1 = np.loadtxt(box1_name)
inputFile2 = np.loadtxt(box2_name)
inputFile3 = np.loadtxt(box3_name)
inputFile4 = np.loadtxt(box4_name)
inputFile5 = np.loadtxt(box5_name)

F4_kval_1=inputFile1[:,0]
F4_psl0_1=inputFile1[:,1]
F4_psl2_1=inputFile1[:,2]
F4_psl4_1=inputFile1[:,3]
F4_psl02_1=inputFile1[:,4]
F4_nkval_1=F4_kval_1.shape[0]
print(F4_nkval_1)

F4_kval_2=inputFile2[:,0]
F4_psl0_2=inputFile2[:,1]
F4_psl2_2=inputFile2[:,2]
F4_psl4_2=inputFile2[:,3]
F4_psl02_2=inputFile2[:,4]
F4_nkval_2=F4_kval_2.shape[0]
print(F4_nkval_2)

F4_kval_3=inputFile3[:,0]
F4_psl0_3=inputFile3[:,1]
F4_psl2_3=inputFile3[:,2]
F4_psl4_3=inputFile3[:,3]
F4_psl02_3=inputFile3[:,4]
F4_nkval_3=F4_kval_3.shape[0]
print(F4_nkval_3)

F4_kval_4=inputFile4[:,0]
F4_psl0_4=inputFile4[:,1]
F4_psl2_4=inputFile4[:,2]
F4_psl4_4=inputFile4[:,3]
F4_psl02_4=inputFile4[:,4]
F4_nkval_4=F4_kval_4.shape[0]
print(F4_nkval_4)

F4_kval_5=inputFile5[:,0]
F4_psl0_5=inputFile5[:,1]
F4_psl2_5=inputFile5[:,2]
F4_psl4_5=inputFile5[:,3]
F4_psl02_5=inputFile5[:,4]
F4_nkval_5=F4_kval_5.shape[0]
print(F4_nkval_5)

F4_psl0_mean = (F4_psl0_1 + F4_psl0_2 + F4_psl0_3 + F4_psl0_4 + F4_psl0_5)/5
F4_psl2_mean = (F4_psl2_1 + F4_psl2_2 + F4_psl2_3 + F4_psl2_4 + F4_psl2_5)/5
F4_psl4_mean = (F4_psl4_1 + F4_psl4_2 + F4_psl4_3 + F4_psl4_4 + F4_psl4_5)/5
F4_psl02_mean = (F4_psl02_1 + F4_psl02_2 + F4_psl02_3 + F4_psl02_4 + F4_psl02_5)/5

print(F4_psl0_mean.shape[0],F4_psl2_mean.shape[0],F4_psl4_mean.shape[0],F4_psl02_mean.shape[0])

F4_psl0_sqrsum = (F4_psl0_1 - F4_psl0_mean)**2 + (F4_psl0_2 - F4_psl0_mean)**2 + (F4_psl0_3 - F4_psl0_mean)**2 + (F4_psl0_4 - F4_psl0_mean)**2 + (F4_psl0_5 - F4_psl0_mean)**2
F4_psl2_sqrsum = (F4_psl2_1 - F4_psl2_mean)**2 + (F4_psl2_2 - F4_psl2_mean)**2 + (F4_psl2_3 - F4_psl2_mean)**2 + (F4_psl2_4 - F4_psl2_mean)**2 + (F4_psl2_5 - F4_psl2_mean)**2
F4_psl4_sqrsum = (F4_psl4_1 - F4_psl4_mean)**2 + (F4_psl4_4 - F4_psl4_mean)**2 + (F4_psl4_3 - F4_psl4_mean)**2 + (F4_psl4_4 - F4_psl4_mean)**2 + (F4_psl4_5 - F4_psl4_mean)**2
F4_psl02_sqrsum = (F4_psl02_1 - F4_psl02_mean)**2 + (F4_psl02_2 - F4_psl02_mean)**2 + (F4_psl02_3 - F4_psl02_mean)**2 + (F4_psl02_4 - F4_psl02_mean)**2 + (F4_psl02_5 - F4_psl02_mean)**2

F4_psl0_stddev = F4_psl0_sqrsum
F4_psl2_stddev = F4_psl2_sqrsum
F4_psl4_stddev = F4_psl4_sqrsum
F4_psl02_stddev = F4_psl0_sqrsum

for i in range(0,F4_nkval_1):
    F4_psl0_stddev[i] = math.sqrt(F4_psl0_sqrsum[i]/(5-1))
    F4_psl2_stddev[i] = math.sqrt(F4_psl2_sqrsum[i]/(5-1))
    F4_psl4_stddev[i] = math.sqrt(F4_psl4_sqrsum[i]/(5-1))
    F4_psl02_stddev[i] = math.sqrt(F4_psl02_sqrsum[i]/(5-1))


#========= Reading GR files =========

data_path_dir = 'GR/'

sim_name = "z1p000"
ext = ".dat"
model = "GR"

# Reading files
box1_name = data_path_dir+"box1001/pofk_rsd_mu_2d_l_"+model+"_Box1001_"+sim_name+ext
box2_name = data_path_dir+"box2001/pofk_rsd_mu_2d_l_"+model+"_Box2001_"+sim_name+ext
box3_name = data_path_dir+"box3001/pofk_rsd_mu_2d_l_"+model+"_Box3001_"+sim_name+ext
box4_name = data_path_dir+"box4001/pofk_rsd_mu_2d_l_"+model+"_Box4001_"+sim_name+ext
box5_name = data_path_dir+"box5001/pofk_rsd_mu_2d_l_"+model+"_Box5001_"+sim_name+ext

inputFile1 = np.loadtxt(box1_name)
inputFile2 = np.loadtxt(box2_name)
inputFile3 = np.loadtxt(box3_name)
inputFile4 = np.loadtxt(box4_name)
inputFile5 = np.loadtxt(box5_name)

GR_kval_1=inputFile1[:,0]
GR_psl0_1=inputFile1[:,1]
GR_psl2_1=inputFile1[:,2]
GR_psl4_1=inputFile1[:,3]
GR_psl02_1=inputFile1[:,4]
GR_nkval_1=GR_kval_1.shape[0]
print(GR_nkval_1)

GR_kval_2=inputFile2[:,0]
GR_psl0_2=inputFile2[:,1]
GR_psl2_2=inputFile2[:,2]
GR_psl4_2=inputFile2[:,3]
GR_psl02_2=inputFile2[:,4]
GR_nkval_2=GR_kval_2.shape[0]
print(GR_nkval_2)

GR_kval_3=inputFile3[:,0]
GR_psl0_3=inputFile3[:,1]
GR_psl2_3=inputFile3[:,2]
GR_psl4_3=inputFile3[:,3]
GR_psl02_3=inputFile3[:,4]
GR_nkval_3=GR_kval_3.shape[0]
print(GR_nkval_3)

GR_kval_4=inputFile4[:,0]
GR_psl0_4=inputFile4[:,1]
GR_psl2_4=inputFile4[:,2]
GR_psl4_4=inputFile4[:,3]
GR_psl02_4=inputFile4[:,4]
GR_nkval_4=GR_kval_4.shape[0]
print(GR_nkval_4)

GR_kval_5=inputFile5[:,0]
GR_psl0_5=inputFile5[:,1]
GR_psl2_5=inputFile5[:,2]
GR_psl4_5=inputFile5[:,3]
GR_psl02_5=inputFile5[:,4]
GR_nkval_5=GR_kval_5.shape[0]
print(GR_nkval_5)

GR_psl0_mean = (GR_psl0_1 + GR_psl0_2 + GR_psl0_3 + GR_psl0_4 + GR_psl0_5)/5
GR_psl2_mean = (GR_psl2_1 + GR_psl2_2 + GR_psl2_3 + GR_psl2_4 + GR_psl2_5)/5
GR_psl4_mean = (GR_psl4_1 + GR_psl4_2 + GR_psl4_3 + GR_psl4_4 + GR_psl4_5)/5
GR_psl02_mean = (GR_psl02_1 + GR_psl02_2 + GR_psl02_3 + GR_psl02_4 + GR_psl02_5)/5

print(GR_psl0_mean.shape[0],GR_psl2_mean.shape[0],GR_psl4_mean.shape[0],GR_psl02_mean.shape[0])

GR_psl0_sqrsum = (GR_psl0_1 - GR_psl0_mean)**2 + (GR_psl0_2 - GR_psl0_mean)**2 + (GR_psl0_3 - GR_psl0_mean)**2 + (GR_psl0_4 - GR_psl0_mean)**2 + (GR_psl0_5 - GR_psl0_mean)**2
GR_psl2_sqrsum = (GR_psl2_1 - GR_psl2_mean)**2 + (GR_psl2_2 - GR_psl2_mean)**2 + (GR_psl2_3 - GR_psl2_mean)**2 + (GR_psl2_4 - GR_psl2_mean)**2 + (GR_psl2_5 - GR_psl2_mean)**2
GR_psl4_sqrsum = (GR_psl4_1 - GR_psl4_mean)**2 + (GR_psl4_4 - GR_psl4_mean)**2 + (GR_psl4_3 - GR_psl4_mean)**2 + (GR_psl4_4 - GR_psl4_mean)**2 + (GR_psl4_5 - GR_psl4_mean)**2
GR_psl02_sqrsum = (GR_psl02_1 - GR_psl02_mean)**2 + (GR_psl02_2 - GR_psl02_mean)**2 + (GR_psl02_3 - GR_psl02_mean)**2 + (GR_psl02_4 - GR_psl02_mean)**2 + (GR_psl02_5 - GR_psl02_mean)**2

GR_psl0_stddev = GR_psl0_sqrsum
GR_psl2_stddev = GR_psl2_sqrsum
GR_psl4_stddev = GR_psl4_sqrsum
GR_psl02_stddev = GR_psl0_sqrsum

for i in range(0,GR_nkval_1):
    GR_psl0_stddev[i] = math.sqrt(GR_psl0_sqrsum[i]/(5-1))
    GR_psl2_stddev[i] = math.sqrt(GR_psl2_sqrsum[i]/(5-1))
    GR_psl4_stddev[i] = math.sqrt(GR_psl4_sqrsum[i]/(5-1))
    GR_psl02_stddev[i] = math.sqrt(GR_psl02_sqrsum[i]/(5-1))


#========= Computing ratios =========

# Monopole
GR_psl0_up = GR_psl0_mean+GR_psl0_stddev
GR_psl0_down = GR_psl0_mean-GR_psl0_stddev
GR_psl0_up_rel = GR_psl0_up / GR_psl0_mean
GR_psl0_down_rel = GR_psl0_down / GR_psl0_mean

F4_psl0_up = F4_psl0_mean+F4_psl0_stddev
F4_psl0_down = F4_psl0_mean-F4_psl0_stddev
F4_psl0_up_rel = F4_psl0_up / GR_psl0_mean
F4_psl0_down_rel = F4_psl0_down / GR_psl0_mean

# Quadrupole
GR_psl2_up = GR_psl2_mean+GR_psl2_stddev
GR_psl2_down = GR_psl2_mean-GR_psl2_stddev
GR_psl2_up_rel = GR_psl2_up / GR_psl0_mean
GR_psl2_down_rel = GR_psl2_down / GR_psl0_mean

F4_psl2_up = F4_psl2_mean+F4_psl2_stddev
F4_psl2_down = F4_psl2_mean-F4_psl2_stddev
F4_psl2_up_rel = F4_psl2_up / GR_psl0_mean
F4_psl2_down_rel = F4_psl2_down / GR_psl0_mean

# Octupole
GR_psl4_up = GR_psl4_mean+GR_psl4_stddev
GR_psl4_down = GR_psl4_mean-GR_psl4_stddev
GR_psl4_up_rel = GR_psl4_up / GR_psl0_mean
GR_psl4_down_rel = GR_psl4_down / GR_psl0_mean

F4_psl4_up = F4_psl4_mean+F4_psl4_stddev
F4_psl4_down = F4_psl4_mean-F4_psl4_stddev
F4_psl4_up_rel = F4_psl4_up / GR_psl0_mean
F4_psl4_down_rel = F4_psl4_down / GR_psl0_mean

#========= Plotting =========

# Monopole

fig1 = plt.figure()

# print the shot noise subtracted P(k) :: activate!

plt.plot(GR_kval_5, GR_psl0_up_rel,color='lightblue')
plt.plot(GR_kval_5, GR_psl0_down_rel,color='lightblue')
plt.plot(F4_kval_5, F4_psl0_up_rel,color='orange')
plt.plot(F4_kval_5, F4_psl0_down_rel,color='orange')

# format the axes
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_0(k)/P_0^{GR}$")

#plt.ylabel(r"$P_l(k)$ [$h^{-3}\mathrm{Mpc}^3$]")

#plt.xlim(0.005, 0.6)
#plt.xlim(0.01, 10.0)
#plt.ylim(1,100000)
plt.xscale('log')
#plt.yscale('log')

pofk_rsd_pdf_name = 'pofk_rsd_'+sim_combined_name+'_l0'+'.pdf'

fig1.savefig(pofk_rsd_pdf_name)


# Quadrupole

fig2 = plt.figure()

# print the shot noise subtracted P(k) :: activate!

plt.plot(GR_kval_5, GR_psl2_up_rel,color='lightblue')
plt.plot(GR_kval_5, GR_psl2_down_rel,color='lightblue')
plt.plot(F4_kval_5, F4_psl2_up_rel,color='orange')
plt.plot(F4_kval_5, F4_psl2_down_rel,color='orange')

# format the axes
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_2(k)/P_0^{GR}$")

#plt.xlim(0.005, 0.6)
#plt.xlim(0.01, 10.0)
#plt.ylim(1,100000)
plt.xscale('log')
#plt.yscale('log')

pofk_rsd_pdf_name = 'pofk_rsd_'+sim_combined_name+'_l2'+'.pdf'

fig2.savefig(pofk_rsd_pdf_name)


# Octupole

fig4 = plt.figure()

# print the shot noise subtracted P(k) :: activate!

plt.plot(GR_kval_5, GR_psl4_up_rel,color='lightblue')
plt.plot(GR_kval_5, GR_psl4_down_rel,color='lightblue')
plt.plot(F4_kval_5, F4_psl4_up_rel,color='orange')
plt.plot(F4_kval_5, F4_psl4_down_rel,color='orange')

# format the axes
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P_4(k)/P_0^{GR}$")

#plt.xlim(0.005, 0.6)
#plt.xlim(0.01, 10.0)
#plt.ylim(1,100000)
plt.xscale('log')
#plt.yscale('log')

pofk_rsd_pdf_name = 'pofk_rsd_'+sim_combined_name+'_l4'+'.pdf'

fig4.savefig(pofk_rsd_pdf_name)

