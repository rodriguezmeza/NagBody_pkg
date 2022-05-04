
# source activate nbodykit-env

#========= Importing libs =========

from nbodykit.lab import *
#from nbodykit import setup_logging, style

import os
import matplotlib.pyplot as plt
#plt.style.use(style.notebook)

import numpy as np
import math

#========= Setting parameters =========

# Tune kmin and kmax according to particle positions ranges (length unit of the simulation box)
# Default kmax = kNyq = pi Nmesh/Lbox, and 2pi/Lbox its the minimum value of dk.
myNmesh=512
#mykmin=0.000001
#mykmax=0.003
mykmin=0.001
mykmax=8.0
#nkval=500
#mydk=(mykmax-mykmin)/(nkval - 1)
Lbox = 200.0
mydk = 2 * 3.14159265 / Lbox
nkval = 1.0 + (mykmax-mykmin)/mydk
xr1=0.01
xr2=10
yr1=1
yr2=100000
mylos=[0,0,1]
print('dk = ', mydk)
print('nkval = ', nkval)

# To have units in Mpc : otherwise let them have the value 1.0
#kunit=1e3
#psunit=1e-9
kunit=1.0
psunit=1.0
# Change below also (mpc_unit)

sim_box_name = "F4_Box1001_z0p000"

# define download dir for the data to the current directory
download_dir = "./output"

# NOTE: change this path if you downloaded the data somewhere else!
data_path = os.path.join(download_dir, 'snapshot_z0p000.0')

# initialize the Gadget1 catalog objects for data and randoms
data = Gadget1Catalog(data_path)

print('data columns = ', data.columns)
print(' ')
print('data = ', data.attrs)
print(' ')
Redshift = data.attrs['Redshift']
az = 1./(Redshift+1)
print('Redshift =', Redshift)
print('a =', az)
aHofa = az*math.sqrt(data.attrs['Omega0']/(az*az*az) + data.attrs['OmegaLambda'])
print('Hubble (z) =', aHofa/az)
print('aHofa =', aHofa)
#mpc_unit = 1000.0 # use 1.0 if units are Mpc.
mpc_unit = 1.0 # use 1.0 if units are Mpc.

pos = data['Position']
nsnap = pos.shape[0]
print(nsnap)
print(pos)
vel=data['GadgetVelocity']
print(vel)

# convert to a MeshSource, using TSC interpolation on 256^3 mesh:: tsc or cic
mesh = data.to_mesh(resampler='cic', Nmesh=myNmesh, compensated=True, position='Position')
#mesh = data.to_mesh(resampler='tsc', Nmesh=myNmesh, compensated=True, position='Position',interlaced=True)
print("mesh = ", mesh)

#=============== Begin: PS =================
# compute the power, specifying desired linear k-binning
r = FFTPower(mesh, mode='1d', dk=mydk, kmin=mykmin, kmax=mykmax)

# the result is stored at "power" attribute
Pk = r.power
print(Pk)

print(Pk.coords)

# Saving to file PS vs k
kvec=Pk['k'] * kunit
pkvec=Pk['power'].real * psunit
pkvec2=(Pk['power'].real - Pk.attrs['shotnoise']) * psunit

nk = kvec.shape[0]

print(nk)

pofka=np.array([kvec,pkvec,pkvec2])
pofkb=np.transpose(pofka)
print pofkb

pofk_txt_name = 'pofk_'+sim_box_name+'.dat'

np.savetxt(pofk_txt_name,pofkb,delimiter='\t',newline='\r\n')

# print out the meta-data
for k in Pk.attrs:
    print("%s = %s" %(k, str(Pk.attrs[k])))
    
    
#--------------- Begin: Fig1 ----------------

fig1 = plt.figure()

# print the shot noise subtracted P(k)
plt.loglog(Pk['k']*kunit, (Pk['power'].real - Pk.attrs['shotnoise'])*psunit)

# format the axes
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$P(k)$ [$h^{-3}\mathrm{Mpc}^3$]")
#plt.xlim(0.001, 0.6)
#plt.xlim(0.01, 10.0)
#plt.ylim(1,100000)

pofk_pdf_name = 'pofk_'+sim_box_name+'.pdf'

fig1.savefig(pofk_pdf_name)


#--------------- End: Fig1 ----------------

#=============== End: PS =================

