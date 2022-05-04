
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

#=============== Begin: RSD =================

# add RSD
# Using the formula: s = x + v_p.n_z / (a H) = x + sqrt(a) V_gadget*n_z /(a H)
#                      = x + V_gadget*n_z /(sqrt(a) H), H = 100 h, in km/s/kpc.
vnorm = mpc_unit /(100.0*aHofa)/math.sqrt(1.0 + Redshift)
print('vnorm =', vnorm)
data['RSDPosition'] = data['Position'] + vnorm * data['GadgetVelocity'] * mylos

# convert to a MeshSource, using TSC (or CIC) interpolation on 256^3 mesh
mesh_rsd = data.to_mesh(resampler='cic', Nmesh=myNmesh, compensated=True, position='RSDPosition')

# compute the 2D power AND ell=0,2,4 multipoles
r = FFTPower(mesh_rsd, mode='2d', dk=mydk, kmin=mykmin, kmax=mykmax, Nmu=5, los=mylos, poles=[0,2,4])

poles = r.poles
print(poles)
print("variables = ", poles.variables)

print(poles)

# Save to a file
kvec=poles['k'] * kunit
pk0vec=poles['power_0'].real * psunit
pk2vec=poles['power_2'].real * psunit
pk4vec=poles['power_4'].real * psunit
pk0vec2=(poles['power_0'].real - poles.attrs['shotnoise']) * psunit

nk = kvec.shape[0]

print(nk)
print(pk0vec)

#pofka=np.array([kvec,pk0vec,pk2vec,pk4vec,kvec*pk0vec,kvec*pk2vec,kvec*pk4vec])
pofka=np.array([kvec,pk0vec,pk2vec,pk4vec,pk0vec2])
pofkb=np.transpose(pofka)
print pofkb

pofk_rsd_txt_name = 'pofk_rsd_mu_2d_l_'+sim_box_name+'.dat'

np.savetxt(pofk_rsd_txt_name,pofkb,delimiter='\t',newline='\r\n')

#--------------- Begin: Fig5 ----------------
fig5 = plt.figure()

for ell in [0, 2, 4]:
    label = r'$\ell=%d$' % (ell)
    P = poles['power_%d' %ell].real
    if ell == 0: P = P - poles.attrs['shotnoise']
    plt.plot(poles['k']*kunit, poles['k']*kunit * P*psunit, label=label)


# format the axes
plt.legend(loc=0)
plt.xlabel(r"$k$ [$h \ \mathrm{Mpc}^{-1}$]")
plt.ylabel(r"$k \ P_\ell$ [$h^{-2} \mathrm{Mpc}^2$]")
plt.xlim(0.0, 1)
#plt.ylim(-1000,3000)
#plt.ylim(-2000,3000)

pofk_rsd_pdf_name = 'pofk_rsd_mu_2d_l_'+sim_box_name+'.pdf'

fig5.savefig(pofk_rsd_pdf_name)

#--------------- End: Fig5 ----------------

#=============== End: RSD =================
