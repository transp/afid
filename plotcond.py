#!/usr/bin/env python3

# Activate the transp_testing virtual environment to access the matplotlib and
# netCDF4 packages

import sys
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from math import sqrt

if len(sys.argv) < 2:
    print("""Usage: %s <netCDF condensed particle filename>\n"""%sys.argv[0])
    exit()

ds = nc.Dataset(sys.argv[1], mode="r")
if 'nptcl' in ds.dimensions.keys():
    nparts = len(ds.dimensions['nptcl'])
else:
    print("Error: particle count dimension not found.\n")
    exit(1)

print("%d particles in file."%nparts)

# Check for required variables
req = ('pphi', 'mu', 'E', 'weight', 'sgn_v')
for v in req:
    if not v in ds.variables.keys():
        print("Error: %s variable not found.\n"%v)
        exit(1)

# Read data into memory
KE = np.copy(ds.variables['E'][:,0])
emin = min(KE);  emax = max(KE)
mu = np.copy(ds.variables['mu'][:,0])
mmin = min(mu);  mmax = max(mu)
pphi = np.copy(ds.variables['pphi'][:,0])
pmin = min(pphi);  pmax = max(pphi)
wt = np.copy(ds.variables['weight'][:,0])
ds.close()

nbins = 128
maxbins = int(min([nparts/100.0, nbins]))

print("Plotting histograms...")

# 1D Histograms
plt.subplot(2,3,1)
dxip = maxbins/(pmax - pmin)
p = np.arange(maxbins+1)/dxip + pmin
phist = np.zeros([maxbins+1])
for j in range(nparts):
    ibinp = int(dxip*(pphi[j] - pmin))
    phist[ibinp] += wt[j]
plt.xlabel('Pphi');  plt.ylabel('weighted particle count')
#plt.hist(pphi, bins=maxbins)
plt.plot(p,phist)

plt.subplot(2,3,2)
dxim = maxbins/(mmax - mmin)
m = np.arange(maxbins+1)/dxim + mmin
mhist = np.zeros([maxbins+1])
for j in range(nparts):
    ibinm = int(dxim*(mu[j] - mmin))
    mhist[ibinm] += wt[j]
plt.xlabel('mu');  plt.ylabel('weighted particle count')
plt.plot(m,mhist)

plt.subplot(2,3,3)
dxie = maxbins/(emax - emin)
E = np.arange(maxbins+1)/dxie + emin
ehist = np.zeros([maxbins+1])
for j in range(nparts):
    ibine = int(dxie*(KE[j] - emin))
    ehist[ibine] += wt[j]
plt.xlabel('K.E.');  plt.ylabel('weighted particle count')
plt.plot(E,ehist)

# 2D Histograms
maxbins = int(min([96.0, 0.05*sqrt(nparts)]))
print("2D histograms have %d ^ 2 bins"%maxbins)

dxip = maxbins/(pmax - pmin)
p = np.arange(maxbins+1)/dxip + pmin
dxim = maxbins/(mmax - mmin)
m = np.arange(maxbins+1)/dxim + mmin
dxie = maxbins/(emax - emin)
E = np.arange(maxbins+1)/dxie + emin

thistpm = np.zeros([maxbins+1,maxbins+1])
thistpe = np.zeros([maxbins+1,maxbins+1])
thistem = np.zeros([maxbins+1,maxbins+1])

for j in range(nparts):
    ibinp = int(dxip*(pphi[j] - pmin))
    ibinm = int(dxim*(mu[j] - mmin))
    ibine = int(dxie*(KE[j] - emin))
    thistpm[ibinm][ibinp] += wt[j]
    thistpe[ibine][ibinp] += wt[j]
    thistem[ibinm][ibine] += wt[j]

plt.subplot(2,3,4)
plt.xlabel('Pphi');  plt.ylabel('mu')
plt.contourf(p,m,thistpm,255)

plt.subplot(2,3,5)
plt.xlabel('Pphi');  plt.ylabel('K.E.')
plt.contourf(p,E,thistpe,255)

plt.subplot(2,3,6)
plt.xlabel('K.E.');  plt.ylabel('mu')
plt.contourf(E,m,thistem,255)

plt.show()
