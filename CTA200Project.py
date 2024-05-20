#!/usr/bin/env python
# coding: utf-8

# In[3]:


#importing the necessary packages
from astropy.io import fits
import matplotlib.pyplot as plt


# In[ ]:


#!ls /raid/DESI/spectra
#this is me looking in the desi raid folder to see where my data is, I want the data in the
# 'sp_x_apogee_x_spspectra_rvtab' folder


# In[ ]:


#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits' #this is me 
#finding the data via its path
#with fits.open(fits_file_path) as hdul: #opening it in the hdul 
    #hdul.info() 
    #for idx, hdu in enumerate(hdul):
        #print(f"HDU #{idx + 1}:", hdu.header) #with this code, I'm just printing out the 
        #headers, which I'll use to identify the data I need. 


# I thought I might want to print out all the data, but decided first to check just how much data was in one folder...

# In[ ]:


#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data #the [1] refers to which hdul extension the data I want is in
    #fe_h_sp_column = data['FeH']

#len(fe_h_sp_column), this printed 7910, Safe to say, I will no longer be printing out all the data


# My first goal is to plot the metallicty of these 7910 stars using the APOGEE and then DESI data. I'm starting with the plot that shows the metallicty using APOGEE, or the Fe/H count.

# In[19]:


#now i need to create the fake data
import random

# Create the dictionary
data = {
    'feh_apogee': random.randint(1, 1000),
    'feh_desi': random.randint(1, 1000),
    'desi_logg': random.randint(1, 1000),
    'desi_teff': random.randint(1, 1000),
    'apogee_logg': random.randint(1, 1000),
    'apogee_teff': random.randint(1, 1000),
    'ra': random.randint(0, 360), #ra has to be between 0 and 360
    'dec': random.randint(-90, 90), #dec has to be between 90 and -90
    'blue': random.randint(1, 1000),
    'red': random.randint(1, 1000),
    'nir': random.randint(1, 1000),
    'blue_flux': random.randint(1, 1000),
    'red_flux': random.randint(1, 1000),
    'nir_flux': random.randint(1, 1000)
}

# Print the dictionary
print(data)


# This will be completely innacurate and just for fun to show my plots do indeed plot! 

# In[4]:


#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data #the [1] refers to which hdul extension the data I want is in
    #desi_fe_h_column = data['FeH'] #here, I want the data under 'FeH'

desi_fe_h_column = data['feh_desi']
plt.hist(desi_fe_h_column, color="skyblue")
#plt.hist(desi_fe_h_column, bins=120, color="skyblue", range=(-1.5, 0.5)), this is the range
#observed when working with real data and the number of bins used
plt.xlabel('DESI [Fe/H]')
plt.ylabel('count')
plt.title('Histogram of Fe abundance')
plt.show()


# In[5]:


#to get a better look at the data, I'm going to make a transparent graph where only the top
#is outlined, for this I need step lines at the top of the bins using plt.step
#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data 
    #desi_fe_h_column = data['FeH'] 
    #print(desi_fe_h_column)


#n, bins, _ = plt.hist(desi_fe_h_column, bins=120, range=(-1.5, 0.5), color='white')
n, bins, _ = plt.hist(desi_fe_h_column, color='white')

plt.step(bins[:-1], n, where='post', color='black', linewidth=1)

plt.xlabel('DESI [Fe/H]')
plt.ylabel('Count')
plt.title('Histogram of Fe abundance')
plt.show()


# In[6]:


#I'm doing the same two graohs again for APOGEE's data
#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[2].data
    #apogee_fe_h_column = data['Fe_H']
    #print(apogee_fe_h_column)

apogee_fe_h_column = data['feh_apogee']
#plt.hist(apogee_fe_h_column, bins=120, color="skyblue", range=(-1.5, 0.5))
plt.hist(apogee_fe_h_column, color="skyblue")
plt.xlabel('APOGEE [Fe/H]')
plt.ylabel('count')
plt.title('Histogram of Fe abundance')
plt.show()


# In[7]:


#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[2].data
    #apogee_fe_h_column = data['Fe_H']
    #print(apogee_fe_h_column)
    
#n, bins, _ = plt.hist(apogee_fe_h_column, bins=120, range=(-1.5, 0.5), color='white')
n, bins, _ = plt.hist(apogee_fe_h_column, color='white')

# Plot a step line at the top of the histogram bins
plt.step(bins[:-1], n, where='post', color='black', linewidth=1)

plt.xlabel('APOGEE [Fe/H]')
plt.ylabel('count')
plt.title('Histogram of Fe abundance')
plt.show()


# In[12]:


#here, my goal is to plot the Teff data on the x axis and Log g on the y axis, starting with apogee

#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[2].data
    #LOGG_column = data['LOGG']
    #TEFF_column = data['TEFF']
    #print(LOGG_column)
    #print(TEFF_column)
  
TEFF_column = data['apogee_teff']
LOGG_column = data['apogee_logg']
#plt.figure(figsize=(6,8.5))
plt.scatter(TEFF_column, LOGG_column, s=0.5, c="red")
#plt.xlim(8000, 2000)
#plt.ylim(6, 0)
plt.xlabel('TEFF in K')
plt.ylabel('LOGG in cgs')
plt.title('APOGEE - Stellar Parameter: Teff/log(g)')
plt.show()


# In[16]:


#Same thing, but with DESI
#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data
    #LOGG_column = data['LOGG']
    #TEFF_column = data['TEFF']
    #print(LOGG_column)
    #print(TEFF_column)

TEFF_column = data['desi_teff']
LOGG_column = data['desi_logg']
#plt.figure(figsize=(6,8.5))
plt.scatter(TEFF_column, LOGG_column, s=0.5, c="red")
#plt.xlim(8000, 2000)
#plt.ylim(6, 0)
plt.xlabel('TEFF in K')
plt.ylabel('LOGG in cgs')
plt.title('DESI - Stellar Parameter: Teff/log(g)')
plt.show()


# In[20]:


#Here, to understand how APOGEE collectds its data and better understand its success, 
#I'm looking at the location of all the data collected
#from astropy.io import fits
#import numpy as np
#import matplotlib.pyplot as plt
#from astropy.coordinates import SkyCoord
#import astropy.units as u

#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data
    #RA = data['TARGET_RA']
    #DEC = data['TARGET_DEC']
    #print(RA)
    #print(DEC)
    
# Convert the coordinates properly
#coords = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')

RA = data['ra']
DEC = data['dec']
plt.figure(figsize=(10, 5))
plt.subplot(111, projection='mollweide')
plt.title('Mollweide Projection of RA and DEC')
plt.grid(True)
plt.scatter(RA, DEC, s=1)
#plt.scatter(coords.ra.wrap_at(180*u.deg).radian, coords.dec.radian, s=1), it was necessary before to further convert
plt.show()


# In[23]:


#Now I'm going to look at a spectra of one specific star, I have my choice from 0 to 7909
#There are three wavelengths the data looked at, blue, red and near-infrared
#import matplotlib.pyplot as plt
#import numpy as np

#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[3].data
    
    #blue_wavelength = data['B_WAVELENGTH'][0]
    #blue_flux = data['flx_B'][0]
    
    #red_wavelength = data['R_WAVELENGTH'][0]
    #red_flux = data['flx_R'][0]
    
    #nir_wavelength = data['Z_WAVELENGTH'][0]
    #nir_flux = data['flx_Z'][0]

blue_wavelength = data['blue']
red_wavelength = data['red']
nir_wavelength = data['nir']

blue_flux = data['blue']
red_flux = data['red']
nir_flux = data['nir']


#plt.figure(figsize=(15, 8))
plt.plot(blue_wavelength, blue_flux, label='Blue arm', color='blue', linewidth=1)
plt.xlabel('Wavelength in Angstrom')
plt.ylabel('Flux in cm^-2*s^-1*A^-1')
plt.title('Spectra')
plt.legend()
plt.show()

#plt.figure(figsize=(15, 8))
plt.plot(red_wavelength, red_flux, label='Red arm', color='red', linewidth=1)
plt.xlabel('Wavelength in Angstrom')
plt.ylabel('Flux in cm^-2*s^-1*A^-1')
plt.title('Spectra')
plt.legend()
plt.show()

#plt.figure(figsize=(15, 8))
plt.plot(nir_wavelength, nir_flux, label='NIR arm', color='orange', linewidth=1)
plt.xlabel('Wavelength in Angstrom')
plt.ylabel('Flux in cm^-2*s^-1*A^-1')
plt.title('Spectra')
plt.legend()
plt.show()

#plt.figure(figsize=(15, 8))
plt.plot(blue_wavelength, blue_flux, label='Blue arm', color='blue', linewidth=1)
plt.plot(red_wavelength, red_flux, label='Red arm', color='red', linewidth=1)
plt.plot(nir_wavelength, nir_flux, label='NIR arm', color='orange', linewidth=1)
plt.xlabel('Wavelength in Angstrom')
plt.ylabel('Flux in cm^-2*s^-1*A^-1')
plt.title('Spectra')
plt.legend()
plt.show()


# I was curious which star that is, so I'm looking for the its coordinates

# In[ ]:


#fits_file_path = '/raid/DESI/spectra/sp_x_apogee_x_spspectra_rvtab.fits'
#with fits.open(fits_file_path) as hdul:
    #hdul.info()
    #data = hdul[1].data
    #RA = data['TARGET_RA'][0]
    #DEC = data['TARGET_DEC'][0]
    #print(RA)
    #print(DEC)
#googling them gave me nothing :(


# Thank you for joining along! Hope you liked what you saw and I had an amazing time in this course. Best of luck grading everything!!
