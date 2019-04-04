#!/usr/bin/env python2.7

'''A tool for determining the UVIT FUV and NUV filters.


   Copyright 2018 Prajwel Joseph
  
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
  
       http://www.apache.org/licenses/LICENSE-2.0
  
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License. 

   Changes; when, what
   -------------------
   Jan 10, 2018: Galactic latitude check incorporated.


   The author would like to acknowledge inputs from Dr. Koshy George
   which greatly helped the developement of this script.
   
'''


import os
import sys
import string
import urllib
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.io import ascii
from requests import Session
from bs4 import BeautifulSoup
from astropy import units as u
from matplotlib.colors import LogNorm
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord


# To get the user input. 
instrument = str(sys.argv[1])
RA = str(sys.argv[2])
DEC = str(sys.argv[3])
#working_arena = str(sys.argv[4])

#instrument = 'uvit'
#RA = "7:36:51.396"
#DEC = "65:36:9.170"
working_arena = '.'

# To do all the stuff in a specific directory.
os.chdir(working_arena)

# To read the TD1 catalogue.
try:
    catalogue = 'td1_catalogue.fits'
    hdu = fits.open(catalogue)
except IOError:
    print('Could not find the catalogue file: {}'.format(catalogue))
    sys.exit(1)

# instrument and radius of search in arsec.
field_radius = {'uvit'  : 1200,
                'sxt'   : 1500,
                'czti'  : 1680,
                'laxpc' : 1680}

# Functions to convert magnitude to UVIT count rates.
flux_norm = 2E-13
def countnuv(flux):
    silica = 955.0
    b4 = 218.5
    b13 = 275.8
    b15 = 59.6
    n2 = 50.6
    flux_ratio = flux / flux_norm
    cr1 =  silica * flux_ratio
    cr2 =  b4 * flux_ratio
    cr3 =  b13 * flux_ratio
    cr4 =  b15 * flux_ratio
    cr5 =  n2 * flux_ratio
    return flux, cr1, cr2, cr3, cr4, cr5

def countfuv(flux):
    caf2 = 74.5
    baf2 = 60.0
    sapphire = 50.0
    silica = 17.3
    flux_ratio = flux / flux_norm
    cr1 =  caf2 * flux_ratio    
    cr2 =  baf2 * flux_ratio
    cr3 =  sapphire * flux_ratio
    cr4 =  silica * flux_ratio 
    return flux, cr1, cr2, cr3, cr4

# Function to find seperation in celestial coordinates.
cc = SkyCoord(RA, DEC, unit = (u.hourangle, u.deg))
def cel_separation(a, b):
    coo = SkyCoord(a, b, frame = 'icrs', unit = 'deg')
    return coo.separation(cc)

# To check if Galactic latitude is between -30 to 30.
gal_lat = cc.galactic.b.value
gal_plane = 'no'
if -30.0 <= gal_lat <= 30.0:
    gal_plane_warning = 'The galactic latitude is between -30 to 30. \
                        \nYour field cannot be checked using TD1 catalogue!'
    print('\n{}\n'.format(gal_plane_warning))
    with open('gal_plane_warning.txt', 'w') as gal_warn:
        gal_warn.write(gal_plane_warning)
    sys.exit(1)

# Reading coordinates and fluxes from catalogue.
hdu = fits.open(catalogue)
alpha = hdu[1].data['ra']
delta = hdu[1].data['dec']
nuv_flux = hdu[1].data['flux_2365_a']
fuv_flux = hdu[1].data['flux_1565_a']

# NUV 
refined_set = [(al, de, nf) for al, de, nf 
                            in zip(alpha, delta, nuv_flux)
                            if (cc.ra.value - 5) <= al <= (cc.ra.value + 5)
                                and (cc.dec.value - 5) <= de <= (cc.dec.value + 5)]

nalpha, ndelta, nuv_flux = zip(*refined_set)

confined_set = [nf for al, de, nf 
                   in zip(nalpha, ndelta, nuv_flux) 
                   if cel_separation(al, de) <= field_radius[instrument] * u.arcsec]

# If list is empty, normal value need to be taken.
if len(confined_set) == 0:
    confined_set.append(flux_norm)

nd = sorted(confined_set)[-1]
flux, ta, tb, tc, td, te = countnuv(nd)
nuv_res = Table([[flux], [ta], [tb], [tc], [td], [te]],
               names = ('flux_2365_a',
                        'silica',
                        'b4',
                        'b13',
                        'b15',
                        'n2'), 
               meta = {'name': 'NUV counts'})

#nuv_res['flux_2365_a'].format = '
nuv_res['silica'].format = '4.1f'
nuv_res['b4'].format = '4.1f'
nuv_res['b13'].format = '4.1f'
nuv_res['b15'].format = '4.1f'
nuv_res['n2'].format = '4.1f'

print('\n\n### NUV\n\n{}\n'.format(nuv_res))

# To select NUV safe filters.
nuv_filter_dict = {0: 'Silica', 1: 'NUV-B4', 2: 'NUV-B13', 3: 'NUV-B15', 4: 'NUV-N2'}
i = 0
nuv_safe = []
for Filter in zip(*nuv_res['silica','b4','b13','b15','n2']):
     if sum(np.array(Filter) > 1500) == 0: 
         nuv_safe.append(nuv_filter_dict[i])
     if i == 0:
         if sum(np.array(Filter) > 1133) == 0:
             nuv_safe.append('NUV-grating')
     i = i + 1

nuv_declaration = 'Safe filters in NUV: {}'.format(nuv_safe)
print('\n\n{}\n'.format(nuv_declaration))

# To write to file.
nuv_table = 'NUV_td1-nd-int.txt'
ascii.write(nuv_res, nuv_table, format = 'csv', overwrite = True)

with open('safe_NUV_filters.txt', 'w') as safe_file:
    safe_file.write(nuv_declaration)

# FUV 
refined_set = [(al, de, ff) for al, de, ff 
                            in zip(alpha, delta, fuv_flux)
                            if (cc.ra.value - 5) <= al <= (cc.ra.value + 5)
                                and (cc.dec.value - 5) <= de <= (cc.dec.value + 5)]

nalpha, ndelta, fuv_flux = zip(*refined_set)

confined_set = [ff for al, de, ff 
                   in zip(nalpha, ndelta, fuv_flux) 
                   if cel_separation(al, de) <= field_radius[instrument] * u.arcsec]

# If list is empty, normal value need to be taken.
if len(confined_set) == 0:
    confined_set.append(flux_norm)

fd = sorted(confined_set)[-1]
flux, ta, tb, tc, td = countfuv(fd)
fuv_res = Table([[flux], [ta], [tb], [tc], [td]],
               names = ('flux_1565_a',
                        'caf2',
                        'baf2',
                        'sapphire',
                        'silica'), 
               meta = {'name': 'NUV counts'})

#fuv_res['flux_1565_a'].format = '
fuv_res['caf2'].format = '4.1f'
fuv_res['baf2'].format = '4.1f'
fuv_res['sapphire'].format = '4.1f'
fuv_res['silica'].format = '4.1f'

print('\n### FUV \n\n{}\n\n'.format(fuv_res))

# To select FUV safe filters.
fuv_filter_dict = {0: 'CaF2', 1: 'BaF2', 2: 'Sapphire', 3: 'Silica'}
j = 0
fuv_safe = []
for Filter in zip(*fuv_res['caf2','baf2','sapphire','silica']):
     if sum(np.array(Filter) > 1500) == 0: 
         fuv_safe.append(fuv_filter_dict[j])
     if j == 0:
         if sum(np.array(Filter) > 892) == 0:
             fuv_safe.append('FUV-grating')
     j = j + 1

fuv_declaration = 'Safe filters in FUV: {}'.format(fuv_safe)
print('\n\n{}\n'.format(fuv_declaration))

# To write to file.
fuv_table = 'FUV_td1-fd-int.txt'
ascii.write(fuv_res, fuv_table, format = 'csv', overwrite = True)

with open('safe_FUV_filters.txt', 'w') as safe_file:
    safe_file.write(fuv_declaration)

print('Done!\n')





