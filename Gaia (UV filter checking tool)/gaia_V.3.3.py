#!/usr/bin/env python2.7

'''A tool for determining the UVIT FUV and NUV filters.


   Copyright 2017 Prajwel Joseph
  
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
   Dec 22, 2017: bug fixes.
   Dec 23, 2017: deals with cases where FUV is not present.
   Dec 26, 2017: the images are plotted with marked positions.
   Dec 27, 2017: galex search radius changed to .001 degree.
                 Safe filters are explicitely declared. 
   Jan 03, 2018: Galex calibration equation bounded. 
   Jan 04, 2018: Bug Fixes. 
   Jan 05, 2018: Absence of Galex tiles are declared explicitly. 
                 Galactic plane warning included.
   Jan 16, 2018: Error outputs directed to file.
   Jan 22, 2018: GALEX catalogue, GALEX image, and TD1 catalogue
                 methods merged into one script.
   Feb 12, 2018: A bug related to masking of the image fixed.  
   Mar 14, 2018: Bug fixes.  
   Mar 18, 2018: The circles were too big, this has been corrected.  


   The author would like to acknowledge inputs from Dr. Koshy George
   and Dr. C. S. Stalin which greatly helped the developement
   of this script.
   
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

# An error file.
error_file = open('error.txt', 'w')

# To read the TD1 catalogue.
try:
    td1_catalogue = 'td1_catalogue.fits'
    td1_hdu = fits.open(td1_catalogue)
except IOError:
    where_is_td1 = 'Could not find the catalogue file: {}'.format(td1_catalogue)
    print(where_is_td1)
    error_file.write(where_is_td1)
    sys.exit(1)

# To check the RA, DEC user input.
if (DEC.count(':') == RA.count(':') == 2):
    pass
else:
    check_input = 'Check your RA DEC input.'
    print('\n{}\n'.format(check_input))
    error_file.write(check_input)
    sys.exit(1)

# instrument and radius of search in arsec.
field_radius = {'uvit'  : 1200,
                'sxt'   : 1500,
                'czti'  : 1680,
                'laxpc' : 1680}

# instrument and required window size.
field_radius_im = {'uvit'  : 800,
                   'sxt'   : 1500,
                   'czti'  : 1120,
                   'laxpc' : 1120}

# Functions to convert GALEX CPS to AB magnitude.
NUV_ZPMAG = 20.08
FUV_ZPMAG = 18.82
def magfuv(fl):
    return (-2.5 * np.log10(np.array(fl)) + FUV_ZPMAG)

def magnuv(fl):
    return (-2.5 * np.log10(np.array(fl)) + NUV_ZPMAG)

# Functions to convert GALEX magnitude to UVIT count rates.
def countfuv(mg):
    caf2 = 1.0
    baf2 = 0.85
    sapphire = 0.63
    silica = 0.22
    mg1 = 18.22

    if 10.511 <= mg <= 15.0:
        mg_c = 5.371 + (20.0 * mg - 210.2) ** 0.5
    elif mg < 10.511:
        mg_c = 5.371
    else: 
        mg_c = mg

    cr1 =  caf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr2 =  baf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr3 =  sapphire *10.0 ** ((mg1 - mg_c) * 0.4)    
    cr4 =  silica * 10.0 ** ((mg1 - mg_c) * 0.4)    
    return mg, mg_c, cr1, cr2, cr3, cr4

countfuv = np.vectorize(countfuv)

def countfuv_abs(mg): # for cases where FUV is absent.
    caf2 = 1.0
    baf2 = 0.85
    sapphire = 0.63
    silica = 0.22
    mg1 = 18.22

    if 9.323 <= mg <= 15.0:
        mg_c = 2.634 + (26.316 * mg - 245.329) ** 0.5
    elif mg < 9.323:
        mg_c = 2.634
    else: 
        mg_c = mg

    mg_c = mg_c - 1.65   
    cr1 =  caf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr2 =  baf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr3 =  sapphire *10.0 ** ((mg1 - mg_c) * 0.4)    
    cr4 =  silica * 10.0 ** ((mg1 - mg_c) * 0.4)    
    return mg, mg_c, cr1, cr2, cr3, cr4

countfuv_abs = np.vectorize(countfuv_abs)

def countnuv(mg):
    silica = 1.0
    b4 = 0.22
    b13 = 0.27
    b15 = 0.074
    n2 = 0.055
    mg1 = 20.0

    if 9.323 <= mg <= 15.0:
        mg_c = 2.634 + (26.316 * mg - 245.329) ** 0.5
    elif mg < 9.323:
        mg_c = 2.634
    else: 
        mg_c = mg 

    cr1 =  silica * 10.0 ** ((mg1 - mg_c) * 0.4)
    cr2 =  b4 * 10 ** ((mg1-mg_c) * 0.4)    
    cr3 =  b13 * 10 ** ((mg1 - mg_c) * 0.4)    
    cr4 =  b15 * 10 ** ((mg1 - mg_c) * 0.4)    
    cr5 =  n2 * 10 ** ((mg1 - mg_c) * 0.4)    
    return mg, mg_c, cr1, cr2, cr3, cr4, cr5

countnuv = np.vectorize(countnuv)

# Functions to TD1 flux to UVIT count rates.
flux_norm = 2E-13
def td1_countnuv(flux):
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

def td1_countfuv(flux):
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

# Function to detect blobs, estimate fluxes, and do some more. 
def im_flux(int_map):

    # To read data.
    try:
        fitsf = fits.open(int_map)[0].data
    except IOError:
        incomplete_fits = 'Incomplete FITS file. Check if Galex Servers are working properly.'
        print('\n{}\n'.format(incomplete_fits))
        error_file.write(incomplete_fits)
        sys.exit(1)

    # To convert RA & DEC to pixel coordinates.
    w = WCS(int_map)
    cor = w.all_world2pix(cc.ra.degree, cc.dec.degree, 1)
    try:
        selcen = [int(round(x)) for x in cor]
    except ValueError:
        ra_dec_outside = 'The provided RA DEC values fell outside the GALEX image.' 
        print('\n{}\n'.format(filter_confusion))
        error_file.write(filter_confusion)
        sys.exit(1)

    # The pixel scale of GALEX taken here is 1.5 arcsec/pixel.
    # To select a rectangular region centered on the provided positions.
    xfi = int(selcen[0]) - field_radius_im[instrument]
    yfi = int(selcen[1]) - field_radius_im[instrument]
    xse = int(selcen[0]) + field_radius_im[instrument]
    yse = int(selcen[1]) + field_radius_im[instrument]
    xfi = max(0, xfi)
    yfi = max(0, yfi)
    fitsf = fitsf[yfi:yse, xfi:xse]

    # To put a circular mask on the image.
    yo,xo = np.ogrid[-field_radius_im[instrument]: field_radius_im[instrument],
                     -field_radius_im[instrument]: field_radius_im[instrument]]

    mask = xo*xo + yo*yo > field_radius_im[instrument] ** 2
   
    # To see if mask needs to be reduced in size.
    mask_xshape, mask_yshape = mask.shape
    fitsf_xshape, fitsf_yshape = fitsf.shape
    x_mask_start = mask_xshape - fitsf_xshape
    y_mask_start = mask_yshape - fitsf_yshape   
    mask = mask[x_mask_start:, y_mask_start:]   
    fitsf[mask] = 0
    
    # To mark positions of bright objects
    det_pos = w.all_world2pix(m_ra, m_dec, 1)
    x_det = det_pos[0] - xfi
    y_det = det_pos[1] - yfi

    # Using a 7x7 box to estimate the flux of detected objects.
    # the units are counts/sec/pixel.
    fluxes = []
    for blob in zip(y_det, x_det): 
        ydet, xdet = blob
        xfirst = int(xdet) - 3
        yfirst = int(ydet) - 3
        xsecond = int(xdet) + 4
        ysecond = int(ydet) + 4
        fluxes.append(fitsf[yfirst:ysecond, xfirst:xsecond].sum())

    # To plot the image.
    plt.imshow(fitsf,
               cmap = 'gray',
               norm = LogNorm(),
               interpolation = 'none')

    plt.title('Detected bright sources marked')
    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.yticks([])

    # To plot detected sources.
    plt.scatter(x_det, y_det, 
                color = 'r',
                marker = 'o',
                alpha = 0.2)

    # To annotate positions.
    anno = np.arange(len(x_det)) + 1
    for q, txt in enumerate(anno):
        plt.annotate(txt, (x_det[q], y_det[q]))
    
    # To save the image.
    figure_name = int_map.replace('.fits.gz', '.png')
    plt.savefig(figure_name,
                format = 'png',
                bbox_inches = 'tight',
                dpi = 300)

    plt.clf()  
    return fluxes

# Function to do all the work on TD1_catalogue.
def td1_estimate(td1_hdu):
    alpha = td1_hdu[1].data['ra']
    delta = td1_hdu[1].data['dec']
    nuv_flux = td1_hdu[1].data['flux_2365_a']
    fuv_flux = td1_hdu[1].data['flux_1565_a']
    
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
    flux, ta, tb, tc, td, te = td1_countnuv(nd)
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
    flux, ta, tb, tc, td = td1_countfuv(fd)
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

# Function to convert ra_deg and dec_deg to ra_hms and dec_dms.
def deg_to_hms(al, dl):
    fuv_coord = SkyCoord(zip(al, dl),
                     frame = 'icrs',
                     unit = 'deg')
    
    RAhms_DECdms = fuv_coord.to_string('hmsdms', sep = ':')
    ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])
    sl_no = np.arange(len(ra_hms)) + 1
    xy_tab = Table([sl_no, ra_hms, dec_dms], names = ('sl_no', 'ra_hms', 'dec_dms'))
    return xy_tab

# Function to format NUV data.
def format_nuv(nuv_counts):
    ntab = Table(nuv_counts,
                 names = ('Mag',
                          'Mag_corrected', 
                          'silica',
                          'b4',
                          'b13',
                          'b15',
                          'n2'), 
                 meta = {'name': 'NUV counts'})

    ntab['Mag'].format = '4.2f'
    ntab['Mag_corrected'].format = '4.2f'
    ntab['silica'].format = '4.2f'
    ntab['b4'].format = '4.2f'
    ntab['b13'].format = '4.2f'
    ntab['b15'].format = '4.2f'
    ntab['n2'].format = '4.2f'
    return ntab

# Function to format FUV data.
def format_fuv(fuv_counts):
    ftab = Table(fuv_counts,
                 names = ('Mag',
                          'Mag_corrected',
                          'caf2',
                          'baf2',
                          'sapphire',
                          'silica'), 
                 meta = {'name': 'FUV counts'})

    ftab['Mag'].format = '4.2f'
    ftab['Mag_corrected'].format = '4.2f'
    ftab['caf2'].format = '4.2f'
    ftab['baf2'].format = '4.2f'
    ftab['sapphire'].format = '4.2f'
    ftab['silica'].format = '4.2f'
    return ftab


# To check if Galactic latitude is between -30 to 30.
gal_lat = cc.galactic.b.value
gal_plane = 'no'
if -30.0 <= gal_lat <= 30.0:
    gal_plane = 'yes'

# To get on with MAST website queries.
ra = string.replace(RA, ':', '+')
dec = string.replace(DEC, ':', '+')

# Mast website form data.
mastdata = {'__EVENTTARGET': '""',
            '__EVENTARGUMENT' : '""',
            '__VIEWSTATE' : '/wEPDwUKMTUwNjg2NDc5Ng8WAh4TVmFsaWRhdGVSZXF1ZXN0TW9kZQIBFgQCAQ8WAh4JaW5uZXJodG1sBRNNQVNULkdhbGV4LlRpbGVMaXN0ZAIDD2QWAgIBDxYKHgtjZWxsc3BhY2luZwUBMB4LY2VsbHBhZGRpbmcFATAeBXdpZHRoBQM3NjAeBmJvcmRlcgUBMB4FYWxpZ24FBmNlbnRlchYGZg9kFgJmD2QWAmYPZBYOAgEPPCsABQEDFCsAARAWCB4GSXRlbUlEBRZfY3RsMl9NQVNULW1lbnVJdGVtMDAxHghJdGVtVGV4dAUETUFTVB4HSXRlbVVSTAUYaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1HhBNZW51SXRlbUNzc0NsYXNzBQt0b3BuYXZjb2xvcmRkZAIDDzwrAAUBAxQrAAEQFggfBwUXX2N0bDJfU1RTY0ktbWVudUl0ZW0wMDEfCAUFU1RTY0kfCQUUaHR0cDovL3d3dy5zdHNjaS5lZHUfCgULdG9wbmF2Y29sb3JkZGQCBQ88KwAFAQMUKwABEBYKHwcFIF9jdGwyX1NlYXJjaGVzX1Rvb2xzLW1lbnVJdGVtMDAxHwgFBVRvb2xzHwkFJmh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9zZWFyY2hlcy5odG1sHg1JdGVtTGVmdEltYWdlBREuLi9NZW51cy9kb3duLmdpZh4SSXRlbUxlZnRJbWFnZUFsaWduCyokU3lzdGVtLldlYi5VSS5XZWJDb250cm9scy5JbWFnZUFsaWduAhQrAAkQFgYfBwU0X2N0bDJfU2VhcmNoZXNfVG9vbHMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMB8IBQZBbGFkaW4fCQU5aHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2NnaS1iaW4vbnBoLWFsYWRpbi5wbD9mcm9tPVNUU2NJZGQQFgYfBwU0X2N0bDJfU2VhcmNoZXNfVG9vbHMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMR8IBQlTY3JhcGJvb2sfCQUmaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L3NjcmFwYm9vay5waHBkZBAWBh8HBTRfY3RsMl9TZWFyY2hlc19Ub29scy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDAyHwgFEFZpemllUi9NQVNUIFhjb3IfCQUjaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L3Zpemllci5waHBkZBAWBh8HBTRfY3RsMl9TZWFyY2hlc19Ub29scy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDAzHwgFDU5FRC9NQVNUIFhjb3IfCQUgaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L25lZC5waHBkZBAWBh8HBTRfY3RsMl9TZWFyY2hlc19Ub29scy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDA0HwgFCUNvcGxvdHRlch8JBSlodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvbWFzdF9jb3Bsb3QuaHRtbGRkEBYGHwcFNF9jdGwyX1NlYXJjaGVzX1Rvb2xzLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDUfCAUIU3BlY3ZpZXcfCQU5aHR0cDovL3d3dy5zdHNjaS5lZHUvcmVzb3VyY2VzL3NvZnR3YXJlX2hhcmR3YXJlL3NwZWN2aWV3ZGQQFgYfBwU0X2N0bDJfU2VhcmNoZXNfVG9vbHMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwNh8IBQhTdGFyVmlldx8JBR9odHRwOi8vc3RhcnZpZXcuc3RzY2kuZWR1L2h0bWwvZGQQFgYfBwU0X2N0bDJfU2VhcmNoZXNfVG9vbHMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwNx8IBQlBYnN0cmFjdHMfCQUnaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2Fic3RyYWN0cy5odG1sZGQQFgYfBwU0X2N0bDJfU2VhcmNoZXNfVG9vbHMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwOB8IBQdtb3JlLi4uHwkFJmh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9zZWFyY2hlcy5odG1sZGRkZAIHDzwrAAUBAxQrAAEQFgofBwUaX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEfCAUOTWlzc2lvbiBTZWFyY2gfCQUmaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L21pc3Npb25zLmh0bWwfCwURLi4vTWVudXMvZG93bi5naWYfDAsrBAIUKwAdEBYGHwcFLl9jdGwyX01pc3Npb25zLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDAfCAURIDxiPiBIdWJibGUgPC9iPiAfCQUnaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2hzdC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMR8IBSAgPGI+IEh1YmJsZSBMZWdhY3kgQXJjaGl2ZSA8L2I+IB8JBSFodHRwOi8vaGxhLnN0c2NpLmVkdS9obGF2aWV3Lmh0bWxkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDAyHwgFFCA8Yj4gSFNUb25saW5lIDwvYj4gHwkFLWh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9oc3RvbmxpbmUvc2VhcmNoLnBocGRkEBYGHwcFLl9jdGwyX01pc3Npb25zLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDMfCAUjIDxiPiBIU1QgUHJlc3MgUmVsZWFzZSBJbWFnZXMgPC9iPiAfCQUoaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L3N0cHIvc2VhcmNoLnBocGRkEBYGHwcFLl9jdGwyX01pc3Npb25zLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDQfCAUOIDxiPiBEU1MgIDwvYj4fCQUqaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2NnaS1iaW4vZHNzX2Zvcm0vZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwNR8IBRQgPGI+IEdBTEVYVmlldyAgPC9iPh8JBQsvR2FsZXhWaWV3L2RkEBYGHwcFLl9jdGwyX01pc3Npb25zLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDYfCAUQIDxiPiBHQUxFWCAgPC9iPh8JBRMvR1I2Lz9wYWdlPW1hc3Rmb3JtZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwNx8IBRsgPGI+IEpXU1QgU0lEIEFyY2hpdmUgIDwvYj4fCQUzaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2p3c3Qvc2lkYXJjaGl2ZS9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwOB8IBRUgPGI+IEtlcGxlciBEYXRhIDwvYj4fCQU2aHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2tlcGxlci9kYXRhX3NlYXJjaC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwOR8IBRggPGI+IEtlcGxlciBUYXJnZXRzIDwvYj4fCQU1aHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2tlcGxlci9rZXBsZXJfZm92L3NlYXJjaC5waHBkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDEwHwgFEyA8Yj4gU3dpZnRVVk9UIDwvYj4fCQUtaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L3N3aWZ0dXZvdC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxMR8IBREgPGI+IFhNTS1PTSAgPC9iPh8JBSpodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUveG1tLW9tL3NlYXJjaC5waHBkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDEyHwgFDiBCRUZTIChPUkZFVVMpHwkFKGh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9iZWZzL3NlYXJjaC5waHBkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDEzHwgFDyBDb3Blcm5pY3VzLXJhdx8JBS5odHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvY29wZXJuaWN1cy9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxNB8IBREgQ29wZXJuaWN1cy1jb2FkZB8JBTRodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvY29wZXJuaWN1cy9jb2FkZC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxNR8IBQYgRVBPQ0gfCQU4aHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2Vwb2NoL2Vwb2NoX21hc3RfZGlyZWN0b3J5Lmh0bWxkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDE2HwgFBiBFVVZFIB8JBShodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvZXV2ZS9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxNx8IBRFGVVNFIE9ic2VydmF0aW9ucx8JBShodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvZnVzZS9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxOB8IBQ5GVVNFIEV4cG9zdXJlcx8JBTFodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvZnVzZS9leHBvc3VyZS9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAxOR8IBQUgR1NDIB8JBTdodHRwOi8vZ3Nzcy5zdHNjaS5lZHUvd2Vic2VydmljZXMvR1NDMi9HU0MyV2ViRm9ybS5hc3B4ZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyMB8IBQYgSFBPTCAfCQUoaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2hwb2wvc2VhcmNoLnBocGRkEBYGHwcFLl9jdGwyX01pc3Npb25zLW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMjEfCAUFIEhVVCAfCQUnaHR0cDovL2FyY2hpdmUuc3RzY2kuZWR1L2h1dC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyMh8IBRAgSU1BUFMgKE9SRkVVUykgHwkFKWh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9pbWFwcy9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyMx8IBQUgSVVFIB8JBSdodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvaXVlL3NlYXJjaC5waHBkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDI0HwgFDyBUVUVTIChPUkZFVVMpIB8JBShodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvdHVlcy9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyNR8IBQUgVUlUIB8JBSdodHRwOi8vYXJjaGl2ZS5zdHNjaS5lZHUvdWl0L3NlYXJjaC5waHBkZBAWBh8HBS5fY3RsMl9NaXNzaW9ucy1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDI2HwgFCyBWTEEtRklSU1QgHwkFLGh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS92bGFmaXJzdC9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyNx8IBQcgV1VQUEUgHwkFKWh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS93dXBwZS9zZWFyY2gucGhwZGQQFgYfBwUuX2N0bDJfTWlzc2lvbnMtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAyOB8IBQdtb3JlLi4uHwkFL2h0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS9zZWFyY2hlcy5odG1sI21pc3Npb25zZGRkZAIJDzwrAAUBAxQrAAEQFgYfBwUbX2N0bDJfVHV0b3JpYWxzLW1lbnVJdGVtMDAxHwgFCFR1dG9yaWFsHwkFLGh0dHA6Ly9hcmNoaXZlLnN0c2NpLmVkdS90dXRvcmlhbC9pbmRleC5odG1sZGRkAgsPPCsABQEDFCsAARAWCB8HBRxfY3RsMl9TaXRlU2VhcmNoLW1lbnVJdGVtMDAxHwgFC1NpdGUgU2VhcmNoHwkFEi4vP3BhZ2U9c2l0ZXNlYXJjaB8KBQt0b3BuYXZjb2xvcmRkZAINDzwrAAUBAxQrAAEQFggfBwUaX2N0bDJfRm9sbG93VXMtbWVudUl0ZW0wMDAfCAUJRm9sbG93IFVzHwsFES4uL01lbnVzL2Rvd24uZ2lmHwwLKwQCFCsAAhAWBh8HBS5fY3RsMl9Gb2xsb3dVcy1tZW51SXRlbTAwMC1zdWJNZW51LW1lbnVJdGVtMDAwHwgFCiBGYWNlYm9vayAfCQUjaHR0cDovL3d3dy5mYWNlYm9vay5jb20vTUFTVEFyY2hpdmVkZBAWBh8HBS5fY3RsMl9Gb2xsb3dVcy1tZW51SXRlbTAwMC1zdWJNZW51LW1lbnVJdGVtMDAxHwgFCSBUd2l0dGVyIB8JBR5odHRwczovL3R3aXR0ZXIuY29tL01BU1RfTmV3cy9kZGRkAgIPZBYEZg9kFgJmD2QWAgIBDzwrAAUBAxQrAAoQFggfBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEfCAUSU2VhcmNoICYgUmV0cmlldmFsHwsFEi4uL01lbnVzL2Fycm93LmdpZh8MCysEAhQrAAMQFgYfBwUuX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMB8IBRVTb3VyY2UgQ2F0YWxvZyBTZWFyY2gfCQUQLi8/cGFnZT1tYXN0Zm9ybWRkEBYGHwcFLl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDEfCAUKU1FMIFNlYXJjaB8JBQ8uLz9wYWdlPXNxbGZvcm1kZBAWCB8HBS5fY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDAyHwgFC1RpbGUgU2VhcmNoHwsFEi4uL01lbnVzL2Fycm93LmdpZh8MCysEAhQrAAgQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMi1zdWJNZW51LW1lbnVJdGVtMDAwHwgFGjxiPkFJUzwvYj46IEFsbCBTa3kgU3VydmV5Hg9JdGVtQ29tbWFuZE5hbWUFA2Fpc2RkEBYGHwcFQl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDItc3ViTWVudS1tZW51SXRlbTAwMR8IBR88Yj5ESVM8L2I+OiBEZWVwIEltYWdpbmcgU3VydmV5Hw0FA2Rpc2RkEBYGHwcFQl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDItc3ViTWVudS1tZW51SXRlbTAwMh8IBSE8Yj5NSVM8L2I+OiBNZWRpdW0gSW1hZ2luZyBTdXJ2ZXkfDQUDbWlzZGQQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMi1zdWJNZW51LW1lbnVJdGVtMDAzHwgFIjxiPk5HUzwvYj46IE5lYXJieSBHYWxheGllcyBTdXJ2ZXkfDQUDbmdzZGQQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMi1zdWJNZW51LW1lbnVJdGVtMDA0HwgFIzxiPkdJSTwvYj46IEd1ZXN0IEludmVzdGlnYXRvciBEYXRhHw0FA2dpaWRkEBYGHwcFQl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDItc3ViTWVudS1tZW51SXRlbTAwNR8IBR48Yj5DQUk8L2I+OiBDYWxpYnJhdGlvbiBTdXJ2ZXkfDQUDY2FpZGQQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDEtc3ViTWVudS1tZW51SXRlbTAwMi1zdWJNZW51LW1lbnVJdGVtMDA2HwgFKDxiPlNQRUNUUkE8L2I+OiBGcm9tIEFsbCBBdmFpbGFibGUgVGlsZXMfDQUHc3BlY3RyYWRkEBYGHwcFQl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDItc3ViTWVudS1tZW51SXRlbTAwNx8IBRI8Yj5BTEwgU1VSVkVZUzwvYj4fDQUKYWxsc3VydmV5c2RkZGQQFgofBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDMfCAUTR3Vlc3QgSW52ZXN0aWdhdG9ycx8JBQ4uLz9wYWdlPWdpbGlzdB8LBRMuLi9NZW51cy9zZWN1cmUuZ2lmHwwLKwQCZGQQFggfBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDUfCAUNRG9jdW1lbnRhdGlvbh8LBRIuLi9NZW51cy9hcnJvdy5naWYfDAsrBAIUKwACEBYIHwcFLl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDA1LXN1Yk1lbnUtbWVudUl0ZW0wMDAfCAULPGI+TUFTVDwvYj4fCwUSLi4vTWVudXMvYXJyb3cuZ2lmHwwLKwQCFCsABRAWBh8HBUJfY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAwNS1zdWJNZW51LW1lbnVJdGVtMDAwLXN1Yk1lbnUtbWVudUl0ZW0wMDAfCAUKSGlnaCBMZXZlbB8JBRIuLz9wYWdlPWdlbmVyYWxmYXFkZBAWBh8HBUJfY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAwNS1zdWJNZW51LW1lbnVJdGVtMDAwLXN1Yk1lbnUtbWVudUl0ZW0wMDEfCAUHUXVlcmllcx8JBQ4uLz9wYWdlPXNxbGZhcWRkEBYGHwcFQl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDA1LXN1Yk1lbnUtbWVudUl0ZW0wMDAtc3ViTWVudS1tZW51SXRlbTAwMh8IBRBEYXRhIERlc2NyaXB0aW9uHwkFDS4vP3BhZ2U9ZGRmYXFkZBAWBh8HBUJfY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAwNS1zdWJNZW51LW1lbnVJdGVtMDAwLXN1Yk1lbnUtbWVudUl0ZW0wMDMfCAUOVXNlciBTdWJtaXR0ZWQfCQUPLi8/cGFnZT11c2VyZmFxZGQQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDUtc3ViTWVudS1tZW51SXRlbTAwMC1zdWJNZW51LW1lbnVJdGVtMDA1HwgFCFR1dG9yaWFsHwkFEC4vP3BhZ2U9dHV0b3JpYWxkZGQQFggfBwUuX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDUtc3ViTWVudS1tZW51SXRlbTAwMR8IBRM8Yj5DYWx0ZWNoIEZBUXM8L2I+HwsFEi4uL01lbnVzL2Fycm93LmdpZh8MCysEAhQrAAIQFgYfBwVCX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDUtc3ViTWVudS1tZW51SXRlbTAwMS1zdWJNZW51LW1lbnVJdGVtMDAwHwgFEENhbHRlY2ggTWV0YWRhdGEfCQULLi8/cGFnZT1mYXFkZBAWBh8HBUJfY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAwNS1zdWJNZW51LW1lbnVJdGVtMDAxLXN1Yk1lbnUtbWVudUl0ZW0wMDEfCAURQ2FsdGVjaCBUZWNoIERvY3MfCQU1aHR0cDovL3d3dy5nYWxleC5jYWx0ZWNoLmVkdS9yZXNlYXJjaGVyL3RlY2hkb2NzLmh0bWxkZGRkEBYGHwcFGl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDA3HwgFDURhdGFiYXNlIEluZm8fCQUMP3BhZ2U9ZGJpbmZvZGQQFgYfBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMDkfCAUUQ29udHJpYnV0ZWQgU29mdHdhcmUfCQUQLi8/cGFnZT1zb2Z0d2FyZWRkEBYIHwcFGl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDExHwgFF0d1ZXN0IEludmVzdGlnYXRvciBTaXRlHwsFEi4uL01lbnVzL2Fycm93LmdpZh8MCysEAhQrAAMQFgYfBwUuX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMTEtc3ViTWVudS1tZW51SXRlbTAwMB8IBQlIb21lIFBhZ2UfCQUdaHR0cDovL2dhbGV4Z2kuZ3NmYy5uYXNhLmdvdi9kZBAWBh8HBS5fY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAxMS1zdWJNZW51LW1lbnVJdGVtMDAxHwgFD0luc3RydW1lbnRhdGlvbh8JBUFodHRwOi8vZ2FsZXhnaS5nc2ZjLm5hc2EuZ292L0RvY3VtZW50cy9FUk9fZGF0YV9kZXNjcmlwdGlvbl8yLmh0bWRkEBYGHwcFLl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDExLXN1Yk1lbnUtbWVudUl0ZW0wMDIfCAUNRGF0YSBQaXBlbGluZR8JBUFodHRwOi8vZ2FsZXhnaS5nc2ZjLm5hc2EuZ292L0RvY3VtZW50cy9FUk9fZGF0YV9kZXNjcmlwdGlvbl8zLmh0bWRkZBAWBh8HBRpfY3RsOF9sZWZ0TWVudS1tZW51SXRlbTAxMx8IBQ1SZWxhdGVkIFNpdGVzHwkFFC4vP3BhZ2U9cmVsYXRlZHNpdGVzZGQQFgYfBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMTUfCAUPQWNrbm93bGVkZ21lbnRzHwkFFy4vP3BhZ2U9YWNrbm93bGVkZ21lbnRzZGQQFhQfCQULL0dhbGV4Vmlldy8eD01lbnVJdGVtVG9vbFRpcAUdR2FsZXhWaWV3IChRdWljayBTZWFyY2gpIFRvb2weCkl0ZW1UYXJnZXQFBl9ibGFuax4QSXRlbUltYWdlQWx0VGV4dAUdR2FsZXhWaWV3IChRdWljayBTZWFyY2gpIFRvb2wfCAUKR2FsZXhWaWV3Oh4NTWVudUl0ZW1XaWR0aBsAAAAAAMBiQAEAAAAfBwUaX2N0bDhfbGVmdE1lbnUtbWVudUl0ZW0wMTYeDkl0ZW1SaWdodEltYWdlBRlpbWFnZXMvR2FsZXhWaWV3VGh1bWIucG5nHhNJdGVtUmlnaHRJbWFnZUFsaWduCysEAR4OTWVudUl0ZW1IZWlnaHQbAAAAAADAYkABAAAAZGQQFhQfCQUJL2Nhc2pvYnMvHw4FG0Nhc0pvYnMgKERhdGFiYXNlIFNRTCkgVG9vbB8PBQZfYmxhbmsfEAUbQ2FzSm9icyAoRGF0YWJhc2UgU1FMKSBUb29sHwgFCENhc0pvYnM6HxEbAAAAAADAYkABAAAAHwcFGl9jdGw4X2xlZnRNZW51LW1lbnVJdGVtMDE3HxIFF2ltYWdlcy9DYXNKb2JzVGh1bWIucG5nHxMLKwQBHxQbAAAAAABAYEABAAAAZGRkAgEPZBYCZg8PFgQeCXNvcnRPcmRlcgUHcmFfY2VudB4Mc2hvd0FsbFRpbGVzBQVmYWxzZWQWBAIBDw8WAh4EVGV4dAXKAjxiPlRoZXJlIGFyZSA0NTE5NCB0b3RhbCB0aWxlcyBpbiBhbGwgdGhlIEdBTEVYIHN1cnZleXMuPC9iPjxicj48Zm9udCBzaXplPSctMScgY29sb3I9J2dyYXknPlBsZWFzZSBub3RlOiBTZWFyY2hlcyBpbiB0aGlzIHBhZ2UgYXBwbHkgb25seSB0byBUSUxFIGxldmVsIHByb2R1Y3RzLjxicj5JZiB5b3Ugd2FudCB0byBzZWFyY2ggR0FMRVggb2JqZWN0IGNhdGFsb2dzLCBwbGVhc2UgdXNlIGVpdGhlciB0aGUgPGEgaHJlZj0nP3BhZ2U9bWFzdGZvcm0nPkNhdGFsb2cgT2JqZWN0IFNlYXJjaDwvYT4gb3IgdGhlIDxhIGhyZWY9Jz9wYWdlPXNxbGZvcm0nPlNRTCBTZWFyY2g8L2E+LmRkAhUPDxYCHgdWaXNpYmxlaGQWAgIDDzwrAAsAZAIDD2QWAmYPZBYCZg9kFgICBQ8PFgIfFwUrTGFzdCBNb2RpZmllZCBEYXRlOjxicj4xMi81LzIwMTYgMTo1MTozOSBQTWRkZBOt2pbUX66uvUbSuy3q9kQU8fEC',
           '__VIEWSTATEGENERATOR': 'C84C2718',
           '__EVENTVALIDATION': '/wEdAAue+6xrb6xgp2ityzurA/pfWsTF2CBs9ziYHlDmus7EnHXVqisK/ch+FuYDN4RJj9bNygAwoalISibjyjYgoB7/Pb1PMsXU2LG7o+i6/zoft2ZmqVWZEJyWTGlJer/5/ymk9SeG9Y8RLbkbyiuf4BcRXP2SoyGCMZyu6LfyUjL5ZgAB13huDNxtBirRDFLR6zW3raPnQUy5sK21W/3eiEs/KUQOVtp9GallVy/IsFMIp4yMEruOYx0KrV7GUndYi0m5y40+',
           '_ctl10:txtTargetName': '',
           '_ctl10:resolverDropList': 'SIMBAD',
           '_ctl10:txtRadius': '0.001',
           '_ctl10:txtRA': ra,
           '_ctl10:txtDec': dec,
           '_ctl10:btnSearch': 'Search'}

# Header for the site.
header = {'Host': 'galex.stsci.edu',
          'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64; rv:45.0) Gecko/20100101 Firefox/45.0',
          'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
          'Accept-Language': 'en-US,en;q=0.5',
          'Accept-Encoding': 'gzip, deflate',
          'Referer': 'http://galex.stsci.edu/GR6/?page=tilelist&survey=allsurveys',
          'Connection': 'keep-alive'}

# Post request.
session = Session()
response = session.post(url = 'http://galex.stsci.edu/GR6/?page=tilelist&survey=allsurveys',
                        data = mastdata,
                        headers = header)

# To make sense of the mess that is MAST. 
soup = BeautifulSoup(response.text, 'html.parser')
try:
    ss = soup.find(id = '_ctl10_TileGrid_imgLink_0').get('href')
except AttributeError:
    no_galex_tiles = '0 Galex tiles found. Galex observations around \
                      \nthe given target is not available. Using TD1\
                      \ncatalogue to estimate UVIT count rates.'
    print('\n\n{}\n\n'.format(no_galex_tiles))
    with open('zero_tiles.txt', 'w') as tile_f:
        tile_f.write(no_galex_tiles)
    if gal_plane == 'yes':
        gal_plane_warning = 'The galactic latitude is between -30 to 30. \
                            \nYour field cannot be checked using TD1 catalogue!'
        print('\n\n{}\n\n'.format(gal_plane_warning))
        error_file.write(gal_plane_warning)
        sys.exit(1)
    else:
        td1_estimate(td1_hdu)
        sys.exit(1)

parent = 'http://galex.stsci.edu/GR6'
son = parent + ss[1:] 

# To make sense of the download page.
response2 = session.get(son)
soup2 = BeautifulSoup(response2.text, 'html.parser')
pp = soup2.find_all("a")

# To download the mcat (galex catalogue).
catalogue_link = []
for link in pp:
    somel = link.get('href')
    try:
        if somel[-12:] == 'mcat.fits.gz':
            catalogue_link.append(somel)
    except:
        pass

if len(catalogue_link) != 0:
    catalogue = catalogue_link[0].split('/')[-1]
    urllib.urlretrieve(catalogue_link[0], catalogue)
else:
    no_catalogue = 'Could not find the catalogue for this region.'
    print('\n{}\n'.format(no_catalogue))
    error_file.write(no_catalogue)
    sys.exit(1)

# Reading coordinates from catalogue.
try:
    hdu = fits.open(catalogue)
except IOError:
    incomplete_fits = 'Incomplete FITS file. Check if Galex Servers are working properly.'
    print('\n{}\n'.format(incomplete_fits))
    error_file.write(incomplete_fits)
    sys.exit(1)

alpha = hdu[1].data['alpha_j2000_merged']
delta = hdu[1].data['delta_j2000_merged']


# NUV 
nuv_mag = hdu[1].data['nuv_mag']
nuv_fwhm = hdu[1].data['nuv_fwhm_world']
refined_set = [(al, de, nm, nf) for al, de, nm, nf
                            in zip(alpha, delta, nuv_mag, nuv_fwhm)
                            if int(nm) != -999 and nm <= 22.]

nalpha, ndelta, nuv_mag, nuv_fwhm = zip(*refined_set)

confined_set = [(nm, al, de, nf) for al, de, nm, nf 
                             in zip(nalpha, ndelta, nuv_mag, nuv_fwhm) 
                             if cel_separation(al, de) <= field_radius[instrument] * u.arcsec]

nd = np.array(sorted(confined_set))[0:5]
cat_nuv_counts = countnuv(nd[:,0])
cat_nuv_res = format_nuv(cat_nuv_counts)

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
xy_tab = deg_to_hms(nd[:,1], nd[:,2])
cat_nuv_res = hstack([xy_tab, cat_nuv_res])

balance = cat_nuv_res['ra_hms', 'dec_dms','Mag']
balance.rename_column('Mag', 'CAT_Mag')
balance['fwhm'] = nd[:,3]
balance['fwhm'].format = '4.4f' 

# FUV 
fuv_mag = hdu[1].data['fuv_mag']

fuv_absent = 'no'
if len(np.unique(fuv_mag)) == 1:  # when FUV data is absent.
    fd = nd
    warning = '\nFUV observations seem to be absent! Using M_fuv = M_nuv - 1.65.'
    print(warning)
    fuv_absent = 'yes'
    with open('absent_FUV.txt', 'w') as the_file:
        the_file.write(warning)
else:
    refined_set = [(al, de, fm) for al, de, fm 
                                in zip(alpha, delta, fuv_mag)
                                if int(fm) != -999 and fm <= 22.]
    
    falpha, fdelta, fuv_mag = zip(*refined_set)
    
    confined_set = [(fm, al, de) for al, de, fm 
                                 in zip(falpha, fdelta, fuv_mag) 
                                 if cel_separation(al, de) <= field_radius[instrument] * u.arcsec]
    
    fd = np.array(sorted(confined_set))[0:5]

if fuv_absent == 'no':
    cat_fuv_counts = countfuv(fd[:,0])
else:
    cat_fuv_counts = countfuv_abs(fd[:,0])

cat_fuv_res = format_fuv(cat_fuv_counts)

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
xy_tab = deg_to_hms(fd[:,1], fd[:,2])
cat_fuv_res = hstack([xy_tab, cat_fuv_res])

# To download the galex images.
images = []
for link in pp:
    somel = link.get('href')
    try:
        if somel[-11:] == 'int.fits.gz':
            images.append(somel)
    except:
        pass

if len(images) != 0:
    for image in images:
        image_file = image.split('/')[-1]
        urllib.urlretrieve(image, image_file)
        detector = str(image_file[-14:-13]).lower()
        if detector == 'f':
            fuv_intmap = image_file
        elif detector == 'n':
            nuv_intmap = image_file
        else: 
            filter_confusion = 'Cannot determine filter! exiting.'
            print('\n{}\n'.format(filter_confusion))
            error_file.write(filter_confusion)
            sys.exit(1) 

# NUV
if 'nuv_intmap' in locals():
    m_ra = nd[:,1]
    m_dec = nd[:,2]
    n_fluxes = im_flux(nuv_intmap)
    n_mags = magnuv(n_fluxes)
    n_counts = countnuv(n_mags)
    im_nuv_res = format_nuv(n_counts)
    xy_tab = deg_to_hms(m_ra, m_dec)
    im_nuv_res = hstack([xy_tab, im_nuv_res])
    balance['IM_Mag'] = im_nuv_res['Mag']

# FUV 
if 'fuv_intmap' not in locals() and 'nuv_intmap' in locals():
    f_counts = countfuv_abs(n_mags)
    im_fuv_res = format_fuv(f_counts)
    im_fuv_res = hstack([xy_tab, im_fuv_res])
    
elif 'fuv_intmap' in locals():
    m_ra = fd[:,1]
    m_dec = fd[:,2]
    f_fluxes = im_flux(fuv_intmap)  
    f_mags = magfuv(f_fluxes)
    f_counts = countfuv(f_mags)
    im_fuv_res = format_fuv(f_counts)
    xy_tab = deg_to_hms(m_ra, m_dec)
    im_fuv_res = hstack([xy_tab, im_fuv_res])

# To decide between catalogue or image.
balance['diff'] = balance['IM_Mag'] - balance['CAT_Mag']

if sum(balance['diff'] > 1.2) == 0 and sum(balance['fwhm'] > 0.0043) == 0:
    nuv_res = cat_nuv_res
    fuv_res = cat_fuv_res
else:
    nuv_res = im_nuv_res
    fuv_res = im_fuv_res

# To select NUV safe filters.
print('\n\n### NUV\n\n{}\n'.format(nuv_res))
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
nuv_table = 'NUV_' + catalogue.replace('.fits.gz', '-nd-int.txt')
ascii.write(nuv_res, nuv_table, format = 'csv', overwrite = True)

with open('safe_NUV_filters.txt', 'w') as safe_file:
    safe_file.write(nuv_declaration)


# To select FUV safe filters.
print('\n### FUV \n\n{}\n\n'.format(fuv_res))
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
fuv_table = 'FUV_' + catalogue.replace('.fits.gz', '-fd-int.txt')
ascii.write(fuv_res, fuv_table, format = 'csv', overwrite = True)

with open('safe_FUV_filters.txt', 'w') as safe_file:
    safe_file.write(fuv_declaration)

print('Done!\n')





