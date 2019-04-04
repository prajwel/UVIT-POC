#!/usr/bin/env python3

# Make sure you've given the code all the permissions using chmod

import os
import re
import fnmatch
import tarfile
import numpy as np

from glob import glob
from astropy.io import fits
from shutil import copy, rmtree, move


# Apparently, this is the way to find things in python
def finD(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in dirs:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

# To take care of the multiple folders arising out of split in VIS channel.
def largest_exposure(list_of_images):
    exp_value = 0
    for image in list_of_images:
        where_to_look = os.path.dirname(os.path.dirname(image))
        exposure_map = find('*exp_regAvg.fits', where_to_look)
        if len(exposure_map) == 1:
            HDU_array = fits.open(exposure_map[0])[0].data
            median_value = np.median(HDU_array[HDU_array > 0])
            print('{} has median value of exposure = {}'.format(exposure_map[0], median_value)) 
            if median_value > exp_value:
                exp_value = median_value
                keep_dir = where_to_look
        elif len(exposure_map) == 0:
            print("\n No exposure map inside {}".format(where_to_look))
        else:
            print("\nMore than one exposure map inside the split folders!\n")

    asi = find('*as.fits', keep_dir)
#    Aexpi = find('*ra-decexp.fits', keep_dir)
    snri = find('*l2_radec.fits', keep_dir)
    sigi = find('*sig_regAvg.fits', keep_dir)
    Iexpi = find('*exp_regAvg.fits', keep_dir)
    Ierri = find('*noise_map_sig.fits', keep_dir)

    if len(snri) == 1:
        move(snri[0], '.')

    if len(Iexpi) == 1:
        move(Iexpi[0], '.')

#    if len(Aexpi) == 1:
#        move(Aexpi[0], '.')

    if len(asi) == 1:
        move(asi[0], '.')

    if len(Ierri) == 1:
        move(Ierri[0], '.')

    if len(sigi) == 1:
        move(sigi[0], '.')
        dat_dir = glob('uvit')
        if len(dat_dir) != 0:
            rmtree(dat_dir[0])

    snr = glob('*l2_radec.fits')

    # To get the corresponding RAS file.
    if len(sigi) == 1: 
        hdu = fits.open(snr[0])
        history = hdu[0].header['history']
        re_result = re.search('V/uvtV\.\d{2}/uvtC', str(history))

        if type(re_result.group()) is str:    
            ras_num = re_result.group()[2:9]
        else:
            print('\nmultiple RAS?! Big problem!')
            exit()
    
        ras_file = currdir + '/' + ras_dict[ras_num]
        copy(ras_file, '.')

    return keep_dir
          



#To find output FUV & NUV directories.
try:
    fuvd = glob('F_[0-9]*')[0]
except IndexError:
    print("\nCheck if files in the format \"F_*\"  are present\n")
    exit()
try:
    nuvd = glob('N_[0-9]*')[0]
except IndexError:
    print("\nCheck if files in the format \"N_*\"  are present\n")
    exit()

# To make a dictionary of all the VIS RAS files.
rasdlist = find('*_dr.fits', '../')
ras_dict = dict([(rasd[-75: -68], rasd) for rasd in rasdlist])

# To find and rename the FUV as, sig, and snr files
currdir = os.getcwd()
os.chdir(currdir + '/' + fuvd)
Fasi = find('*as.fits', '.')
Fsigi = find('*sig_regAvg.fits', '.')

if len(Fasi) == 0 and len(Fsigi) > 0:
    print('No frame with astrometry in {}'.format(fuvd))

if len(Fsigi) == 0:
    fuvdir = glob('uvit')
    if len(fuvdir) != 0:
        print('directory {} is empty'.format(fuvd))
        rmtree(fuvdir[0])

if len(Fsigi) > 0:
    selected_dir = largest_exposure(Fsigi)
    print("The selected directory is {}".format(selected_dir))

Fas = glob('*as.fits')
Fsnr = glob('*l2_radec.fits')
Fsig = glob('*sig_regAvg.fits')
FIexp = glob('*exp_regAvg.fits')
FIerr = glob('*noise_map_sig.fits')
ras = glob('*dr.fits')

if len(Fas) != 0:
    os.rename(Fas[0], (Fas[0][0:36] + 'A_l2wcs.fits'))
if len(Fsnr) != 0:
    os.rename(Fsnr[0], (Fsnr[0][0:36] + '_l2ce.fits'))
if len(Fsig) != 0:
    os.rename(Fsig[0], (Fsig[0][0:36] + 'I_l2sig.fits'))
if len(FIexp) != 0:
    os.rename(FIexp[0], (FIexp[0][0:36] + 'I_l2exp.fits'))
if len(FIerr) != 0:
    os.rename(FIerr[0], (FIerr[0][0:36] + 'I_l2err.fits'))
if len(ras) != 0:
    os.rename(ras[0], (ras[0][0:36] + '_l2dr.fits')) 

ras = glob('*dr.fits')
VD = '../V' + fuvd[1:]
try:
    os.mkdir(VD)
except FileExistsError:
    pass
try:
    copy(ras[0], VD)
    os.remove(ras[0])
except IndexError:
    pass

# To find and rename the NUV as, sig, and snr files
os.chdir(currdir)
os.chdir(currdir + '/' + nuvd)
Nasi = find('*as.fits', '.')
Nsigi = find('*sig_regAvg.fits', '.')


if len(Nasi) == 0 and len(Nsigi) > 0:
    print('No frame with astrometry in {}'.format(nuvd))

if len(Nsigi) == 0:
    nuvdir = glob('uvit')
    if len(nuvdir) != 0:
        print('directory {} is empty'.format(nuvd))
        rmtree(nuvdir[0])

if len(Nsigi) > 0 :
    selected_dir = largest_exposure(Nsigi)
    print('\nThe selected directory is {}\n'.format(selected_dir))

Nas = glob('*as.fits')
Nsnr = glob('*l2_radec.fits')
Nsig = glob('*sig_regAvg.fits')
NIexp = glob('*exp_regAvg.fits')
NIerr = glob('*noise_map_sig.fits')
ras = glob('*dr.fits')

if len(Nas) != 0:
    os.rename(Nas[0], (Nas[0][0:36] + 'A_l2wcs.fits'))
if len(Nsnr) != 0:
    os.rename(Nsnr[0], (Nsnr[0][0:36] + '_l2ce.fits'))
if len(Nsig) != 0:
    os.rename(Nsig[0], (Nsig[0][0:36] + 'I_l2sig.fits'))
if len(NIexp) != 0:
    os.rename(NIexp[0], (NIexp[0][0:36] + 'I_l2exp.fits'))
if len(NIerr) != 0:
    os.rename(NIerr[0], (NIerr[0][0:36] + 'I_l2err.fits'))
if len(ras) != 0:
    os.rename(ras[0], (ras[0][0:36] + '_l2dr.fits')) 

ras = glob('*dr.fits')
VD = '../V' + nuvd[1:]
try:
    os.mkdir(VD)
except FileExistsError:
    pass
try:
    copy(ras[0], VD)
    os.remove(ras[0])
except IndexError:
    pass

os.chdir(currdir)

print("\nDone\n")














