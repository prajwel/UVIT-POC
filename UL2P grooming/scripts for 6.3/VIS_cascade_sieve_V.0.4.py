#!/usr/bin/env python3

'''This program is written to help arrange UL2P output as per ICD. 


   Copyright 2019 Prajwel Joseph
  
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

'''


import os
import re
import sys
import fnmatch
import tarfile
import numpy as np

from glob import glob
from astropy.io import fits
from shutil import copy, rmtree, move, copytree


dir_to_sieve = sys.argv[1]

if '/' in dir_to_sieve:
    print("\nplease provide the directory name without '/'\n")
    exit()

print("\nWorking on {}".format(dir_to_sieve))


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

# To take care of the multiple folders.
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

    asi = find('*as_Sig.fits', keep_dir)
#    Aexpi = find('*ra-decexp.fits', keep_dir)
    snri = find('*l2_radec.fits', keep_dir)
    sigi = find('*sig_regAvg.fits', keep_dir)
    Iexpi = find('*exp_regAvg.fits', keep_dir)
    Ierri = find('*noiseMap_regAvg.fits', keep_dir)

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
          
def sortout_files(uvchannel_dir):
    # To find and rename the as, sig, exp, and snr files
    os.chdir(uvchannel_dir)
    asi = find('*as_Sig.fits', '.')
    sigi = find('*sig_regAvg.fits', '.')

    if len(glob('*fits')) > 0:
        rmtree(currdir + '/' + new_dir + '/' + uvchannel_dir)
        os.mkdir(currdir + '/' + new_dir + '/' + uvchannel_dir)
    
    if len(asi) == 0 and len(sigi) > 0:
        print('No frame with astrometry in {}'.format(uvchannel_dir))
    
    sieve = True
    if len(sigi) == 0:
        sieve = False
        uvdir = glob('uvit')
        if len(uvdir) != 0:
            print('directory {} is empty'.format(uvchannel_dir))
            rmtree(uvdir[0])
    
    if len(sigi) > 0:
        selected_dir = largest_exposure(sigi)
        print("The selected directory is {}".format(selected_dir))
        rmtree(currdir + '/' + new_dir + '/' + uvchannel_dir + selected_dir[1:])

    ass = glob('*as_Sig.fits')
    snr = glob('*l2_radec.fits')
    sig = glob('*sig_regAvg.fits')
    Iexp = glob('*exp_regAvg.fits')
    Ierr = glob('*noiseMap_regAvg.fits')
    ras = glob('*dr.fits')

    if len(ass) != 0:
        os.rename(ass[0], (ass[0][0:36] + 'A_l2wcs.fits'))
    if len(snr) != 0:
        os.rename(snr[0], (snr[0][0:36] + '_l2ce.fits'))
    if len(sig) != 0:
        os.rename(sig[0], (sig[0][0:36] + 'I_l2sig.fits'))
    if len(Iexp) != 0:
        os.rename(Iexp[0], (Iexp[0][0:36] + 'I_l2exp.fits'))
    if len(Ierr) != 0:
        os.rename(Ierr[0], (Ierr[0][0:36] + 'I_l2err.fits'))
    if len(ras) != 0:
        os.rename(ras[0], (ras[0][0:36] + '_l2dr.fits')) 

    ras = glob('*dr.fits')
    VD = '../V' + uvchannel_dir[1:]
    try:
        os.mkdir(VD)
    except FileExistsError:
        pass
    try:
        copy(ras[0], VD)
        os.remove(ras[0])
    except IndexError:
        pass

    return sieve
    

if len(finD('uvit', dir_to_sieve)) == 0:
    print("This directory does not require sieving")
    exit()

if dir_to_sieve == 'uvt_ci':
    print("uvit_ci does not require sieving")
    exit()


# To make a dictionary of all the VIS RAS files.
rasdlist = find('*_dr.fits', '../')
ras_dict = dict([(rasd[-75: -68], rasd) for rasd in rasdlist])


do_sieve = True
while do_sieve == True:
    all_data = glob('uvt_[0-9][0-9]')
    uvt_dir_numbers = [i[-2:] for i in all_data]
    if int(max(uvt_dir_numbers)) > 8:
        suffix = str(int(max(uvt_dir_numbers)) + 1)
    else:
        suffix = '0' + str(int(max(uvt_dir_numbers)) + 1)
    new_dir = 'uvt_' + suffix
    copytree(dir_to_sieve, new_dir)
    currdir = os.getcwd()
    os.chdir(currdir + '/' + dir_to_sieve)


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

    fuv_do_sieve = sortout_files(fuvd)
    os.chdir(currdir + '/' + dir_to_sieve)
    nuv_do_sieve = sortout_files(nuvd)
    os.chdir(currdir + '/' + new_dir)

    os.rename(fuvd, fuvd[:2] + suffix)
    os.rename(nuvd, nuvd[:2] + suffix)
    os.rename('V_' + fuvd[-2:], 'V_' + suffix)

    if (fuv_do_sieve == False) and (nuv_do_sieve == False):
        do_sieve = False
    else:
        do_sieve = True

    if do_sieve == True:
        print("\nRemainder copied to {}\n".format(new_dir))
        dir_to_sieve = new_dir

    os.chdir(currdir)



if do_sieve == False:
    rmtree(dir_to_sieve)
    rmtree(new_dir)
    print("\n{} have been deleted\n".format(dir_to_sieve))
    













