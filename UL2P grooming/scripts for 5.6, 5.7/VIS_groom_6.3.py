#!/usr/bin/env python3


import re
import os
import fnmatch
import tarfile

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


# To put into place the driver module param file.
L1_filename = glob('LEVL1AS1UVT*tar_V*.2*')[0]
L1_datepid = list(L1_filename[:26])
L1_datepid[4] = '2'
ICD_driver_name = "".join(L1_datepid) + '_L2_DM_params.txt'
os.rename('UVIT_DriverModule.par', ICD_driver_name)
os.mkdir('pipeline')
move(ICD_driver_name, './pipeline/')

# To find output FUV & NUV directories.
fuvdlist = glob('output_FUV_[0-9]*')
nuvdlist = glob('output_NUV_[0-9]*')
if len(fuvdlist) == 0 and len(fuvdlist) == 0:
    print('\nCheck if output files in the format \"output_FUV*\" are present\n')

# To make a dictionary of all the VIS RAS files.
rasdlist = find('*dr.fits', '.')
ras_dict = dict([(rasd[-75: -68], rasd[1:]) for rasd in rasdlist])

# To find and rename the orbit-wise FUV, NUV files
currdir = os.getcwd()
for fuvd in fuvdlist:
    os.chdir(currdir + '/' + fuvd)
    Fasi = find('*as_Sig.fits', '.')
#    FAexpi = find('*ra-decexp.fits', '.')
    Fsnri = find('*l2_radec.fits', '.')
    Fsigi = find('*sig_regAvg.fits', '.')
    FIexpi = find('*exp_regAvg.fits', '.')
    FIerri = find('*noiseMap_regAvg.fits', '.')

    if len(Fsnri) == 1:
        move(Fsnri[0], '.')

    if len(FIexpi) == 1:
        move(FIexpi[0], '.')

#    if len(FAexpi) == 1:
#        move(FAexpi[0], '.')

    if len(Fasi) == 1:
        move(Fasi[0], '.')

    if len(FIerri) == 1:
        move(FIerri[0], '.')

    if len(Fsigi) == 1:
        move(Fsigi[0], '.')
        fuvdir = glob('uvit')
        if len(fuvdir) != 0:
            rmtree(fuvdir[0])

    if len(Fasi) == 0 and len(Fsigi) > 0:
        print('No frame with astrometry in {}'.format(fuvd))

    if len(Fsigi) == 0:
        fuvdir = glob('uvit')
        if len(fuvdir) != 0:
            print('directory {} is empty'.format(fuvd))
            rmtree(fuvdir[0])
            continue

    if len(Fsigi) > 1:
        print('\nExists more than one FUV image inside {}'.format(fuvd))
        print('Formatting not carried out.\n')
        continue

    Fas = glob('*as_Sig.fits')
#    FAexp = glob('*ra-decexp.fits')
    Fsnr = glob('*l2_radec.fits')
    Fsig = glob('*sig_regAvg.fits')
    FIexp = glob('*exp_regAvg.fits')
    FIerr = glob('*noiseMap_regAvg.fits')

    # To get the corresponding RAS file.
    if len(Fsnr) == 1: 
        Fhdu = fits.open(Fsnr[0])
        history = Fhdu[0].header['history']
        re_result = re.search('V/uvtV\.\d{2}/uvtC', str(history))

        if type(re_result.group()) is str:    
            ras_num = re_result.group()[2:9]
        else:
            print('\nmultiple RAS?! Big problem!')
            exit()
    
        ras_file = currdir + ras_dict[ras_num]
        copy(ras_file, '.')
    else:
        print('\nFUV events file is absent at {}!\n'.format(fuvd))
        continue

    ras = glob('*dr.fits')
    
    if len(Fas) != 0:
        os.rename(Fas[0], (Fas[0][0:36] + 'A_l2wcs.fits'))
#    if len(FAexp) != 0:
#        os.rename(FAexp[0], (FAexp[0][0:36] + 'A_l2exp.fits'))
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


os.chdir(currdir)
for nuvd in nuvdlist:
    os.chdir(currdir + '/' + nuvd)
    Nasi = find('*as_Sig.fits', '.')
#    NAexpi = find('*ra-decexp.fits', '.')
    Nsnri = find('*l2_radec.fits', '.')
    Nsigi = find('*sig_regAvg.fits', '.')
    NIexpi = find('*exp_regAvg.fits', '.')
    NIerri = find('*noiseMap_regAvg.fits', '.')

    if len(Nsnri) == 1:
        move(Nsnri[0], '.')

    if len(NIexpi) == 1:
        move(NIexpi[0], '.')

#    if len(NAexpi) == 1:
#        move(NAexpi[0], '.')

    if len(Nasi) == 1:
        move(Nasi[0], '.')

    if len(NIerri) == 1:
        move(NIerri[0], '.')

    if len(Nsigi) == 1:
        move(Nsigi[0], '.')
        nuvdir = glob('uvit')
        if len(nuvdir) != 0:
            rmtree(nuvdir[0])

    if len(Nasi) == 0 and len(Nsigi) > 0:
        print('No frame with astrometry in {}'.format(nuvd))

    if len(Nsigi) == 0:
        nuvdir = glob('uvit')
        if len(nuvdir) != 0:
            print('directory {} is empty'.format(nuvd))
            rmtree(nuvdir[0])
            continue

    if len(Nsigi) > 1:
        print('\nExists more than one NUV image inside {}'.format(nuvd))
        print('Formatting not carried out.\n')
        continue

    Nas = glob('*as_Sig.fits')
#    NAexp = glob('*ra-decexp.fits')
    Nsnr = glob('*l2_radec.fits')
    Nsig = glob('*sig_regAvg.fits')
    NIexp = glob('*exp_regAvg.fits')
    NIerr = glob('*noiseMap_regAvg.fits')

    # To get the corresponding RAS file.
    if len(Nsnr) == 1: 
        Nhdu = fits.open(Nsnr[0])
        history = Nhdu[0].header['history']
        re_result = re.search('V/uvtV\.\d{2}/uvtC', str(history))

        if type(re_result.group()) is str:    
            ras_num = re_result.group()[2:9]
        else:
            print('\nmultiple RAS?! Big problem!')
            exit()
    
        ras_file = currdir + ras_dict[ras_num]
        copy(ras_file, '.')
    else:
        print('\nNUV events file is absent at {}!\n'.format(nuvd))
        continue

    ras = glob('*dr.fits')

    if len(Nas) != 0:
        os.rename(Nas[0], (Nas[0][0:36] + 'A_l2wcs.fits'))
#    if len(NAexp) != 0:
#        os.rename(NAexp[0], (NAexp[0][0:36] + 'A_l2exp.fits'))
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

os.chdir(currdir)

#renaming the FUV & NUV folders. 
for fuvd in fuvdlist:
    try:
        t = int(fuvd[-2:])
        os.rename(fuvd, "F_" + str(t))
    except ValueError:
        os.rename(fuvd, "F_0" + fuvd[-1:])

for nuvd in nuvdlist:
    try:
        t= int(nuvd[-2:])
        os.rename(nuvd, "N_" + str(t))
    except ValueError:
        os.rename(nuvd, "N_0" + nuvd[-1:])

#Creating structure as per ICD
fcount = max(len(glob("F_*[0-9]")), len(glob("N_*[0-9]")))

for fc in range(1, fcount+1):
    if fc < 10:
        try:
            uvtD = 'uvt_0' + str(fc)
            VD = 'V_0' + str(fc)
            os.mkdir(uvtD)
            os.mkdir(VD)
        except OSError:
            print("Check user permissions")

        try:
            move("F_0"+str(fc), uvtD)
        except IOError:
            pass

        try:
            move("N_0"+str(fc), uvtD)
        except IOError:
            pass

        for dras in find('*dr.fits', './' + uvtD):
            copy(dras, VD)
            os.remove(dras)

        move(VD, uvtD)            

    else:
        try:
            uvtD = 'uvt_' + str(fc)
            VD = 'V_' + str(fc)
            os.mkdir(uvtD)
            os.mkdir(VD)
        except OSError:
            print("Check user permissions")

        try:
            move("F_"+str(fc), uvtD)
        except IOError:
            pass

        try:
            move("N_"+str(fc), uvtD)
        except IOError:
            pass

        for dras in find('*dr.fits', './' + uvtD):
            copy(dras, VD)
            os.remove(dras)

        move(VD, uvtD)            

# To find, move, and rename the combined FUV, NUV files.
os.mkdir('uvt_ci/')
os.mkdir('uvt_ci/F_ci/')
os.mkdir('uvt_ci/N_ci/')
for ciD in finD('*outputFUV*', '.'):
    os.chdir(currdir + ciD[1:])
    cis = find('*.fits', '.')
    for ci in cis:
        move(ci, currdir + '/uvt_ci/F_ci/')

os.chdir(currdir)

for ciD in finD('*outputNUV*', '.'):
    os.chdir(currdir + ciD[1:])
    cis = find('*.fits', '.')
    for ci in cis:
        move(ci, currdir + '/uvt_ci/N_ci/')

os.chdir(currdir + '/uvt_ci/')

for ci_as in find('*as_Sig.fits', '.'):
    as_hdu = fits.open(ci_as)
    prefx = as_hdu[0].header['NAMEPRFX']
    os.rename(ci_as, ci_as[:7] + prefx[:36] + 'A_l2wcs.fits')

for ci_as_exp in find('*as_Exp.fits', '.'):
    os.remove(ci_as_exp)

for ci_as_err in find('*as_NoiseMap.fits', '.'):
    os.remove(ci_as_err)

for ci_sig in find('*FinalImage_Sig.fits', '.'):
    sig_hdu = fits.open(ci_sig)
    prefx = sig_hdu[0].header['NAMEPRFX']
    os.rename(ci_sig, ci_sig[:7] + prefx[:36] + 'A_l2sig.fits')

for ci_exp in find('*FinalImage_Exp.fits', '.'):
    exp_hdu = fits.open(ci_exp)
    prefx = exp_hdu[0].header['NAMEPRFX']
    os.rename(ci_exp, ci_exp[:7] + prefx[:36] + 'A_l2exp.fits')

for ci_err in find('*FinalImage_NoiseMap.fits', '.'):
    err_hdu = fits.open(ci_err)
    prefx = err_hdu[0].header['NAMEPRFX']
    os.rename(ci_err, ci_err[:7] + prefx[:36] + 'A_l2err.fits')             

print("#######\nFini!\n#######")

os.chdir(currdir)
rm1dir = glob("outputFUV*")
rm2dir = glob("outputNUV*")
for rm1 in rm1dir:
    rmtree(rm1)
for rm2 in rm2dir:
    rmtree(rm2)












