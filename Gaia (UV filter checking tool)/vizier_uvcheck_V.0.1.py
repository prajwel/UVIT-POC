#!/usr/bin/env python3

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
'''


import os
import sys
import numpy as np

from astropy import units as u
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack


# To get the user input. 
#instrument = str(sys.argv[1])
#RA = str(sys.argv[2])
#DEC = str(sys.argv[3])

instrument = 'uvit'
RA = "7:36:51.396"
DEC = "65:36:9.170"


# Functions to convert magnitude to UVIT count rates.
def countfuv(mg):
    caf2 = 1.0
    baf2 = 0.85
    sapphire = 0.63
    silica = 0.22
    mg1 = 18.22

    if mg <= 15.0:
        mg_c = 5.371 + (20.0 * mg - 210.2) ** 0.5
    else: 
        mg_c = mg 

    cr1 =  caf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr2 =  baf2 * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr3 =  sapphire * 10.0 ** ((mg1 - mg_c) * 0.4)    
    cr4 =  silica * 10.0 ** ((mg1 - mg_c) * 0.4)    
    return mg, mg_c, cr1, cr2, cr3, cr4

countfuv = np.vectorize(countfuv)

def countnuv(mg):
    silica = 1.0
    b4 = 0.22
    b13 = 0.27
    b15 = 0.074
    n2 = 0.055
    mg1 = 20.0
    if mg <= 15.0:
        mg_c = 2.634 + (26.316 * mg - 245.329) ** 0.5
    else: 
        mg_c = mg 

    cr1 =  silica * 10.0 ** ((mg1 - mg_c) * 0.4)
    cr2 =  b4 * 10 ** ((mg1-mg_c) * 0.4)    
    cr3 =  b13 * 10 ** ((mg1 - mg_c) * 0.4)    
    cr4 =  b15 * 10 ** ((mg1 - mg_c) * 0.4)    
    cr5 =  n2 * 10 ** ((mg1 - mg_c) * 0.4)    
    return mg, mg_c, cr1, cr2, cr3, cr4, cr5

countnuv = np.vectorize(countnuv)

# To check the RA, DEC input.
if (DEC.count(':') == RA.count(':') == 2):
     pass
else:
   sys.exit(1)

# To do all the stuff in a specific directory.
working_arena = '.'
os.chdir(working_arena)

# instrument and radius of search in arsec.
field_radius = {'uvit'  : 1200,
                'sxt'   : 1500,
                'czti'  : 1680,
                'laxpc' : 1680}

# the row limit of astroquery is set to 50 by default.
Vizier.ROW_LIMIT = -1

# so that astropy can understand the coordinates. 
radec = RA + " " + DEC
cc = SkyCoord(radec, unit = (u.hourangle, u.deg))

# Searching the GALEX catalogue.
result = Vizier.query_region(cc,
                             radius = field_radius[instrument] * u.arcsec,
                             catalog = 'II/312')

# Stacking the results (there could both AIS and MIS).
result = vstack(result, metadata_conflicts = 'silent')

if len(result) == 0:
    sys.exit()

# FUV 
fd = result['RAJ2000', 'DEJ2000', 'FUV']
fd = fd.filled(99).group_by('FUV')[:5]
ma, ma_c, ta, tb, tc, td = countfuv(fd['FUV'])
restab = Table([ma, ma_c, ta, tb, tc, td],
               names = ('Mag',
                        'Mag_corrected',
                        'caf2',
                        'baf2',
                        'sapphire',
                        'silica'), 
               meta = {'name': 'FUV counts'})

restab['Mag'].format = '4.2f'
restab['Mag_corrected'].format = '4.2f'
restab['caf2'].format = '4.2f'
restab['baf2'].format = '4.2f'
restab['sapphire'].format = '4.2f'
restab['silica'].format = '4.2f'

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
coord = SkyCoord(list(zip(fd['RAJ2000'].data.data, fd['DEJ2000'].data.data)),
                 frame = 'icrs',
                 unit = 'deg')

RAhms_DECdms = coord.to_string('hmsdms', sep = ':')
ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])
xy_tab = Table([ra_hms, dec_dms], names = ('ra_hms', 'dec_dms'))
restab = hstack([xy_tab, restab])

print('\n\n{}\n'.format(restab))

#NUV
nd = result['RAJ2000', 'DEJ2000', 'NUV']
nd = nd.filled(99).group_by('NUV')[:5]
ma, ma_c, ta, tb, tc, td, te = countnuv(nd['NUV'])
restab = Table([ma, ma_c, ta, tb, tc, td, te],
               names = ('Mag',
                        'Mag_corrected', 
                        'silica',
                        'b4',
                        'b13',
                        'b15',
                        'n2'), 
               meta = {'name': 'NUV counts'})

restab['Mag'].format = '4.2f'
restab['Mag_corrected'].format = '4.2f'
restab['silica'].format = '4.2f'
restab['b4'].format = '4.2f'
restab['b13'].format = '4.2f'
restab['b15'].format = '4.2f'
restab['n2'].format = '4.2f'

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
coord = SkyCoord(list(zip(nd['RAJ2000'].data.data, nd['DEJ2000'].data.data)),
                 frame = 'icrs',
                 unit = 'deg')

RAhms_DECdms = coord.to_string('hmsdms', sep = ':')
ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])
xy_tab = Table([ra_hms, dec_dms], names = ('ra_hms', 'dec_dms'))
restab = hstack([xy_tab, restab])

print('\n{}\n\nDone!\n'.format(restab))




