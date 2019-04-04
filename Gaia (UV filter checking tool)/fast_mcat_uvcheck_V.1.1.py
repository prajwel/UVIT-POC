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

   The author would like to acknowledge inputs from Dr. Koshy George
   which greatly helped the developement of this script.
   
'''


import os
import sys
import string
import urllib
import numpy as np

from astropy.io import fits
from requests import Session
from bs4 import BeautifulSoup
from astropy import units as u
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord


# To get the user input. 
instrument = str(sys.argv[1])
RA = str(sys.argv[2])
DEC = str(sys.argv[3])

#instrument = 'uvit'
#RA = "7:36:51.396"
#DEC = "65:36:9.170"

# To check the RA, DEC user input.
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

    if mg <= 15.0:
        mg_c = 2.634 + (26.316 * mg - 245.329) ** 0.5
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

# Function to find seperation in celestial coordinates.
cc = SkyCoord(RA, DEC, unit = (u.hourangle, u.deg))
def cel_separation(a, b):
    coo = SkyCoord(a, b, frame = 'icrs', unit = 'deg')
    return coo.separation(cc)


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
           '_ctl10:txtRadius': '0.4',
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
ss = soup.find(id = '_ctl10_TileGrid_imgLink_0').get('href')
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
    sys.exit(1)

# Reading coordinates from catalogue.
hdu = fits.open(catalogue)
alpha = hdu[1].data['alpha_j2000_merged']
delta = hdu[1].data['delta_j2000_merged']


# NUV 
nuv_mag = hdu[1].data['nuv_mag']
refined_set = [(al, de, nm) for al, de, nm 
                            in zip(alpha, delta, nuv_mag)
                            if int(nm) != -999 and nm <= 22.]

nalpha, ndelta, nuv_mag = zip(*refined_set)

confined_set = [(nm, al, de) for al, de, nm 
                             in zip(nalpha, ndelta, nuv_mag) 
                             if cel_separation(al, de) <= field_radius[instrument] * u.arcsec]

nd = np.array(sorted(confined_set))[0:5]
ma, ma_c, ta, tb, tc, td, te = countnuv(nd[:,0])
nuv_res = Table([ma, ma_c, ta, tb, tc, td, te],
               names = ('Mag',
                        'Mag_corrected', 
                        'silica',
                        'b4',
                        'b13',
                        'b15',
                        'n2'), 
               meta = {'name': 'NUV counts'})

nuv_res['Mag'].format = '4.2f'
nuv_res['Mag_corrected'].format = '4.2f'
nuv_res['silica'].format = '4.2f'
nuv_res['b4'].format = '4.2f'
nuv_res['b13'].format = '4.2f'
nuv_res['b15'].format = '4.2f'
nuv_res['n2'].format = '4.2f'

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
coord = SkyCoord(zip(nd[:,1], nd[:,2]),
                 frame = 'icrs',
                 unit = 'deg')

RAhms_DECdms = coord.to_string('hmsdms', sep = ':')
ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])
xy_tab = Table([ra_hms, dec_dms], names = ('ra_hms', 'dec_dms'))
nuv_res = hstack([xy_tab, nuv_res])

print('\n\n### NUV\n\n{}\n'.format(nuv_res))

# FUV 
fuv_mag = hdu[1].data['fuv_mag']

fuv_absent = 'no'
if len(np.unique(fuv_mag)) == 1:  # when FUV data is absent.
    fd = nd
    print('\nFUV observations seem to be absent! Using M_fuv = M_nuv - 1.65.')
    fuv_absent = 'yes'
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
    ma, ma_c, ta, tb, tc, td = countfuv(fd[:,0])
else:
    ma, ma_c, ta, tb, tc, td = countfuv_abs(fd[:,0])

fuv_res = Table([ma, ma_c, ta, tb, tc, td],
               names = ('Mag',
                        'Mag_corrected',
                        'caf2',
                        'baf2',
                        'sapphire',
                        'silica'), 
               meta = {'name': 'FUV counts'})

fuv_res['Mag'].format = '4.2f'
fuv_res['Mag_corrected'].format = '4.2f'
fuv_res['caf2'].format = '4.2f'
fuv_res['baf2'].format = '4.2f'
fuv_res['sapphire'].format = '4.2f'
fuv_res['silica'].format = '4.2f'

# To convert ra_deg and dec_deg to ra_hms and dec_dms.
coord = SkyCoord(zip(fd[:,1], fd[:,2]),
                 frame = 'icrs',
                 unit = 'deg')

RAhms_DECdms = coord.to_string('hmsdms', sep = ':')
ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])
xy_tab = Table([ra_hms, dec_dms], names = ('ra_hms', 'dec_dms'))
fuv_res = hstack([xy_tab, fuv_res])

print('\n### FUV \n\n{}\n\nDone!\n'.format(fuv_res))



