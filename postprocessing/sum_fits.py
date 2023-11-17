#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:29:03 2020

@author: flaminia
"""

#CREATE FINAL pre-POST PROCESSING FITS FILE

print ("\n\nUSAGE: python sum_fits.py IMAGE PSF [-options]")
print ("OPTIONS: INI_FILE, FILTER, MODULE, INPUT_PATH, ADD_TO_FILENAME ")
print ("NOTE: ...\n")


import glob
import numpy as np
from astropy.io import fits
from astropy.stats import mad_std
import sys
from astropy.io import ascii

INI_FILE="file.ini"
FILTER="no_filter"
MODULE="igm"
INPUT_PATH="./"
ADD=""

p=1
while p>0:
    try:
        if '-INI_FILE' in sys.argv[p]:
            INI_FILE=str(sys.argv[p+1])
            inp=ascii.read(sys.argv[p+1])
            snap=inp['snap']
            pl=inp['pl']-1
            z_pl=inp['z']
        elif '-FILTER' in sys.argv[p]:
            FILTER=str(sys.argv[p+1])       
        elif '-MODULE' in sys.argv[p]:
            MODULE=str(sys.argv[p+1])
        elif '-INPUT_PATH' in sys.argv[p]:
            INPUT_PATH=str(sys.argv[p+1])
        elif '-ADD_TO_FILENAME' in sys.argv[p]:
            ADD="."+str(sys.argv[p+1])
        p+=2
    except:
        p=-1


print("-INI_FILE", INI_FILE,"-FILTER", FILTER, "-MODULE", MODULE, "-INPUT_PATH", INPUT_PATH, "-ADD_TO_FILENAME",ADD)
image_list = [ INPUT_PATH+FILTER+'.'+MODULE+'.'+str(n)+'_'+str(m)+'.fits' for n,m in zip(snap,pl) ]
header_list=[]


image_concat = []
for image in image_list:
    image_concat.append(fits.getdata(image))
    header = fits.getheader(image)
    header_list.append(header)

z_first = header_list[0]["REDSHIFT"]
z_last = header_list[-1]["REDSHIFT"]
truenpix = header_list[0]["N_PIXEL_BY_SIDE"]
bufferpix = header_list[0]["BUFFER_PIXELS"]
pixu = header_list[0]["PIXELUNIT"]
ntotxy4 = sum(header["N_STELLAR_PARTICLES"] for header in header_list)
fov = header_list[0]["PHYSICALSIZE_BY_SIDE"]
h0 = header_list[0]["HUBBLE"]
dl_first = float(header_list[0]["Dl_UP"])*h0
dl_last = float(header_list[-1]["Dl_UP"])*h0

print(z_first,z_last,truenpix,bufferpix,pixu,ntotxy4,fov,h0,dl_first,dl_last)

final_image = np.zeros(shape=image_concat[0].shape)

for image in image_concat:
    final_image += image



outfile = 'stacked.'+FILTER+"."+MODULE+ADD+'.fits'

header = fits.Header()
header["HIERARCH N_PIXEL_BY_SIDE"] = truenpix
header["HIERARCH BUFFER_PIXELS"] = bufferpix
header["HIERARCH PHYSICALSIZE_BY_SIDE_DEG"] = fov
header["HIERARCH PIXEL_UNIT"] = pixu
header["HIERARCH N_STELLAR_PARTICLES"] = ntotxy4
header["HIERARCH FILTER"] = FILTER
header["HIERARCH LOWER_REDSHIFT"] = z_first
header["HIERARCH UPPER_REDSHIFT"] = z_last
header["HIERARCH LOWER_Dl_Mpc"] = dl_first
header["HIERARCH UPPER_Dl_Mpc"] = dl_last
header["h0"] = h0



hdu = fits.PrimaryHDU(data=final_image, header=header)
hdu.writeto(outfile, overwrite=True)

print("\n:::::FINAL IMAGE", outfile, ":::::")
