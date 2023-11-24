#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 11:13:43 2020

@author: flaminia
"""
#merge outcat (galaxy cat) to create final light-cone cat

print ("\n\nUSAGE: python merge_cat.py [-options]\n")
print ("OPTIONS: INI_FILE, MODULE, INPUT_PATH, ADD_TO_FILENAME ")
print ("NOTE: ...INI_FILE must be created from /lc/planes_list(#cols 1, 6, 7) with following column names: #pl,snap,z")
print ("      ...MODULE is df,dc or igm.")
print ("      ...INPUT_PATH is ./ by default.")
print ("      ...ADD_TO_FILENAME is something you want to add to the name of the final catalog.\n")


import numpy as np 
import pandas as pd
import sys
from astropy.io import ascii

INI_FILE="file.ini"
INPUT_PATH="./"
MODULE="igm"
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
        elif '-MODULE' in sys.argv[p]:
            MODULE=str(sys.argv[p+1])
        elif '-INPUT_PATH' in sys.argv[p]:
            INPUT_PATH=str(sys.argv[p+1])
        elif '-ADD_TO_FILENAME' in sys.argv[p]:
            ADD="."+str(sys.argv[p+1])
        p+=2
    except:
        p=-1

print("-INI_FILE", INI_FILE,"-MODULE", MODULE, "-INPUT_PATH", INPUT_PATH, "-ADD_TO_FILENAME",ADD)



filenames = [INPUT_PATH+'outcat.'+MODULE+'.'+str(n)+'_'+str(m)+'.txt' for n,m in zip(snap,pl) ]

write_header = True
final_file = INPUT_PATH+'final_cat.'+MODULE+ADD+'.dat'
with open(final_file, 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            if write_header:
                header = next(infile)
                outfile.write(header)
                write_header = False
            else:
                next(infile)
            outfile.write(infile.read())

print("\n FINAL GALAXY CATALOG WRITTEN", final_file)
