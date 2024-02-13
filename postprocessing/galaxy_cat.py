#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:07:23 2020

@author: flaminia
"""

#MAKE GALAXY CATALOG FROM PARTICLE CATALOG

print ("\n\nUSAGE: python galaxy_cat.py [-options]\n")
print ("OPTIONS: INI_FILE, MODULE, INPUT_PATH, FILTER, FOV_DEG, RES_AS")
print ("NOTE: ...INI_FILE must be created from /lc/planes_list(#cols 1, 6, 7) with following column names: #pl,snap,z .")
print ("      ...MODULE is the module of the input particles catalogs (df,dc or igm).")
print ("      ...INPUT_PATH is ./ by default.")
print ("      ...FILTER_POS is the position (from 0 to N-1) of the filter in the filter file used to weight the flux in the computation of the CM.")
print ("      ...FOV_DEG is the size of the field of view by side in degrees (the same given in image.ini).")
print ("      ...RES_AS is the resolution of the raw image in arcsecs (the same given in image.ini).\n")

import numpy as np 
import pandas as pd
import sys
from astropy.io import ascii
import math


def load_data(file_path):
    with open(file_path, 'r') as file:
        fline = file.readline()
        ncols = len(fline.split())
    data = np.loadtxt(file_path, dtype="double", unpack=True, usecols=range(ncols))
    return data.tolist()


INI_FILE="file.ini"
MODULE="igm"
INPUT_PATH="./"
FILTER=0
FOV=0.
RES=0.


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
        elif '-FILTER' in sys.argv[p]:
            FILTER=int(sys.argv[p+1])
        elif '-FOV_DEG' in sys.argv[p]:
            FOV=float(sys.argv[p+1])
        elif '-RES_AS' in sys.argv[p]:
            RES=float(sys.argv[p+1]) 
        p+=2
    except:
        p=-1


       
truenpix = int(FOV * 3600 / RES)
bufferpix = int(((truenpix + 1) * 20 / 14142)) #int(math.ceil((truenpix + 1) * 20 / 14142)) 
npix = truenpix + bufferpix

print("-INI_FILE", INI_FILE, "-MODULE", MODULE, "-INPUT_PATH", INPUT_PATH, "-FILTER", FILTER,"-FOV_DEG",FOV, "-RES_AS", RES)
print(" COMPUTED N_PIXELS", truenpix, "BUFFER_PIXELS", bufferpix, "\n")



for i,j in zip(snap,pl): 

    out_path = INPUT_PATH+'outcat.'+MODULE+'.%s_%s.txt'%(i,j)
    in_path = INPUT_PATH+'flux.'+MODULE+'.%s_%s.txt'%(i,j)
    data = load_data(in_path) 
    ncols = len(data)
    flux = []

    arrays = [[] for _ in range(ncols)]
    for j in range(len(data[0])):
        for k in range(ncols):
            arrays[k].append(data[k][j])
    #del data
    

    if(MODULE == 'dc' or MODULE == 'igm'):
        idsh = [int(value) for value in arrays[0]]
        x = [int(int(arrays[1][j]*npix+1)-bufferpix/2.+1) for j in range(len(data[0]))]
        y = [int(int(arrays[2][j]*npix+1)-bufferpix/2.+1) for j in range(len(data[0]))]
        zred = arrays[3]
        mass = arrays[4]
        met = arrays[6]
        age = arrays[7]
        ZHI = arrays[8]
        NHI = arrays[9]
        npF = [int(value) for value in arrays[10]]
        npTNG = [int(value) for value in arrays[11]]

        nfilters=ncols-12
        if(ncols>12):
            for k in range(12, nfilters + 12):
                flux.append(arrays[k])

        #dataframe
        fixed_cols = {'id': idsh, 'x': x, 'y': y, 'zr': zred, 'm': mass, 'age': age, 'met': met, 'Npsh': npF, 'NpTNG': npTNG, 'ZHI': ZHI, 'NHI': NHI}
        df = pd.DataFrame(fixed_cols)
        for k in range(nfilters):
            df[f'f{k}'] = flux[k]

        
    elif(MODULE == 'df'):
        idsh = [int(value) for value in arrays[0]]
        x = [int(int(arrays[1][j]*npix+1)-bufferpix/2.+1) for j in range(len(data[0]))]
        y = [int(int(arrays[2][j]*npix+1)-bufferpix/2.+1) for j in range(len(data[0]))]
        zred = arrays[3]
        mass = arrays[4]
        met = arrays[6]
        age = arrays[8]
        npTNG = [int(value) for value in arrays[9]]

        nfilters=ncols-10
        if(ncols>10):
            for k in range(10, nfilters + 10):
                flux.append(arrays[k])

        #dataframe
        fixed_cols = {'id': idsh, 'x': x, 'y': y, 'zr': zred, 'm': mass, 'age': age, 'met': met, 'NpTNG': npTNG}
        df = pd.DataFrame(fixed_cols)
        for k in range(nfilters):
            df[f'f{k}'] = flux[k]
    
    
    #mass-weighted metallicity    
    al=df.groupby(['id']).apply(lambda dfx: (dfx["met"] * dfx["m"]).sum() / dfx["m"].sum()).reset_index()
    al.columns=['id', 'w_met']

    #mass-weighted age    
    al2=df.groupby(['id']).apply(lambda dfx: (dfx["age"] * dfx["m"]).sum() / dfx["m"].sum()).reset_index()
    al2.columns=['id', 'w_age']

    #max age
    h3=df.groupby(['id']).max()['age'].reset_index()
    h3.columns=['id', 'o_age']

    #center weigthed with flux
    xct = df.groupby(['id']).apply(lambda dfx: (dfx["x"] * dfx[f'f{FILTER}']).sum() / dfx[f'f{FILTER}'].sum()).reset_index()
    yct = df.groupby(['id']).apply(lambda dfx: (dfx["y"] * dfx[f'f{FILTER}']).sum() / dfx[f'f{FILTER}'].sum()).reset_index()
    xct.rename(columns = {'id':'idx', 0:'x_c'}, inplace = True)
    yct.rename(columns = {'id':'idy', 0:'y_c'}, inplace = True)
    dfc=pd.concat([xct['idx'],xct['x_c'],yct['y_c']], axis=1)
    dfc.rename(columns = {'idx':'id'}, inplace = True)
    df = df.merge(dfc, on=['id',], how='left')
    
    #radius
    df['dy']=(df['y']-df['y_c'])**2
    df['dx']=(df['x']-df['x_c'])**2    
    df['R'] =(df['dx']+df['dy'])**0.5
    #max radius
    h2=df.groupby(['id']).max()['R'].reset_index()
    h2.columns=['id', 'Rmax']  

    #half-light radius
    df = df.sort_values(['id', 'R'])
    for i in range(nfilters):
        df[f'flux_cumsum{i}'] = df.groupby(['id'])[f'f{i}'].transform(pd.Series.cumsum)
    #half flux  
        column_name = f'half_flux{i}'
        df_half_flux = df.groupby(['id']).apply(lambda x: x[f'f{i}'].sum() / 2).reset_index().rename(columns={0: column_name})
        df = pd.merge(df, df_half_flux, how="left", on=['id'])       
    #discrepancy
        df[f'flux_diff{i}'] = abs(df[f'half_flux{i}'] - df[f'flux_cumsum{i}'])
    #R_hlf
        df_new = df.loc[df.groupby('id')[f'flux_diff{i}'].idxmin()].reset_index()
        df_R = df_new[['id', 'R']]
        df_R.columns = ['id', f'R_hl{i}']
        df = pd.merge(df, df_R, how="left", on=['id'])
       
    ###############
    
    
    #sum columns for galaxy database
    if(MODULE == 'dc' or MODULE == 'igm'):
        sumcols = ['m', 'zr', 'Npsh', 'NpTNG', 'ZHI', 'NHI'] + [f'f{i}' for i in range(nfilters)] + [f'R_hl{i}' for i in range(nfilters)] + [f'half_flux{i}' for i in range(nfilters)]

    elif(MODULE == 'df'):
        sumcols = ['m', 'zr', 'NpTNG'] + [f'f{i}' for i in range(nfilters)] + [f'R_hl{i}' for i in range(nfilters)] + + [f'half_flux{i}' for i in range(nfilters)]
        
    g = df.groupby(['id'])[sumcols].sum().reset_index()
    
    #f in R=3 pixels
    for i in range(nfilters):
        filter_name = f'f{i}'
        df_filter = df.groupby('id').apply(lambda x: x[x['R'] < 3][filter_name].sum()).reset_index()
        df_filter.columns = ['id', f'{filter_name}_3']
        g = pd.merge(g, df_filter, how="left", on=['id'])
              
    #final galaxy database
    g['w_met'] = al['w_met']
    g['w_age']=al2['w_age']
    g['older_age']= h3['o_age']
    g['x_c']=dfc['x_c']
    g['y_c']=dfc['y_c']    
    g['Rmax']=h2['Rmax']
    g['count'] = g['id'].map(df['id'].value_counts())
    g['av_zr'] = g['zr']*1./g['count']
    g['NPTNG'] = (g['NpTNG'] / g['count']).astype(int)
    if(MODULE == 'dc' or MODULE == 'igm'):
        g['NPsh'] = (g['Npsh'] / g['count']).astype(int)
        g['Zhi'] = g['ZHI']*1. / g['count']
	g['Nhi'] = g['NHI']*1. / g['count']
    elif(MODULE=='df'):
        g['NPsh'] = g['count']
        g['Zhi'] = 0.
	g['Nhi'] = 0.

    

    

    for i in range(nfilters):
        g[f'r_hlf{i}'] = g[f'R_hl{i}'] / g['count']
        g[f'hlf{i}'] = g[f'half_flux{i}'] / g['count']

    cols=['id', 'x_c', 'y_c', 'av_zr', 'm', 'w_met', 'w_age', 'older_age','Rmax', 'NPsh', 'NPTNG', 'Zhi', 'Nhi' ]+[f'f{i}' for i in range(nfilters)] + [f'f{i}_3' for i in range(nfilters)]+[f'hlf{i}' for i in range(nfilters)] + [f'r_hlf{i}' for i in range(nfilters)]

    with open(out_path, 'w') as outfile:
        outfile.write("#")
    outfile.close()
    
    #g[cols].to_csv(out_path, sep=' ', index=False, header=lambda x: "#" + ' '.join(x) if x.name in cols else x)
    g[cols].to_csv(out_path, sep=' ', mode='a', index=False)

    print("  ", out_path, "WRITTEN")

print("\n DONE. PLEASE MERGE CATALOGS WITH merge_cat.py FOR THE FINAL GALAXY CATALOG.")

