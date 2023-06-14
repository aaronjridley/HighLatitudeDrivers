#!/usr/bin/env python

import numpy as np
from amie_routines import *
import sys
from netCDF4 import Dataset
import datetime as dt

# ----------------------------------------------------------------------
# read netcdf file
# ----------------------------------------------------------------------

def read_gang_file(filename):

    data = {}

    data["Vars"] = ['Potential (kV)', \
                    'Electron Energy Flux (ergs/cm2/s)', \
                    'Electron Mean Energy (keV)']

    
    ncfile = Dataset(filename, 'r')
    data['lats'] = np.array(ncfile.variables['lat'])
    data['mlts'] = np.array(ncfile.variables['lon'])/15.0
    data[data["Vars"][0]] = np.array(ncfile.variables['pot'])
    data[data["Vars"][1]] = np.array(ncfile.variables['efx'])
    data[data["Vars"][2]] = np.array(ncfile.variables['ekv'])

    nTimes, nLats, nMlts = data[data["Vars"][0]].shape
    data["nLats"] = nLats
    data["nMlts"] = nMlts
    data["nTimes"] = nTimes

    year = np.array(ncfile.variables['year'])
    month = np.array(ncfile.variables['month'])
    day = np.array(ncfile.variables['day'])
    ut = np.array(ncfile.variables['ut'])

    hour = ut.astype(int)
    minute = ((ut - hour)*60.0).astype(int)
    second = (ut*3600.0 % 60.0).astype(int)
    time = []

    imf = [0.0, 0.0, 0.0, 0.0]
    ae = [0.0, 0.0, 0.0, 0.0]
    dst = [0.0, 0.0]
    hp = [0.0, 0.0]
    hpi = np.array(ncfile.variables['hpi'])
    pcp = np.array(ncfile.variables['pcp'])

    data['imf'] = []
    data['ae'] = []
    data['dst'] = []
    data['hp'] = []
    data['cpcp'] = []
    
    for i in range(len(year)):
        s = second[i]
        t = dt.datetime(year[i],
                        month[i],
                        day[i],
                        hour[i],
                        minute[i],
                        second[i])
        if (s == 59):
            t = t + dt.timedelta(seconds = 1)
        time.append(t)

        data['imf'].append(imf)
        data['ae'].append(ae)
        data['dst'].append(dst)
        hp[0] = hpi[i]
        data['hp'].append(hp)
        data['cpcp'].append(pcp[i])

    
    ncfile.close()

    data['times'] = time
    data["nVars"] = len(data["Vars"])
    data["version"] = 0.0
    
    return data

north_file = 'mar15_20_2015_amp_nh.nc'
south_file = 'mar15_20_2015_amp_sh.nc'

nData = read_gang_file(north_file)
sData = read_gang_file(south_file)

sTime = nData['times'][0].strftime('%Y%m%d')
eTime = nData['times'][-1].strftime('%Y%m%d')

print(sTime, eTime)

base = 'amie_' + sTime + '_to_' + eTime

amie_write_binary(base + '_north.bin', nData)
amie_write_binary(base + '_south.bin', sData)
