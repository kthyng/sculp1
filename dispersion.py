'''
Script to run drifters backward from Galveston Bay to examine the Bay's 
connectivity with the shelf region.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob
from matplotlib.mlab import find


def run():
    '''
    Run code to save dispersion calculations.
    '''

    # Location of TXLA model output
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # grid = tracpy.inout.readgrid(loc)
    grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
    grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

    Files = glob.glob('tracks/*gc.nc')

    # Dnameoverall = os.path.join(test, 'D2overall.npz')
    # D2 = []; nnans = [];
    for File in Files: # loop through all the runs of that type
        D2name = 'calcs/' + File[:-5].split('/')[-1] + 'D2.npz'
        if not os.path.exists(D2name):
            d = netCDF.Dataset(File)
            xg = d.variables['xg'][:]; yg = d.variables['yg'][:]
            # eliminate entries equal to -1
            ind = xg==-1
            xg[ind] = np.nan; yg[ind] = np.nan
            xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')
            d.close()
            # D^2 is already averaged
            D2_temp, nnans_temp, pairs = tracpy.calcs.rel_dispersion(xp, yp, r=1.05, squared=True)
            np.savez(D2name, D2=D2_temp, nnans=nnans_temp) 
        else:
            d = np.load(D2name)
            D2_temp = d['D2']; nnans_temp = d['nnans'];
            d.close()
            # pdb.set_trace()
        # D2.append(D2_temp)
        # nnans.append(nnans_temp)

    # # After I have run through all the times for this type of run, do average and save
    # D2 = np.nansum(np.asarray(D2), axis=0) # sum the individual sums (already squared)
    # nnans = np.nansum(np.asarray(nnans), axis=0) # sum non-nans for averages

    # D2 = D2.squeeze()/nnans
    # np.savez(Dnameoverall, D2=D2, t=t)


if __name__ == "__main__":
    run()    
