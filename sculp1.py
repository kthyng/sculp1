'''
Script to run drifters within the grid of SCULP 1 drifters and at the same timing
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
from datetime import datetime, timedelta
import glob
from tracpy.tracpy_class import Tracpy


grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
# currents_filename = list(np.sort(glob.glob('/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_????.nc')))

# grid_filename = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
# can't aggregate years between 2012 and before with 2013 and 2014 bc they have different variables
# years = np.arange(2011,2013)
# currents_filename = []
# for year in years:
#     currents_filename.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_????.nc')))

# years = np.arange(2013,2015)
# currents_filename = []
# for year in years:
#     currents_filename.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_*.nc')))

years = np.arange(2003,2015)
currents_filename = []
for year in years:
    currents_filename.extend(np.sort(glob.glob('/home/kthyng/shelf/' + str(year) + '/ocean_his_*.nc')))

# grid = tracpy.inout.readgrid(grid_filename, usebasemap=True)
grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)


def init(name):
    '''
    Initialization for the simulation.
    '''

    # loc = 'http://barataria.tamu.edu:6060/thredds/dodsC/NcML/txla_nesting6.nc'

    time_units = 'seconds since 1970-01-01'

    # horizontal_diffusivity project showed that relative dispersion did not
    # change between nsteps=25 and 50, but does between nsteps=5 and 25, and
    # interim numbers have not been tested yet.
    nsteps = 25 # in-between tracks: 12 # old tracks: 25 

    # Number of steps to divide model output for outputting drifter location
    N = 5

    # Number of days
    ndays = 50

    # This is a forward-moving simulation
    ff = 1 

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0. # old tracks: 5.
    av = 0. # m^2/s

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines.
    dostream = 0

    # Initialize Tracpy class
    tp = Tracpy(currents_filename, grid_filename=grid_filename, name=name, tseas=tseas, ndays=ndays, nsteps=nsteps, dostream=dostream, savell=False, doperiodic=0, 
                N=N, ff=ff, ah=ah, av=av, doturb=doturb, do3d=do3d, z0=z0, zpar=zpar, 
                time_units=time_units, usebasemap=True, grid=grid, vert_filename=vert_filename)

    # tp._readgrid()

    if os.path.exists('calcs/seeds.npz'):
        seeds = np.load('calcs/seeds.npz')
        lon0 = seeds['lon0']; lat0 = seeds['lat0']
        seeds.close()
    else:
        llcrnrlon = -93.8; urcrnrlon = -92.2; llcrnrlat = 28; urcrnrlat = 29.2; # New
        xcrnrs, ycrnrs = grid['basemap']([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
        X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], 700), 
                            np.arange(ycrnrs[0], ycrnrs[1], 700))
        lon0, lat0 = grid['basemap'](X, Y, inverse=True)

        # Eliminate points that are outside domain or in masked areas
        lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grid)
        
        # save starting locations for future use
        np.savez('calcs/seeds.npz', lon0=lon0, lat0=lat0)

    # # equal weightings for drifters for transport.
    # T0 = np.ones(lon0.size, order='F')

    # U = np.ma.zeros(tp.grid['xu'].shape, order='F')
    # V = np.ma.zeros(tp.grid['xv'].shape, order='F')

    # pdb.set_trace()
       
    return tp, lon0, lat0


def run():

    year = 2013

    # Weekly Oct, Nov, Dec; biweekly Jan, Feb, Mar; monthly Apr, May, Jun, Jul
    startdates = np.array([datetime(year, 10, 1, 0, 1), datetime(year, 10, 8, 0, 1),
                            datetime(year, 10, 15, 0, 1), datetime(year, 10, 22, 0, 1),
                            datetime(year, 11, 1, 0, 1), datetime(year, 11, 8, 0, 1),
                            datetime(year, 11, 15, 0, 1), datetime(year, 11, 22, 0, 1),
                            datetime(year, 12, 1, 0, 1), datetime(year, 12, 8, 0, 1),
                            datetime(year, 12, 15, 0, 1), datetime(year, 12, 22, 0, 1),
                            datetime(year+1, 1, 1, 0, 1), datetime(year+1, 1, 15, 0, 1),
                            datetime(year+1, 2, 1, 0, 1), datetime(year+1, 2, 15, 0, 1),
                            datetime(year+1, 3, 1, 0, 1), datetime(year+1, 3, 15, 0, 1),
                            datetime(year+1, 4, 1, 0, 1), datetime(year+1, 5, 1, 0, 1),
                            datetime(year+1, 6, 1, 0, 1), datetime(year+1, 7, 1, 0, 1)])

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')
        
    # loop through state dates
    for startdate in startdates:

        date = startdate

        name = date.isoformat()[0:13]

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc') and \
            not os.path.exists('tracks/' + name + 'gc.nc'):

            # Read in simulation initialization
            tp, lon0, lat0 = init(name)

            # Run tracpy
            # Save directly to grid coordinates
            lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)


        # Increment by 24 hours for next loop, to move through more quickly
        date = date + timedelta(hours=24)


if __name__ == "__main__":
    run()    
