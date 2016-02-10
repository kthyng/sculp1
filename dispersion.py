'''
Script to run calculate drifters from the SCULP 1 numerical drifters.
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
from matplotlib.mlab import find
import tracpy.calcs

mpl.rcParams.update({'font.size': 16})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'

def plot():
    '''
    Plot relative dispersion for LaCasce data vs. different years I have available.
    '''

    colors = ['.5', '.2', '.4'] # no diff, doturb=1, doturb=2,3
    symbols = [':', '-.', '-'] 
    lacasce = np.loadtxt('lacasce_dispersion_Points.txt')
    lacasce50 = np.loadtxt('lacasce_disperion_50days.txt')

    fig2 = plt.figure(figsize=(14,8))
    ax2 = fig2.add_subplot(111)
    ax2.set_xlabel('Days'); ax2.set_ylabel('Number of drifters in domain');
    fig = plt.figure(figsize=(14,8))
    ax = fig.add_subplot(111)
    ax.set_xlabel('Days'); ax.set_ylabel('Relative Dispersion [km$^2\!$]');
    ax.text(0.58, 0.165, 'Data', color='r', transform=ax.transAxes)
    ax.text(0.58, 0.12, 'Different years', color='0.5', transform=ax.transAxes)
    ax.text(0.58, 0.075, 'Mean of different years', color='0.3', transform=ax.transAxes)
    ax.semilogy(np.array([7,7]), np.array([10**(-1), 10**(5)]), '-', color='lightgrey', linewidth=2)
    # ax2.plot(np.array([7,7]), np.array([1500000, 4000000]), '-', color='lightgrey', linewidth=2)
    # # Skip 2012 because of model output not matching variable set
    # years = np.array([2003,2004,2005,2006,2007,2008,2009,2010,2011,2013])
    years = np.arange(2003, 2014)
    d = netCDF.Dataset('tracks/2003-10-01T00gc.nc')
    t = d.variables['tp']
    days = (t[0]-t[0,0])/(3600.*24) # time vector in days
    d.close()
    D2mean = 0; nnansmean = 0;
    for year in years:
        Dname = 'calcs/' + str(year) + 'D2.npz'
        D = np.load(Dname)
        D2 = D['D2']; nnans = D['nnans']; D.close()
        D2mean += D2*nnans
        nnansmean += nnans
        ax.semilogy(days, D2, color='0.5', linewidth=4, alpha=0.5)#, symbols[doturb], color=colors[doturb], linewidth=4)
        ax2.plot(days,nnans, color='0.5', alpha=0.5, lw=4)
        
    ax.semilogy(days, D2mean/nnansmean, color='0.3', linewidth=6)#, symbols[doturb], color=colors[doturb], linewidth=4)
    ax.semilogy(days[0:301], 7*np.exp(0.55*days[0:301]), 'r--', linewidth=5, alpha=0.5)#.5) # Exponential growth for first 10 days
    ax.semilogy(days[320:], 11*days[320:]**2.2, 'r--', linewidth=5, alpha=0.5)#.5)
    # ax.semilogy(days[302:], 10*days[302:]**2.2, 'r--', linewidth=5, alpha=0.7)#.5)
    ax.semilogy(lacasce[:,0], lacasce[:,1], 'ro', label='LaCasce', markersize=13, alpha=0.9)#, mew=1.5, mec='grey')
    # ax.semilogy(lacasce50[:,0], lacasce50[:,1], 'r*', label='LaCasce', markersize=15, alpha=0.6, markeredgewidth=1.5, mec='grey')
    # labels for analytic functions
    ax.text(0.7, 50, '$e^{0.55t}$', color='r', fontsize=25, fontweight='bold', alpha=0.7)
    ax.text(20, 17000, '$t^{2.2}$', color='r', fontsize=25, fontweight='bold', alpha=0.7)
    ax.set_ylim(.5,10**5)
    ax.set_xlim(0, 25)
    ax2.plot(days, nnansmean, 'k')
    fig.savefig('figures/relative_dispersion_comp.pdf', bbox_inches='tight')
    fig2.savefig('figures/relative_dispersion_comp-nnans.pdf', bbox_inches='tight')


def run():
    '''
    Run code to save dispersion calculations.
    '''

    # Skip 2012 because of model output not matching variable set
    years = np.array([2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013])

    # Location of TXLA model output
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # grid = tracpy.inout.readgrid(loc)
    grid_filename = '/atch/raid1/zhangxq/Projects/txla_nesting6/txla_grd_v4_new.nc'
    vert_filename='/atch/raid1/zhangxq/Projects/txla_nesting6/ocean_his_0001.nc'
    grid = tracpy.inout.readgrid(grid_filename, vert_filename=vert_filename, usebasemap=True)

    for year in years:

        Files = glob.glob('tracks/' + str(year) + '-1[0,1,2]-*gc.nc')
        Files.extend(glob.glob('tracks/' + str(year+1) + '-0[1-7]-*gc.nc'))

        D2nameyear = os.path.join('calcs', str(year) + 'D2.npz')
        print D2nameyear
        if os.path.exists(D2nameyear): continue;
        D2 = []; nnans = [];
        for File in Files: # loop through all the runs of that type
            D2name = 'calcs/' + File[:-5].split('/')[-1] + 'D2.npz'
            if not os.path.exists(D2name):
                print File
                d = netCDF.Dataset(File)
                xg = d.variables['xg'][:]; yg = d.variables['yg'][:]
                # eliminate entries equal to -1
                ind = xg==-1
                xg[ind] = np.nan; yg[ind] = np.nan
                xp, yp, _ = tracpy.tools.interpolate2d(xg, yg, grid, 'm_ij2ll')
                d.close()
                # D^2 is already averaged
                D2_temp, nnans_temp, pairs = tracpy.calcs.rel_dispersion(xp, yp, r=[0, 1.05], squared=True)
                np.savez(D2name, D2=D2_temp, nnans=nnans_temp) 
            else:
                d = np.load(D2name)
                D2_temp = d['D2']; nnans_temp = d['nnans'];
                d.close()
                # pdb.set_trace()
            D2.append(D2_temp*nnans_temp) # un-average
            nnans.append(nnans_temp)

        # After I have run through all the times for this year, do average and save
        D2 = np.nansum(np.asarray(D2), axis=0) # sum the individual sums (already squared)
        nnans = np.nansum(np.asarray(nnans), axis=0) # sum non-nans for averages
        D2 = D2.squeeze()/nnans
        np.savez(D2nameyear, D2=D2, nnans=nnans)


if __name__ == "__main__":
    run()    
