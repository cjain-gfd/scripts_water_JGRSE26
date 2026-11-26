# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# updated on 26 Nov 2026 with legends
# to plot rms surface velocities, mobility and inferred tectonic regimes
# for models presented in the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 2. Contrastive effects of water
# submitted to JGR: Solid Earth

# This script uses StagPy library to post-process StagYY data.
# Please cite as:
# Adrien Morison, Martina Ulvrova, Stephane Labrosse, B4rsh, theofatou, & tfrass49. (2024).
# StagPython/StagPy: v0.20.1 (v0.20.1). Zenodo. https://doi.org/10.5281/zenodo.13684760

# Other python packages required: os, math, numpy, h5py, pandas, matplotlib

# Once required packages are installed, execute the script from terminal with command "python plot_vrms_mobility_tectonics.py"

import os
import math
import numpy as np
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (1/3/5/7/11) as batch
batch=11
### end of input ###


if batch == 1:
   model_folder = 'Table4_models/'
   model_list = ["03_AT19r1/"]
elif batch == 3:
   model_folder = 'Table5_models/'
   model_list = ["07_QW5r1/", "13_QW10hr1/"]   
elif batch == 5:
   model_folder = 'Table6_models/'   
   model_list = ["06_Qaw20r3/", "07_Qaw20sat20r1/"]
elif batch == 7:
   model_folder = 'Table6_models/'   
   model_list = ["19_Qsat20H2Oetar3/", "21_Qsat25H2Oetar2/"]
elif batch == 11:
   model_folder = 'Table6_models/'   
   model_list = ["12_Qaw20sat25r3/", "16_Qf02aw20sat40r2/"]
else:
   print('Wrong figure')

here = os.getcwd()

if os.path.isdir(here+'/../figures/'):
   print('Directory for figures exists')
else:
   print('Making directory for figures')
   os.mkdir(here+'/../figures')

total_runs = len(model_list)
print('# models being processed: '+str(total_runs))

# customising matplotlib and figure
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Helvetica"
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.titleweight'] = "bold"
mpl.rcParams['axes.labelweight'] = "normal"
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [4.5,1] # (1 panel, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams['ytick.minor.visible'] = False
mpl.rcParams["pdf.use14corefonts"] = True

box_ratio=8
evolved_time = 4500.
time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]
max_mobility = 3.0
mobility_ticks=[0.0,1.0,2.0,3.0]
min_vrms = 1e-3
max_vrms = 1e3
vrms_ticks=[min_vrms,1,max_vrms]

ML=['#024b79']
PSL=['#43b7c2']
RO=['#ffad48']
EL=['#c25b60']
SL=['#822b42']

velocity_rms_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_surface_vrms_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
steps_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_plateness_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_time_avg = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
tmean_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_step_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
TotalTTGmass_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

seconds_in_year = 31556926
seconds_in_Myear = seconds_in_year*1e6
earth_radius = 6371.e3
core_radius = 3481.e3

# empty lists
timesteps = []
time_TTGmass   = []
TTG1   = []
TTG2   = []
TTG3   = []
frames = []
steps  = []
time = []

# define size of moving window for rolling average
window_size_mobility  = 10
window_size_plateness = 10


run = int(0)
for path in model_list:

   volume_mantle = (4 * math.pi*(earth_radius**3-core_radius**3))/3
   # quadrant analytical volume is 1/16 of full annulus, which in itself is pi/ny times actual volume

   if batch == 1:
      nytot = 2048
      numcores_files = 32
      analytical_volume_full = (4 * math.pi*(earth_radius**3-core_radius**3))/3 * (math.pi)/(nytot)
      factor=volume_mantle/analytical_volume_full
      if run == 0:
         print('Scaling factor for full annulus: '+str(factor))      
   else:
      nytot = 512
      numcores_files = 8
      analytical_volume_quad = (4 * math.pi*(earth_radius**3-core_radius**3))/3 * (math.pi)/(nytot*16)
      factor=volume_mantle/analytical_volume_quad
      if run == 0:
         print('Scaling factor for quadrant annulus: '+str(factor))

   print('Run # '+str(run+1))
   print('Model  '+path)

   nz = 128

   # from stagpy import stagyydata
   sdat = stagyydata.StagyyData(here+'/../'+str(model_folder)+path)

   # extracting time from time_botT.h5
   # this is needed when plotting the time evolution of a field
   time_list = []
   f = h5py.File(here+'/../'+str(model_folder)+path+'+hdf5/time_botT.h5', 'r')

   ndatasets = len(f.keys())
   for key in f.keys():
      data = f[key]
      time_list.append(float(data[0]))
   f.close()

   time=np.array(time_list)
   time_Myr=time/seconds_in_Myear
   print('# HDF datasets: '+str(ndatasets))
   time_Myr_run[run]=time_Myr
   time_list.clear()


   # reading in the timesteps when snapshots are saved
   snaps = []
   file = (here+'/../'+str(model_folder)+path+'snaps.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      snaps.append(int(p[0]))
   nframes = int(len(snaps))
   print('# frames (rprof.dat): ' + str(nframes))
   print('Evolution time (Myr):    '+str('{:.2f}'.format(sdat.tseries.time[len(sdat.tseries.time)-1]*1e-6/seconds_in_year)))
   # timesteps for initial and final states
   step1=snaps[1]
   step2=snaps[nframes-1]


   # Using pandas to calculate moving average (with frames) of mobility and plateness
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_plates_analyse.dat')
   df = pd.read_table(file,sep='  ',engine='python', header=None, skiprows=1,usecols=range(6))                                                       
   df = df.set_axis(['istep','time', 'Rtp_vol', 'Rtp_surf', 'mobility', 'f80'], axis=1)

   df['mobility_10'] = df.mobility.rolling(window_size_mobility, min_periods=1).mean()
   moving_average_mobility=np.array(df['mobility_10'])
   moving_average_mobility_run[run] = moving_average_mobility
   df['plateness'] = 1.-(df['f80']/0.6)
   df['plateness_10'] = df.plateness.rolling(window_size_plateness, min_periods=1).mean()
   moving_average_plateness=np.array(df['plateness_10'])
   moving_average_plateness_run[run] = moving_average_plateness


   # loading mobility and plateness from _plates_analyse.dat
   mobility_list = []
   f80_list = []
   plateness_list = []
   time_list = []
   steps_mobility_list = []
   file = (here+'/../'+str(model_folder)+path+'+op/earthDIM_plates_analyse.dat')
   f = open(file,'r')
   lines=f.readlines()[1:]
   f.close
   nlines=int(len(lines))

   for i in np.arange(0,nlines):
      p = lines[i].split()
      steps_mobility_list.append(int(p[0]))
      time_list.append(float(p[1]))
      mobility_list.append(float(p[4]))
      f80_list.append(float(p[5]))
      plateness_list.append(1.-(float(p[5])/0.6))

   # calculating time-averaged mobility
   sum=0.
   for i in np.arange(0,nlines-1):
      sum = sum + (time_list[i+1]-time_list[i]) * mobility_list[i]

   mobility_time_avg[run] = sum/time_list[nlines-1]
   print('Time-averaged mobility:  '+str('{:.2f}'.format(mobility_time_avg[run])))

   # calculating time-averaged plateness
   sum=0.
   for i in np.arange(0,nlines-1):
      sum = sum + (time_list[i+1]-time_list[i]) * plateness_list[i]

   plateness_time_avg[run] = sum/time_list[nlines-1]
   # print('Time-averaged plateness:      '+str('{:.2f}'.format(plateness_time_avg[run])))

   mobility_run[run] = np.array(mobility_list)
   plateness_run[run] = np.array(plateness_list)
   steps_mobility_run[run] = np.array(steps_mobility_list)
   time_Myr_mobility_run[run]=np.array(time_list)/seconds_in_Myear

   mobility_list.clear()
   f80_list.clear()
   plateness_list.clear()
   steps_mobility_list.clear()
   time_list.clear()


   # computing total mass of TTG produced
   current_model=path[0:len(path)-1]
   if current_model == '12_Qaw20sat25r3' or current_model == '16_Qf02aw20sat40r2':
      file = (here+'/../'+str(model_folder)+path+'TTGmass.dat')
      f = open(file,'r')
      lines = f.readlines()
      f.close
      for line in lines:
         p = line.split()
         timesteps.append(int(p[0]))
         time_TTGmass.append(float(p[1]))
         TTG1.append(float(p[2]))
         TTG2.append(float(p[3]))
         TTG3.append(float(p[4]))

      finaltime=time_TTGmass[len(time_TTGmass)-1]*1e-6/seconds_in_year
      time_Myr_step = np.array(time_TTGmass)*1e-6/seconds_in_year

      # importing arrays
      TTG1v = np.array(TTG1)
      TTG2v = np.array(TTG2)
      TTG3v = np.array(TTG3)

      # declaring empty array
      TTGmass = np.array([0])

      # first sum med-P and high-P TTG
      # then add that to low-P TTG
      TTGmass = np.add(TTG1v,np.add(TTG2v,TTG3v))

      # getting cumulative mass
      TTGmass_cumulative = np.cumsum(TTGmass)
      lowP_cumulative = np.cumsum(TTG1v)
      medP_cumulative = np.cumsum(TTG2v)
      highP_cumulative = np.cumsum(TTG3v)

      print('TTG final mass (scaled):     '+str('{:.2e}'.format(TTGmass_cumulative[len(time_TTGmass)-1]*factor)))

      # clean up for the next run
      timesteps.clear()
      time_TTGmass.clear()
      TTG1.clear()
      TTG2.clear()
      TTG3.clear()

      time_Myr_step_run[run] = time_Myr_step
      TotalTTGmass_run[run] = TTGmass_cumulative * factor


   # loading velocities (in the top cell, nz-1) from radial profiles
   velocity_rms_list = [0] # this is vsurf_rms
   for i in np.arange(1,nframes):
      velocity_rms_list.append(float(sdat.steps[snaps[i]].rprofs['vrms'].values[nz-1]*100*seconds_in_year))

   velocity_rms_run[run] = np.array(velocity_rms_list)
   velocity_rms_list.clear()


   # plotting
   fig, ax = plt.subplots(1, 1, squeeze=False)
   ax[0][0].plot(time_Myr_run[run],velocity_rms_run[run], linestyle='-', linewidth=1, alpha=1, color = '#000000')
   

   current_model=path[0:len(path)-1]
   loc1='upper right'
   loc2='lower right'

   if current_model == '03_AT19r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,270.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],270.,520.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],520.,930.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],930.,1100.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1100.,1550.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1550.,1730.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1730.,4500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(231., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(987., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3403., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower center'
      loc2='lower right'
   elif current_model == '07_QW5r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,190.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],190.,440.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],440.,700.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],700.,980.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],980.,1120.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1120.,1330.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1330.,1600.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1600.,1800.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1800.,3450.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3450.,3800.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3800.,4500.,   linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(501., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1954., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3360., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower center'
      loc2='upper right'
   elif current_model == '13_QW10hr1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,350.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],350.,420.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],420.,560.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],560.,750.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],750.,980.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],980.,1100.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1100.,1300.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1300.,1400.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1400.,1760.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1760.,1880.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1880.,2850.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2850.,3350.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3350.,3500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3500.,3850.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3850.,4020.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4020.,4260.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4260.,4500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(437., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1993., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3118., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower center'
      loc2='upper right'
   elif current_model == '06_Qaw20r3':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,200.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],200.,420.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],420.,500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],500.,1500.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1500.,1600.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1600.,2000.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2000.,2300.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2300.,2550.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2550.,2680.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2680.,2980.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2980.,3100.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3100.,3300.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3300.,3400.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3400.,3900.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3900.,4250.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4250.,4500.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].axvline(76., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1937., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3121., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower left'
      loc2='upper right'
   elif current_model == '07_Qaw20sat20r1':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,400.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],400.,500.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],500.,600.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],600.,720.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],720.,800.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],800.,880.,   linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],880.,1350.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1350.,1500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1500.,2850.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],2850.,3400.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3400.,3730.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3730.,3850.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3850.,4300.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4300.,4500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].axvline(379., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1839., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3150., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper center'
      loc2='upper right'
   elif current_model == '19_Qsat20H2Oetar3':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,270.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],270.,400.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],400.,650.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],650.,830.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],830.,1000.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1000.,1100.,   linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1100.,3400.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3400.,4000.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4000.,4350.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4350.,4500.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].axvline(253., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(2111., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3033., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper right'
      loc2='lower right'
   elif current_model == '21_Qsat25H2Oetar2':
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],500.,650.,  linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],650.,790.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],790.,950.,  linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],950.,1050.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1050.,1150.,  linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1150.,1230.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1230.,1300.,  linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1300.,1450.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1450.,1650.,  linewidth=0,color = PSL,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1650.,3550.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],3550.,4100.,  linewidth=0,color = RO,alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],4100.,4500.,  linewidth=0,color = ML,alpha=0.8,edgecolor='none')
      ax[0][0].axvline(281., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(2080., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(3169., color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='lower center'
      loc2='upper right'

   ax[0][0].set_xlim(0., evolved_time)
   ax[0][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)      
   ax[0][0].set_ylim(min_vrms,max_vrms)
   ax[0][0].set_yscale('log')
   ax[0][0].set_yticks(vrms_ticks)
   ax[0][0].set_ylabel('$v_{rms}$ (cm/yr)',labelpad=4)
   ax[0][0].set_xlabel('Time (Myr)')
   ax[0][0].set_box_aspect(1/box_ratio)
   ax[0][0].legend(ncol = 1, labels=['v$_{rms}$'], loc=loc1, \
   labelspacing = 0.5, edgecolor = 'none', facecolor = 'none', fancybox = False)   

   ax2 = ax[0][0].twinx()
   ax2.plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=0.7, color = '#ffffff')
   ax2.set_xlim(0., evolved_time)
   ax2.set_xscale('linear')
   ax2.set_xticks(time_ticks)      
   ax2.set_yscale('linear')
   ax2.set_ylim(0.,max_mobility)
   ax2.set_yticks(mobility_ticks)
   ax2.set_ylabel('Mobility',labelpad=4)
   ax2.set_box_aspect(1/box_ratio)
   ax2.legend(ncol = 1, labels=['mobility'], loc=loc2, \
   labelspacing = 0.5, labelcolor = '#ffffff', edgecolor = 'none', facecolor = 'none', fancybox = False)


   present_day_TTG=2.171e22
   if current_model == '12_Qaw20sat25r3':
      mpl.rcParams['figure.figsize'] = [4.5,1.65] # (1 panel, width x height, 1=600px with dpi=300)
      fig, ax = plt.subplots(2, 1, squeeze=False)
      ax[0][0].plot(time_Myr_run[run],velocity_rms_run[run], linestyle='-', linewidth=1, alpha=1, color = '#000000')
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,180.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],180.,210.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],210.,240.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],240.,290.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],290.,410.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],410.,440.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],440.,490.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],490.,665.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],665.,785.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],785.,815.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],815.,980.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],980.,1410.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],1410.,1500.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].axvline(36.93, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(478.98, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(684.24, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1224.83, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper right'
      loc2='lower right'
      evolved_time=1500.  
      time_ticks=[0,300,600,900,1200,1500]
      ax[0][0].set_xlim(0., evolved_time)
      ax[0][0].set_xscale('linear')
      ax[0][0].set_xticks(time_ticks)
      ax[0][0].set_ylim(min_vrms,max_vrms)
      ax[0][0].set_yscale('log')
      ax[0][0].set_yticks(vrms_ticks)
      ax[0][0].set_ylabel('$v_{rms}$ (cm/yr)',labelpad=4)
      ax[0][0].set_box_aspect(1/box_ratio)
      ax[0][0].legend(ncol = 1, labels=['v$_{rms}$'], loc=loc1, \
      labelspacing = 0.5, edgecolor = 'none', facecolor = 'none', fancybox = False)  

      ax2 = ax[0][0].twinx()
      ax2.plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=0.7, color = '#ffffff')
      ax2.set_xlim(0., evolved_time)
      ax2.set_xscale('linear')
      ax2.set_xticks(time_ticks)
      ax2.set_yscale('linear')
      ax2.set_ylim(0.,max_mobility)
      ax2.set_yticks(mobility_ticks)
      ax2.set_ylabel('Mobility',labelpad=4)
      ax2.set_box_aspect(1/box_ratio)
      ax2.legend(ncol = 1, labels=['mobility'], loc=loc2, \
      labelspacing = 0.5, labelcolor = '#ffffff', edgecolor = 'none', facecolor = 'none', fancybox = False)

      ax[1][0].plot(time_Myr_step_run[run],TotalTTGmass_run[run]/present_day_TTG, linestyle='-', linewidth=1, alpha=1, color = '#000000')
      ax[1][0].fill_betweenx(time_Myr_run[run],0.,469.,  linewidth=0,color = '#AE027E',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],469.,900.,  linewidth=0,color = '#E61D8C',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],900.,1300.,  linewidth=0,color = '#F668B2',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],1300.,1500.,  linewidth=0,color = '#F881B5',alpha=1,edgecolor='none')
      ax[1][0].axvline(36.93, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(478.98, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(684.24, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(1224.83, color='#000000', linestyle='--',linewidth=1, alpha=0.5)       
      ax[1][0].set_xlim(0., evolved_time)
      ax[1][0].set_xscale('linear')
      ax[1][0].set_xticks(time_ticks)
      ax[1][0].set_ylim(0,1.5)
      ax[1][0].set_yscale('linear')
      ax[1][0].set_yticks([0.5,1,1.5])
      ax[1][0].set_ylabel('TTG (CCM)',labelpad=6)
      ax[1][0].set_xlabel('Time (Myr)')
      ax[1][0].set_box_aspect(1/box_ratio)
   elif current_model == '16_Qf02aw20sat40r2':
      mpl.rcParams['figure.figsize'] = [4.5,1.65] # (1 panel, width x height, 1=600px with dpi=300)       
      fig, ax = plt.subplots(2, 1, squeeze=False)
      ax[0][0].plot(time_Myr_run[run],velocity_rms_run[run], linestyle='-', linewidth=1, alpha=1, color = '#000000')         
      ax[0][0].fill_betweenx(time_Myr_run[run],0.,100.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],100.,300.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],300.,330.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],330.,440.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],440.,460.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],460.,670.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],670.,700.,  linewidth=0,color = '#024b79',alpha=0.8,edgecolor='none')
      ax[0][0].fill_betweenx(time_Myr_run[run],700.,1500.,  linewidth=0,color = '#ffad48',alpha=1,edgecolor='none')
      ax[0][0].axvline(37.47, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(201.05, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(688.41, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[0][0].axvline(1234.58, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      loc1='upper right'
      loc2='lower right'
      evolved_time=1500.  
      time_ticks=[0,300,600,900,1200,1500]
      ax[0][0].set_xlim(0., evolved_time)
      ax[0][0].set_xscale('linear')
      ax[0][0].set_xticks(time_ticks)
      ax[0][0].set_ylim(min_vrms,max_vrms)
      ax[0][0].set_yscale('log')
      ax[0][0].set_yticks(vrms_ticks)
      ax[0][0].set_ylabel('$v_{rms}$ (cm/yr)',labelpad=4)
      ax[0][0].set_box_aspect(1/box_ratio)
      ax[0][0].legend(ncol = 1, labels=['v$_{rms}$'], loc=loc1, \
      labelspacing = 0.5, edgecolor = 'none', facecolor = 'none', fancybox = False)

      ax2 = ax[0][0].twinx()
      ax2.plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=0.7, color = '#ffffff')
      ax2.set_xlim(0., evolved_time)
      ax2.set_xscale('linear')
      ax2.set_xticks(time_ticks)
      ax2.set_yscale('linear')
      ax2.set_ylim(0.,max_mobility)
      ax2.set_yticks(mobility_ticks)
      ax2.set_ylabel('Mobility',labelpad=4)
      ax2.set_box_aspect(1/box_ratio)
      ax2.legend(ncol = 1, labels=['mobility'], loc=loc2, \
      labelspacing = 0.5, labelcolor = '#ffffff', edgecolor = 'none', facecolor = 'none', fancybox = False)

      ax[1][0].plot(time_Myr_step_run[run],TotalTTGmass_run[run]/present_day_TTG, linestyle='-', linewidth=1, alpha=1, color = '#000000')
      ax[1][0].fill_betweenx(time_Myr_run[run],0.,469.,  linewidth=0,color = '#AE027E',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],469.,900.,  linewidth=0,color = '#E61D8C',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],900.,1300.,  linewidth=0,color = '#F668B2',alpha=1,edgecolor='none')
      ax[1][0].fill_betweenx(time_Myr_run[run],1300.,1500.,  linewidth=0,color = '#F881B5',alpha=1,edgecolor='none')       
      ax[1][0].axvline(37.47, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(201.05, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(688.41, color='#000000', linestyle='--',linewidth=1, alpha=0.5)
      ax[1][0].axvline(1234.58, color='#000000', linestyle='--',linewidth=1, alpha=0.5)       
      ax[1][0].set_xlim(0., evolved_time)
      ax[1][0].set_xscale('linear')
      ax[1][0].set_xticks(time_ticks)
      ax[1][0].set_ylim(0,1.5)
      ax[1][0].set_yscale('linear')
      ax[1][0].set_yticks([0.5,1,1.5])
      ax[1][0].set_ylabel('TTG (CCM)',labelpad=6)
      ax[1][0].set_xlabel('Time (Myr)')
      ax[1][0].set_box_aspect(1/box_ratio)

   fig.tight_layout()

   file_name='fig'+str(batch)+'_vrms_'+current_model+'.pdf'
   plt.savefig(here+'/../figures/'+file_name)

   run = run + 1
