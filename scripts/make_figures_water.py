# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# updated on 10 Sep 2026
# to make figures for the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 2. Contrastive effects of water
# submitted to JGR: Solid Earth

# This script uses StagPy library to post-process StagYY data.
# Please cite as:
# Adrien Morison, Martina Ulvrova, Stephane Labrosse, B4rsh, theofatou, & tfrass49. (2024).
# StagPython/StagPy: v0.20.1 (v0.20.1). Zenodo. https://doi.org/10.5281/zenodo.13684760

# Other python packages required: os, math, numpy, h5py, pandas, matplotlib

# Once required packages are installed, execute the script from terminal with command "python make_figures_water.py"

import os
import math
import numpy as np
import h5py
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from stagpy import stagyydata


### beginning of input ###
# user only has to input figure number (2/4/6/8) as batch
batch=2
### end of input ###


if batch == 2:
   model_folder = 'Table4_models/'
   model_list = ["01_AfcT19r1/", "02_AfcT19r2/", "03_AT19r1/", "04_AT19r2/", "05_AT19r3/", "06_AT19r4/", "07_AT19r5/", "08_AT19r6/", "09_AT19r7/", \
               "10_AT19r8/", "11_AT19FBr1/", "12_AT19FBr2/", "13_AT19FBr3/", "14_AT18FBr1/", "15_AT18FBr2/", "16_AT18FBr3/", "17_AT17FBr1/", "18_AT17FBr2/", "19_AT17FBr3/"]
   palette = ['#00647A','#C25B60','#C25B60','#FFBF70','#43B7C2','#000000']
   labels_a = ('1900 K, constant $\\mu$ (2x)','1900 K (8x)','1900 K, FB (3x)','1825 K, FB (3x)','1750 K, FB (3x)','AT19r1')
   figure_title = ('Full annulus | different $\\it{T}$$_{P0}$ | initial $\\mu$ = 0.1, which adapts with time')
   location='lower center'
   style = ['solid','solid','dashed','dashed','dashed','dotted']
   symbol = ['^','^','o','o','o']
   legend_columns = 3
   legend_panel=2

elif batch == 4:
   model_folder = 'Table5_models/'
   model_list = ["01_QW1r1/", "02_QW1r2/", "03_QW1r3/", "04_QW2r1/", "05_QW2r2/", "06_QW2r3/", "07_QW5r1/", "08_QW5r2/", "09_QW5r3/", \
               "10_QW10r1/", "11_QW10r2/", "12_QW10r3/", "13_QW10hr1/", "14_QW10hr2/", "15_QW10hr3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#00647A','#00647A','#A800FF','#f19106'] #veronica
   labels_a = ('1 OM, 500 ppm (3x)','2 OM, 500 ppm (3x)','5 OM, 500 ppm (3x)','10 OM, 500 ppm (3x)','10 OM, 5000 ppm (3x)','QW5r1','QW10hr1')
   figure_title = ('Quadrant annulus | different water content (OM) and requirement for TTG formation')
   location='upper center'
   style = ['solid','solid','solid','solid','dashed','dashed','dashed']
   symbol = ['o','o','o','o','D','o','^']
   legend_columns = 4
   legend_panel = 2

elif batch == 6:
   model_folder = 'Table6_models/'
   model_list = ["01_Qaw15r1/", "02_Qaw15r2/", "03_Qaw15r3/", "04_Qaw20r1/", "05_Qaw20r2/", "06_Qaw20r3/", \
               "07_Qaw20sat20r1/", "08_Qaw20sat20r2/", "09_Qaw20sat20r3/", "10_Qaw20sat25r1/", "11_Qaw20sat25r2/", "12_Qaw20sat25r3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#00647A','#A800FF','#f19106'] #veronica
   labels_a = ('a$_{w}$ = 1.5, 100$\\%$ saturation (3x)','a$_{w}$ = 2.0, 100$\\%$ saturation (3x)','a$_{w}$ = 2.0, 20$\\%$ saturation (3x)','a$_{w}$ = 2.0, 25$\\%$ saturation (3x)','Qaw20r3','Qaw20sat20r1')
   figure_title = ('Quadrant annulus | 5 OM, 5000 ppm | water-dependent density')
   location='upper center'
   style = ['solid','solid','solid','solid','dashed','dashed']
   symbol = ['o','o','^','^']
   legend_columns = 3
   legend_panel = 2

elif batch == 8:
   model_folder = 'Table6_models/'
   model_list = ["17_Qsat20H2Oetar1/", "18_Qsat20H2Oetar2/", "19_Qsat20H2Oetar3/", "20_Qsat25H2Oetar1/", "21_Qsat25H2Oetar2/", "22_Qsat25H2Oetar3/", \
                 "06_Qaw20r3/", "07_Qaw20sat20r1/"]
   palette = ['#C25B60','#FFBF70','#522680','#f19106','#43B7C2','#00647A'] #tekhelet
   labels_a = ('Qaw20r3 (no effect on viscosity)','Qaw20sat20r1 (no effect on viscosity)','20$\\%$ saturation (3x)','25$\\%$ saturation (3x)','Qsat20H2Oetar3','Qsat25H2Oetar2')
   figure_title = ('Quadrant annulus | 5 OM, 5000 ppm')
   location='upper center'
   style = ['solid','solid','dashed','dashed','dashed','dashed']
   symbol = ['o','^']
   legend_columns = 3
   legend_panel = 1

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
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.titlesize'] = 11
mpl.rcParams['axes.titleweight'] = "bold"
mpl.rcParams['axes.labelweight'] = "normal"
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['figure.dpi'] = 600
mpl.rcParams['figure.figsize'] = [7,6] # (3 panels, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams["pdf.use14corefonts"] = True
box_ratio=4.5
evolved_time = 4500.
time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]
minTin=1150.
maxTin=1950.
Tin_ticks=[1150,1350,1550,1750,1950]
max_mobility = 3.0
mobility_ticks=[0.0,1.0,2.0,3.0]
time_xaxis_interpolation=4500

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

tp_run      = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
temp200_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
SD200_run   = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
temp400_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
SD400_run   = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

time_Myr_step_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
TotalTTGmass_run = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

# for interpolated values
tp_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_run_even_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
TotalTTGmass_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
tmean_run_interp = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
moving_average_mobility_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
plateness_batch = np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))


seconds_in_year = 31556926
seconds_in_Myear = seconds_in_year*1e6
earth_radius = 6371.e3
core_radius = 3481.e3
density_TTG = 2700.

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

load_internalT = True
compute_TTGmass = True

run = int(0)
for path in model_list:

   volume_mantle = (4 * math.pi*(earth_radius**3-core_radius**3))/3
   # quadrant analytical volume is 1/16 of full annulus, which in itself is pi/ny times actual volume
   if batch == 2:
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

   # from stagpy import stagyydata
   sdat = stagyydata.StagyyData(here+'/../'+str(model_folder)+path)

   # extracting time from time_botT.h5
   # this is needed when plotting the time evolution of a field
   time_list = []
   coreT_list = []
   f = h5py.File(here+'/../'+str(model_folder)+path+'+hdf5/time_botT.h5', 'r')

   ndatasets = len(f.keys())
   for key in f.keys():
      data = f[key]
      time_list.append(float(data[0]))
      coreT_list.append(float(data[1]))
   f.close()

   time=np.array(time_list)
   coreT=np.array(coreT_list)
   time_Myr=time/seconds_in_Myear
   print('Final core T (K):        '+str('{:.0f}'.format(coreT[len(coreT)-1])))
   # print('# HDF datasets: '+str(ndatasets))

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
   # print('# frames (rprof.dat): ' + str(nframes))
   print('Evolution time (Myr):    '+str('{:.2f}'.format(sdat.tseries.time[len(sdat.tseries.time)-1]*1e-6/seconds_in_year)))
   # timesteps for initial and final states
   step1=snaps[1]
   step2=snaps[nframes-1]



   # loading Tmean from time series
   tmean_list = []

   for i in snaps:
      tmean_list.append(float(sdat.tseries['Tmean'].values[i]))
   tmean_run[run] = np.array(tmean_list)
   tmean_list.clear()
   print('Final mantle T (K):      '+str('{:.0f}'.format(tmean_run[run][len(tmean_run[run])-1])))   

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

   mobility_batch[run] = mobility_run[run]
   time_Myr_mobility_batch[run] = time_Myr_mobility_run[run]
   moving_average_mobility_batch[run] = moving_average_mobility_run[run]
   plateness_batch[run] = plateness_run[run]
        
   # interpolating mobility and plateness on a much finer time axis with 4500 elements
   time_interp2 = np.linspace(0,time_xaxis_interpolation,num=time_xaxis_interpolation+1)
   mobility_run_even_interp[run] = np.interp(time_interp2,time_Myr_mobility_batch[run],moving_average_mobility_batch[run])


   # computing total mass of TTG produced
   if compute_TTGmass:
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

      print('TTG mass (scaled, kg):   '+str('{:.2e}'.format(TTGmass_cumulative[len(time_TTGmass)-1]*factor)))

      # getting cumulative volume
      TTGvolume = TTGmass_cumulative * factor * 1e-9 / density_TTG

      timesteps.clear()
      time_TTGmass.clear()
      TTG1.clear()
      TTG2.clear()
      TTG3.clear()

      time_Myr_step_run[run] = time_Myr_step
      TotalTTGmass_run[run] = TTGmass_cumulative * factor

      # interpolating TTG mass on a much finer time axis with 1500/4500 elements
      time_interp3 = np.linspace(0,time_xaxis_interpolation,num=(time_xaxis_interpolation+1)*10)
      TotalTTGmass_run_interp[run] = np.interp(time_interp3,time_Myr_step_run[run],TotalTTGmass_run[run])


   if load_internalT:

      tp_list = []
      t200_list = []
      sd200_list = []
      t400_list = []
      sd400_list = []

      file = (here+'/../'+str(model_folder)+path+'potentialT.dat')
      f = open(file,'r')
      lines=f.readlines()
      f.close
      nlines=int(len(lines))

      for i in np.arange(0,nlines):
         p = lines[i].split()
         tp_list.append(float(p[1]))
         t200_list.append(float(p[2]))
         sd200_list.append(float(p[3]))
         t400_list.append(float(p[4]))
         sd400_list.append(float(p[5]))
      
      tp_run[run]      = np.array(tp_list)
      temp200_run[run] = np.array(t200_list)
      SD200_run[run]   = np.array(sd200_list)
      temp400_run[run] = np.array(t400_list)
      SD400_run[run]   = np.array(sd400_list)

      tp_list.clear()
      t200_list.clear()
      sd200_list.clear()
      t400_list.clear()
      sd400_list.clear()

      # interpolating Tp on a much finer time axis with 4500 elements
      time_interp4 = np.linspace(0,time_xaxis_interpolation,num=time_xaxis_interpolation+1)
      tp_run_interp[run] = np.interp(time_interp4,time_Myr_run[run],(temp200_run[run]+SD200_run[run]))
      print('Final pot. T (deg C):    '+str('{:.0f}'.format(temp200_run[run][len(temp200_run[run])-1]+SD200_run[run][len(SD200_run[run])-1])))


   # plotting potential temperature, mobility (moving average), TTG mass per model 
   fig, ax = plt.subplots(3, 1, figsize = (7,6), squeeze=False, dpi = 600)

   ax[0][0].plot(time_Myr_run[run],temp200_run[run]+SD200_run[run], linestyle='-', linewidth=1, alpha=1, color = '#C25B60')
   ax[0][0].set_xlim(0., evolved_time)
   ax[0][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)
   ax[0][0].set_ylim(minTin,maxTin)
   ax[0][0].set_yscale('linear')
   ax[0][0].set_yticks(Tin_ticks)   
   ax[0][0].set_ylabel('Potential T ($\\degree$C)')
   ax[0][0].set_box_aspect(1/box_ratio)

   ax[1][0].plot(time_Myr_run[run],moving_average_mobility_run[run], linestyle='-', linewidth=1, alpha=1, color = '#FFBF70')
   ax[1][0].set_xlim(0., evolved_time)
   ax[1][0].set_xscale('linear')
   ax[0][0].set_xticks(time_ticks)
   ax[1][0].set_ylim(0., max_mobility)
   ax[1][0].set_yscale('linear')
   ax[1][0].set_yticks(mobility_ticks)   
   ax[1][0].set_ylabel('Mobility',labelpad=23)
   ax[1][0].set_box_aspect(1/box_ratio)

   ax[2][0].plot(time_Myr_step_run[run],TotalTTGmass_run[run], linestyle='-', linewidth=1, alpha=1, color = '#43B7C2')
   ax[2][0].set_xlim(0., evolved_time)
   ax[2][0].set_xscale('linear')
   ax[2][0].set_xticks(time_ticks)
   ax[2][0].set_ylim(0., 8.e22)
   ax[2][0].set_yscale('linear')
   ax[2][0].set_ylabel('TTG mass (kg)',labelpad=23)
   ax[2][0].set_xlabel('Time (Myr)')
   ax[2][0].set_box_aspect(1/box_ratio)

   # saving it in the model folder
   file_name='model_'+str(path[0:len(path)-1])+'_diagnostics.png'
   fig.tight_layout()
   plt.savefig(here+'/../'+str(model_folder)+path+file_name)


   run = run + 1


# calculating range and mean from all models
if batch == 2:
   # potential T
   maxTin1900_cf = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1]])
   minTin1900_cf = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1]])
   meanTin1900_cf = np.mean(np.array([maxTin1900_cf,minTin1900_cf]),axis=0)

   maxTin1900 = np.maximum.reduce([tp_run_interp[2],tp_run_interp[3],tp_run_interp[4], \
                                 tp_run_interp[5],tp_run_interp[6],tp_run_interp[7], \
                                 tp_run_interp[8],tp_run_interp[9]])
   minTin1900 = np.minimum.reduce([tp_run_interp[2],tp_run_interp[3],tp_run_interp[4], \
                                 tp_run_interp[5],tp_run_interp[6],tp_run_interp[7], \
                                 tp_run_interp[8],tp_run_interp[9]])    
   meanTin1900 = np.mean(np.array([maxTin1900,minTin1900]),axis=0)


   maxTin1900FB = np.maximum.reduce([tp_run_interp[10],tp_run_interp[11],tp_run_interp[12]])
   minTin1900FB = np.minimum.reduce([tp_run_interp[10],tp_run_interp[11],tp_run_interp[12]])
   meanTin1900FB = np.mean(np.array([maxTin1900FB,minTin1900FB]),axis=0)

   maxTin1825FB = np.maximum.reduce([tp_run_interp[13],tp_run_interp[14],tp_run_interp[15]])
   minTin1825FB = np.minimum.reduce([tp_run_interp[13],tp_run_interp[14],tp_run_interp[15]])
   meanTin1825FB = np.mean(np.array([maxTin1825FB,minTin1825FB]),axis=0)

   maxTin1750FB = np.maximum.reduce([tp_run_interp[16],tp_run_interp[17],tp_run_interp[18]])
   minTin1750FB = np.minimum.reduce([tp_run_interp[16],tp_run_interp[17],tp_run_interp[18]])
   meanTin1750FB = np.mean(np.array([maxTin1750FB,minTin1750FB]),axis=0)
      
   # Mobility M
   maxM1900_cf = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1]])
   minM1900_cf = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1]])
   meanM1900_cf = np.mean(np.array([maxM1900_cf,minM1900_cf]),axis=0)

   maxM1900 = np.maximum.reduce([mobility_run_even_interp[2],mobility_run_even_interp[3],mobility_run_even_interp[4], \
                                 mobility_run_even_interp[5],mobility_run_even_interp[6],mobility_run_even_interp[7], \
                                 mobility_run_even_interp[8],mobility_run_even_interp[9]])
   minM1900 = np.minimum.reduce([mobility_run_even_interp[2],mobility_run_even_interp[3],mobility_run_even_interp[4], \
                                 mobility_run_even_interp[5],mobility_run_even_interp[6],mobility_run_even_interp[7], \
                                 mobility_run_even_interp[8],mobility_run_even_interp[9]])    
   meanM1900 = np.mean(np.array([maxM1900,minM1900]),axis=0)

   maxM1900FB = np.maximum.reduce([mobility_run_even_interp[10],mobility_run_even_interp[11],mobility_run_even_interp[12]])
   minM1900FB = np.minimum.reduce([mobility_run_even_interp[10],mobility_run_even_interp[11],mobility_run_even_interp[12]])
   meanM1900FB = np.mean(np.array([maxM1900FB,minM1900FB]),axis=0)

   maxM1825FB = np.maximum.reduce([mobility_run_even_interp[13],mobility_run_even_interp[14],mobility_run_even_interp[15]])
   minM1825FB = np.minimum.reduce([mobility_run_even_interp[13],mobility_run_even_interp[14],mobility_run_even_interp[15]])
   meanM1825FB = np.mean(np.array([maxM1825FB,minM1825FB]),axis=0)

   maxM1750FB = np.maximum.reduce([mobility_run_even_interp[16],mobility_run_even_interp[17],mobility_run_even_interp[18]])
   minM1750FB = np.minimum.reduce([mobility_run_even_interp[16],mobility_run_even_interp[17],mobility_run_even_interp[18]])
   meanM1750FB = np.mean(np.array([maxM1750FB,minM1750FB]),axis=0)  

   # TTG mass
   maxTTG1900_cf = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1]])
   minTTG1900_cf = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1]])
   meanTTG1900_cf = np.mean(np.array([maxTTG1900_cf,minTTG1900_cf]),axis=0)

   maxTTG1900 = np.maximum.reduce([TotalTTGmass_run_interp[2],TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4], \
                                 TotalTTGmass_run_interp[5],TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7], \
                                 TotalTTGmass_run_interp[8],TotalTTGmass_run_interp[9]])
   minTTG1900 = np.minimum.reduce([TotalTTGmass_run_interp[2],TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4], \
                                 TotalTTGmass_run_interp[5],TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7], \
                                 TotalTTGmass_run_interp[8],TotalTTGmass_run_interp[9]])    
   meanTTG1900 = np.mean(np.array([maxTTG1900,minTTG1900]),axis=0)

   maxTTG1900FB = np.maximum.reduce([TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11],TotalTTGmass_run_interp[12]])
   minTTG1900FB = np.minimum.reduce([TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11],TotalTTGmass_run_interp[12]])
   meanTTG1900FB = np.mean(np.array([maxTTG1900FB,minTTG1900FB]),axis=0)

   maxTTG1825FB = np.maximum.reduce([TotalTTGmass_run_interp[13],TotalTTGmass_run_interp[14],TotalTTGmass_run_interp[15]])
   minTTG1825FB = np.minimum.reduce([TotalTTGmass_run_interp[13],TotalTTGmass_run_interp[14],TotalTTGmass_run_interp[15]])
   meanTTG1825FB = np.mean(np.array([maxTTG1825FB,minTTG1825FB]),axis=0)

   maxTTG1750FB = np.maximum.reduce([TotalTTGmass_run_interp[16],TotalTTGmass_run_interp[17],TotalTTGmass_run_interp[18]])
   minTTG1750FB = np.minimum.reduce([TotalTTGmass_run_interp[16],TotalTTGmass_run_interp[17],TotalTTGmass_run_interp[18]])
   meanTTG1750FB = np.mean(np.array([maxTTG1750FB,minTTG1750FB]),axis=0)

elif batch == 4:
   # potential T
   maxTin1OM_DE = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   minTin1OM_DE = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   meanTin1OM_DE = np.mean(np.array([maxTin1OM_DE,minTin1OM_DE]),axis=0)

   maxTin2OM_DE = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   minTin2OM_DE = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   meanTin2OM_DE = np.mean(np.array([maxTin2OM_DE,minTin2OM_DE]),axis=0)

   maxTin5OM_DE = np.maximum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   minTin5OM_DE = np.minimum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   meanTin5OM_DE = np.mean(np.array([maxTin5OM_DE,minTin5OM_DE]),axis=0)

   maxTin10OM_DE = np.maximum.reduce([tp_run_interp[9],tp_run_interp[10],tp_run_interp[11]])
   minTin10OM_DE = np.minimum.reduce([tp_run_interp[9],tp_run_interp[10],tp_run_interp[11]])
   meanTin10OM_DE = np.mean(np.array([maxTin10OM_DE,minTin10OM_DE]),axis=0)

   maxTin10OM_5000 = np.maximum.reduce([tp_run_interp[12],tp_run_interp[13],tp_run_interp[14]])
   minTin10OM_5000 = np.minimum.reduce([tp_run_interp[12],tp_run_interp[13],tp_run_interp[14]])
   meanTin10OM_5000 = np.mean(np.array([maxTin10OM_5000,minTin10OM_5000]),axis=0)

   # Mobility M
   maxM1OM_DE = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   minM1OM_DE = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   meanM1OM_DE = np.mean(np.array([maxM1OM_DE,minM1OM_DE]),axis=0)

   maxM2OM_DE = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   minM2OM_DE = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   meanM2OM_DE = np.mean(np.array([maxM2OM_DE,minM2OM_DE]),axis=0)

   maxM5OM_DE = np.maximum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   minM5OM_DE = np.minimum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   meanM5OM_DE = np.mean(np.array([maxM5OM_DE,minM5OM_DE]),axis=0)

   maxM10OM_DE = np.maximum.reduce([mobility_run_even_interp[9],mobility_run_even_interp[10],mobility_run_even_interp[11]])
   minM10OM_DE = np.minimum.reduce([mobility_run_even_interp[9],mobility_run_even_interp[10],mobility_run_even_interp[11]])
   meanM10OM_DE = np.mean(np.array([maxM10OM_DE,minM10OM_DE]),axis=0)

   maxM10OM_5000 = np.maximum.reduce([mobility_run_even_interp[12],mobility_run_even_interp[13],mobility_run_even_interp[14]])
   minM10OM_5000 = np.minimum.reduce([mobility_run_even_interp[12],mobility_run_even_interp[13],mobility_run_even_interp[14]])
   meanM10OM_5000 = np.mean(np.array([maxM10OM_5000,minM10OM_5000]),axis=0)

   # TTG mass
   maxTTG1OM_DE = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   minTTG1OM_DE = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   meanTTG1OM_DE = np.mean(np.array([maxTTG1OM_DE,minTTG1OM_DE]),axis=0)

   maxTTG2OM_DE = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   minTTG2OM_DE = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   meanTTG2OM_DE = np.mean(np.array([maxTTG2OM_DE,minTTG2OM_DE]),axis=0)

   maxTTG5OM_DE = np.maximum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   minTTG5OM_DE = np.minimum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   meanTTG5OM_DE = np.mean(np.array([maxTTG5OM_DE,minTTG5OM_DE]),axis=0)

   maxTTG10OM_DE = np.maximum.reduce([TotalTTGmass_run_interp[9],TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11]])
   minTTG10OM_DE = np.minimum.reduce([TotalTTGmass_run_interp[9],TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11]])
   meanTTG10OM_DE = np.mean(np.array([maxTTG10OM_DE,minTTG10OM_DE]),axis=0)

   maxTTG10OM_5000 = np.maximum.reduce([TotalTTGmass_run_interp[12],TotalTTGmass_run_interp[13],TotalTTGmass_run_interp[14]])
   minTTG10OM_5000 = np.minimum.reduce([TotalTTGmass_run_interp[12],TotalTTGmass_run_interp[13],TotalTTGmass_run_interp[14]])
   meanTTG10OM_5000 = np.mean(np.array([maxTTG10OM_5000,minTTG10OM_5000]),axis=0)

elif batch == 6:
   # potential T
   maxTinrho15 = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   minTinrho15 = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   meanTinrho15 = np.mean(np.array([maxTinrho15,minTinrho15]),axis=0)

   maxTinrho20 = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   minTinrho20 = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   meanTinrho20 = np.mean(np.array([maxTinrho20,minTinrho20]),axis=0)

   maxTinrho20alpha20 = np.maximum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   minTinrho20alpha20 = np.minimum.reduce([tp_run_interp[6],tp_run_interp[7],tp_run_interp[8]])
   meanTinrho20alpha20 = np.mean(np.array([maxTinrho20alpha20,minTinrho20alpha20]),axis=0)  

   maxTinrho20alpha25 = np.maximum.reduce([tp_run_interp[9],tp_run_interp[10],tp_run_interp[11]])
   minTinrho20alpha25 = np.minimum.reduce([tp_run_interp[9],tp_run_interp[10],tp_run_interp[11]])
   meanTinrho20alpha25 = np.mean(np.array([maxTinrho20alpha25,minTinrho20alpha25]),axis=0)

   # Mobility M
   maxMrho15 = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   minMrho15 = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   meanMrho15 = np.mean(np.array([maxMrho15,minMrho15]),axis=0)

   maxMrho20 = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   minMrho20 = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   meanMrho20 = np.mean(np.array([maxMrho20,minMrho20]),axis=0)

   maxMrho20alpha20 = np.maximum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   minMrho20alpha20 = np.minimum.reduce([mobility_run_even_interp[6],mobility_run_even_interp[7],mobility_run_even_interp[8]])
   meanMrho20alpha20 = np.mean(np.array([maxMrho20alpha20,minMrho20alpha20]),axis=0)

   maxMrho20alpha25 = np.maximum.reduce([mobility_run_even_interp[9],mobility_run_even_interp[10],mobility_run_even_interp[11]])
   minMrho20alpha25 = np.minimum.reduce([mobility_run_even_interp[9],mobility_run_even_interp[10],mobility_run_even_interp[11]])
   meanMrho20alpha25 = np.mean(np.array([maxMrho20alpha25,minMrho20alpha25]),axis=0)

   # TTG mass
   maxTTGrho15 = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   minTTGrho15 = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   meanTTGrho15 = np.mean(np.array([maxTTGrho15,minTTGrho15]),axis=0)

   maxTTGrho20 = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   minTTGrho20 = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   meanTTGrho20 = np.mean(np.array([maxTTGrho20,minTTGrho20]),axis=0)

   maxTTGrho20alpha20 = np.maximum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   minTTGrho20alpha20 = np.minimum.reduce([TotalTTGmass_run_interp[6],TotalTTGmass_run_interp[7],TotalTTGmass_run_interp[8]])
   meanTTGrho20alpha20 = np.mean(np.array([maxTTGrho20alpha20,minTTGrho20alpha20]),axis=0)

   maxTTGrho20alpha25 = np.maximum.reduce([TotalTTGmass_run_interp[9],TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11]])
   minTTGrho20alpha25 = np.minimum.reduce([TotalTTGmass_run_interp[9],TotalTTGmass_run_interp[10],TotalTTGmass_run_interp[11]])
   meanTTGrho20alpha25 = np.mean(np.array([maxTTGrho20alpha25,minTTGrho20alpha25]),axis=0)

elif batch == 8:
   # potential T
   maxTinrho0alpha20 = np.maximum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   minTinrho0alpha20 = np.minimum.reduce([tp_run_interp[0],tp_run_interp[1],tp_run_interp[2]])
   meanTinrho0alpha20 = np.mean(np.array([maxTinrho0alpha20,minTinrho0alpha20]),axis=0)

   maxTinrho0alpha25 = np.maximum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   minTinrho0alpha25 = np.minimum.reduce([tp_run_interp[3],tp_run_interp[4],tp_run_interp[5]])
   meanTinrho0alpha25 = np.mean(np.array([maxTinrho0alpha25,minTinrho0alpha25]),axis=0)

   # Mobility M
   maxMrho0alpha20 = np.maximum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   minMrho0alpha20 = np.minimum.reduce([mobility_run_even_interp[0],mobility_run_even_interp[1],mobility_run_even_interp[2]])
   meanMrho0alpha20 = np.mean(np.array([maxMrho0alpha20,minMrho0alpha20]),axis=0)

   maxMrho0alpha25 = np.maximum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   minMrho0alpha25 = np.minimum.reduce([mobility_run_even_interp[3],mobility_run_even_interp[4],mobility_run_even_interp[5]])
   meanMrho0alpha25 = np.mean(np.array([maxMrho0alpha25,minMrho0alpha25]),axis=0)   

   # TTG mass
   maxTTGrho0alpha20 = np.maximum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   minTTGrho0alpha20 = np.minimum.reduce([TotalTTGmass_run_interp[0],TotalTTGmass_run_interp[1],TotalTTGmass_run_interp[2]])
   meanTTGrho0alpha20 = np.mean(np.array([maxTTGrho0alpha20,minTTGrho0alpha20]),axis=0)

   maxTTGrho0alpha25 = np.maximum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   minTTGrho0alpha25 = np.minimum.reduce([TotalTTGmass_run_interp[3],TotalTTGmass_run_interp[4],TotalTTGmass_run_interp[5]])
   meanTTGrho0alpha25 = np.mean(np.array([maxTTGrho0alpha25,minTTGrho0alpha25]),axis=0)

else:
    print('No interpolation done')




fig, ax = plt.subplots(3, 1, squeeze=False)

# Potential temperature   
if batch == 2:
   ax[0][0].plot(time_interp4,meanTin1900_cf, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin1900   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin1900FB , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin1825FB , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin1750FB , linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_Myr_run[2],temp200_run[2]+SD200_run[2], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[0][0].fill_between(time_interp4,minTin1900_cf,maxTin1900_cf,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin1900,maxTin1900      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin1900FB,maxTin1900FB  ,linewidth=0,color = palette[2],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin1825FB,maxTin1825FB  ,linewidth=0,color = palette[4],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin1750FB,maxTin1750FB  ,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 4:
   ax[0][0].plot(time_interp4,meanTin1OM_DE   , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')        
   ax[0][0].plot(time_interp4,meanTin2OM_DE   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTin5OM_DE   , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')  
   ax[0][0].plot(time_interp4,meanTin10OM_DE  , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')            
   ax[0][0].plot(time_interp4,meanTin10OM_5000, linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_Myr_run[6],temp200_run[6]+SD200_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[0][0].plot(time_Myr_run[12],temp200_run[12]+SD200_run[12], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])
   ax[0][0].fill_between(time_interp4,minTin1OM_DE,maxTin1OM_DE      ,linewidth=0,color = palette[0],alpha=0.2)        
   ax[0][0].fill_between(time_interp4,minTin2OM_DE,maxTin2OM_DE      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin5OM_DE,maxTin5OM_DE      ,linewidth=0,color = palette[2],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin10OM_DE,maxTin10OM_DE    ,linewidth=0,color = palette[3],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTin10OM_5000,maxTin10OM_5000,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 6:
   ax[0][0].plot(time_interp4,meanTinrho15       , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTinrho20       , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')            
   ax[0][0].plot(time_interp4,meanTinrho20alpha20, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTinrho20alpha25, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_Myr_run[5],temp200_run[5]+SD200_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[0][0].plot(time_Myr_run[6],temp200_run[6]+SD200_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[0][0].fill_between(time_interp4,minTinrho15,maxTinrho15        ,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTinrho20,maxTinrho20        ,linewidth=0,color = palette[1],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTinrho20alpha20,maxTinrho20alpha20        ,linewidth=0,color = palette[2],alpha=0.2)            
   ax[0][0].fill_between(time_interp4,minTinrho20alpha25,maxTinrho20alpha25        ,linewidth=0,color = palette[3],alpha=0.2)   
elif batch == 8:
   ax[0][0].plot(time_Myr_run[6],temp200_run[6]+SD200_run[6], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[0][0].plot(time_Myr_run[7],temp200_run[7]+SD200_run[7], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])   
   ax[0][0].plot(time_interp4,meanTinrho0alpha20     , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_interp4,meanTinrho0alpha25     , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[0][0].plot(time_Myr_run[2],temp200_run[2]+SD200_run[2], linestyle=style[2], linewidth=1, alpha=1, color = palette[2])
   ax[0][0].plot(time_Myr_run[4],temp200_run[4]+SD200_run[4], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[0][0].fill_between(time_interp4,minTinrho0alpha20,maxTinrho0alpha20          ,linewidth=0,color = palette[0],alpha=0.2)
   ax[0][0].fill_between(time_interp4,minTinrho0alpha25,maxTinrho0alpha25          ,linewidth=0,color = palette[1],alpha=0.2)

ax[0][0].set_xlim(0., evolved_time)
ax[0][0].set_xscale('linear')
ax[0][0].set_xticks(time_ticks)
ax[0][0].set_ylim(minTin,maxTin)
ax[0][0].set_yscale('linear')
ax[0][0].set_yticks(Tin_ticks)
ax[0][0].set_ylabel('Potential T ($\\degree$C)')
ax[0][0].set_box_aspect(1/box_ratio)
labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)


#  Mobility (moving average)
if batch == 2:
   ax[1][0].plot(time_interp2,meanM1900_cf, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM1900   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 400, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM1900FB , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM1825FB , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM1750FB , linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_Myr_run[2],moving_average_mobility_run[2], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[1][0].axvspan(1500.,2000., color='#000000', alpha=0.1,linewidth=0)
   ax[1][0].axvspan(3800.,3900., color='#000000', alpha=0.1,linewidth=0)
   ax[1][0].fill_between(time_interp2,minM1900_cf,maxM1900_cf,linewidth=0,color = palette[0],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM1900,maxM1900      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM1900FB,maxM1900FB  ,linewidth=0,color = palette[2],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM1825FB,maxM1825FB  ,linewidth=0,color = palette[3],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM1750FB,maxM1750FB  ,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 4:
   ax[1][0].plot(time_interp2,meanM1OM_DE   , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')    
   ax[1][0].plot(time_interp2,meanM2OM_DE   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM5OM_DE   , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanM10OM_DE  , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')         
   ax[1][0].plot(time_interp2,meanM10OM_5000, linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_Myr_run[6],moving_average_mobility_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[1][0].plot(time_Myr_run[12],moving_average_mobility_run[12], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])        
   ax[1][0].axvspan(1500.,2000., color='#000000', alpha=0.1,linewidth=0)
   ax[1][0].axvspan(3800.,3900., color='#000000', alpha=0.1,linewidth=0)         
   ax[1][0].fill_between(time_interp2,minM1OM_DE,maxM1OM_DE      ,linewidth=0,color = palette[0],alpha=0.2)    
   ax[1][0].fill_between(time_interp2,minM2OM_DE,maxM2OM_DE      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minM5OM_DE,maxM5OM_DE      ,linewidth=0,color = palette[2],alpha=0.2)    
   ax[1][0].fill_between(time_interp2,minM10OM_DE,maxM10OM_DE    ,linewidth=0,color = palette[3],alpha=0.2)   
   ax[1][0].fill_between(time_interp2,minM10OM_5000,maxM10OM_5000,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 6:
   ax[1][0].plot(time_interp2,meanMrho15       , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')  
   ax[1][0].plot(time_interp2,meanMrho20       , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')        
   ax[1][0].plot(time_interp2,meanMrho20alpha20, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanMrho20alpha25, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_Myr_run[5],moving_average_mobility_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[1][0].plot(time_Myr_run[6],moving_average_mobility_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[1][0].axvspan(1500.,2000., color='#000000', alpha=0.1,linewidth=0)
   ax[1][0].axvspan(3800.,3900., color='#000000', alpha=0.1,linewidth=0)        
   ax[1][0].fill_between(time_interp2,minMrho15,maxMrho15        ,linewidth=0,color = palette[0],alpha=0.2)    
   ax[1][0].fill_between(time_interp2,minMrho20,maxMrho20        ,linewidth=0,color = palette[1],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minMrho20alpha20,maxMrho20alpha20        ,linewidth=0,color = palette[2],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minMrho20alpha25,maxMrho20alpha25        ,linewidth=0,color = palette[3],alpha=0.2)    
elif batch == 8:
   ax[1][0].plot(time_Myr_run[6],moving_average_mobility_run[6], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[1][0].plot(time_Myr_run[7],moving_average_mobility_run[7], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[1][0].plot(time_interp2,meanMrho0alpha20     , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 500, mfc='#ffffff')
   ax[1][0].plot(time_interp2,meanMrho0alpha25     , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 500, mfc='#ffffff')        
   ax[1][0].plot(time_Myr_run[2],moving_average_mobility_run[2], linestyle=style[2], linewidth=1, alpha=1, color = palette[2])
   ax[1][0].plot(time_Myr_run[4],moving_average_mobility_run[4], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[1][0].axvspan(1500.,2000., color='#000000', alpha=0.1,linewidth=0)
   ax[1][0].axvspan(3800.,3900., color='#000000', alpha=0.1,linewidth=0)   
   ax[1][0].fill_between(time_interp2,minMrho0alpha20,maxMrho0alpha20          ,linewidth=0,color = palette[0],alpha=0.2)
   ax[1][0].fill_between(time_interp2,minMrho0alpha25,maxMrho0alpha25          ,linewidth=0,color = palette[1],alpha=0.2)  

ax[1][0].set_xlim(0., evolved_time)
ax[1][0].set_xscale('linear')
ax[1][0].set_xticks(time_ticks)  
ax[1][0].set_ylim(0., max_mobility)
ax[1][0].set_yscale('linear')
ax[1][0].set_yticks(mobility_ticks)
ax[1][0].set_ylabel('Mobility',labelpad=23)
ax[1][0].set_box_aspect(1/box_ratio)
labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)


# TTG mass
if batch == 2:
   ax[2][0].plot(time_interp3,meanTTG1900_cf, linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 4000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG1900   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 4000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG1900FB , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG1825FB , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG1750FB , linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_Myr_step_run[2],TotalTTGmass_run[2], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[2][0].fill_between(time_interp3,minTTG1900_cf,maxTTG1900_cf,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG1900,maxTTG1900      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG1900FB,maxTTG1900FB  ,linewidth=0,color = palette[2],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG1825FB,maxTTG1825FB  ,linewidth=0,color = palette[3],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG1750FB,maxTTG1750FB  ,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 4:
   ax[2][0].plot(time_interp3,meanTTG1OM_DE   , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 5000, mfc='#ffffff')    
   ax[2][0].plot(time_interp3,meanTTG2OM_DE   , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG5OM_DE   , linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTG10OM_DE  , linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 5000, mfc='#ffffff')          
   ax[2][0].plot(time_interp3,meanTTG10OM_5000, linestyle=style[4], linewidth=1, alpha=1, color = palette[4],marker=symbol[4], markersize=4, mec=palette[4], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_Myr_step_run[6],TotalTTGmass_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[2][0].plot(time_Myr_step_run[12],TotalTTGmass_run[12], linestyle=style[6], linewidth=1, alpha=1, color = palette[6])        
   ax[2][0].fill_between(time_interp3,minTTG1OM_DE,maxTTG1OM_DE      ,linewidth=0,color = palette[0],alpha=0.2)    
   ax[2][0].fill_between(time_interp3,minTTG2OM_DE,maxTTG2OM_DE      ,linewidth=0,color = palette[1],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG5OM_DE,maxTTG5OM_DE      ,linewidth=0,color = palette[2],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG10OM_DE,maxTTG10OM_DE    ,linewidth=0,color = palette[3],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTG10OM_5000,maxTTG10OM_5000,linewidth=0,color = palette[4],alpha=0.2)   
elif batch == 6:
   ax[2][0].plot(time_interp3,meanTTGrho15       , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 5000, mfc='#ffffff')  
   ax[2][0].plot(time_interp3,meanTTGrho20       , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 5000, mfc='#ffffff')        
   ax[2][0].plot(time_interp3,meanTTGrho20alpha20, linestyle=style[2], linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 5000, mfc='#ffffff')  
   ax[2][0].plot(time_interp3,meanTTGrho20alpha25, linestyle=style[3], linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 5000, mfc='#ffffff')  
   ax[2][0].plot(time_Myr_step_run[5],TotalTTGmass_run[5], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[2][0].plot(time_Myr_step_run[6],TotalTTGmass_run[6], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[2][0].fill_between(time_interp3,minTTGrho15,maxTTGrho15        ,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTGrho20,maxTTGrho20        ,linewidth=0,color = palette[1],alpha=0.2)        
   ax[2][0].fill_between(time_interp3,minTTGrho20alpha20,maxTTGrho20alpha20        ,linewidth=0,color = palette[2],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTGrho20alpha25,maxTTGrho20alpha25        ,linewidth=0,color = palette[3],alpha=0.2)   
elif batch == 8:
   ax[2][0].plot(time_Myr_step_run[6],TotalTTGmass_run[6], linestyle=style[4], linewidth=1, alpha=1, color = palette[4])
   ax[2][0].plot(time_Myr_step_run[7],TotalTTGmass_run[7], linestyle=style[5], linewidth=1, alpha=1, color = palette[5])
   ax[2][0].plot(time_interp3,meanTTGrho0alpha20     , linestyle=style[0], linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_interp3,meanTTGrho0alpha25     , linestyle=style[1], linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 5000, mfc='#ffffff')
   ax[2][0].plot(time_Myr_step_run[2],TotalTTGmass_run[2], linestyle=style[2], linewidth=1, alpha=1, color = palette[2])
   ax[2][0].plot(time_Myr_step_run[4],TotalTTGmass_run[4], linestyle=style[3], linewidth=1, alpha=1, color = palette[3])
   ax[2][0].fill_between(time_interp3,minTTGrho0alpha20,maxTTGrho0alpha20          ,linewidth=0,color = palette[0],alpha=0.2)
   ax[2][0].fill_between(time_interp3,minTTGrho0alpha25,maxTTGrho0alpha25          ,linewidth=0,color = palette[1],alpha=0.2) 
             
ax[2][0].set_xlim(0., evolved_time)    
ax[2][0].set_xscale('linear')
ax[2][0].set_xticks(time_ticks)
ax[2][0].set_ylim(0., 8.e22)
ax[2][0].set_yscale('linear')
ax[2][0].set_ylabel('TTG mass (kg)',labelpad=23)
ax[2][0].set_xlabel('Time (Myr)')
ax[2][0].set_box_aspect(1/box_ratio)    
labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)

ax[0][0].set_title(figure_title)

fig.text(0.01, 0.96, 'a', fontsize=15, fontweight="bold")
fig.text(0.01, 0.65, 'b', fontsize=15)
fig.text(0.01, 0.34, 'c', fontsize=15)

fig.tight_layout()
file_name='fig'+str('{:01d}'.format(batch))+'.pdf'
plt.savefig(here+'/../figures/'+file_name)


