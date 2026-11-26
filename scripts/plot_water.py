# !/usr/bin/env python

# written on 9 Apr 2026 by Charitra Jain
# to plot water content for models presented in the manuscript
# Jain & Sobolev (2026): Controls on the tectonic evolution of the early Earth: 2. Contrastive effects of water
# submitted to JGR: Solid Earth

# Python packages required: os, numpy, matplotlib

# Once required packages are installed, execute the script from terminal with command "python plot_water.py"

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


### beginning of input ###
# user only has to input figure number (9 or 10) as batch
batch=9
### end of input ###


if batch == 9:
   model_list = ["Table5_models/01_QW1r1/", "Table5_models/04_QW2r1/", "Table5_models/07_QW5r1/", "Table5_models/10_QW10r1/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#A800FF']
   labels_a = ('1 OM (QW1r1)','2 OM (QW2r1)','5 OM (QW5r1)','10 OM (QW10r1)')
   figure_title = ('Water content in models with different size of initial external resevoir')
   location='upper center'
   style = ['solid','solid','solid','solid','solid']
   symbol = ['^','o','s','D','H']
   legend_columns = 4
   legend_panel = 1
elif batch == 10:
   model_list = ["Table5_models/07_QW5r1/", "Table6_models/06_Qaw20r3/", "Table6_models/07_Qaw20sat20r1/", "Table6_models/19_Qsat20H2Oetar3/"]
   palette = ['#C25B60','#FFBF70','#43B7C2','#A800FF']
   labels_a = ('100$\\%$ saturation (QW5r1)','a$_{w}$ = 2.0, 100$\\%$ saturation (Qaw20r3)','a$_{w}$ = 2.0, 20$\\%$ saturation (Qaw20sat20r1)','a$_{w}$ = $f$ (P,C), 20$\\%$ saturation (Qsat20H2Oetar3)')
   figure_title = ('Water content in models with initial external reservoir of 5 OM')
   location='upper center'
   style = ['solid','solid','solid','solid','solid']
   symbol = ['^','o','s','D','H']
   legend_columns = 2
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
mpl.rcParams['figure.figsize'] = [7,8] # (4 panels, width x height, 1=600px with dpi=300)
mpl.rcParams['xtick.direction'] = "in"
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.direction'] = "in"
mpl.rcParams["pdf.use14corefonts"] = True
box_ratio=4.5
evolved_time = 4500.
time_ticks=[0,500,1000,1500,2000,2500,3000,3500,4000,4500]
min_water=-1.
max_water=10.
water_ticks=[0,2.5,5,7.5,10]

time_step_Myr_run    =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
time_Myr_run    =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
water_ocean_run    =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
water_mantle_run     =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
water_surface_run     =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
water_UM_run     =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))
water_TZ_run     =  np.frompyfunc(list, 0, 1)(np.empty((total_runs), dtype=object))

run = int(0)
for path in model_list:

   print('Run # '+str(run+1))
   print('Model  '+path)

   time_step_Myr_list = []
   water_ocean_list = []
   water_mantle_list = []
   file = (here+'/../'+path+'water_ocean_mantle_total.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      time_step_Myr_list.append(float(p[0]))
      water_ocean_list.append(float(p[1]))
      water_mantle_list.append(float(p[2]))      

   time_step_Myr_run[run] = np.array(time_step_Myr_list)
   water_ocean_run[run] = np.array(water_ocean_list)
   water_mantle_run[run] = np.array(water_mantle_list)
   time_step_Myr_list.clear()
   water_ocean_list.clear()
   water_mantle_list.clear()
   print('Run #'+str(run+1)+': Final ocean water (OM):   '+str(water_ocean_run[run][len(water_ocean_run[run])-1]))
   print('Run #'+str(run+1)+': Final mantle water (OM):  '+str(water_mantle_run[run][len(water_mantle_run[run])-1]))

   time_Myr_list = []
   water_surface_list = []
   file = (here+'/../'+path+'water_amount_OM_1km_10km.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      time_Myr_list.append(float(p[0]))
      water_surface_list.append(float(p[1]))
   
   time_Myr_run[run] = np.array(time_Myr_list)
   water_surface_run[run] = np.array(water_surface_list)
   time_Myr_list.clear()
   water_surface_list.clear()
   print('Run #'+str(run+1)+': Final surface water (OM): '+str(water_surface_run[run][len(water_surface_run[run])-1]))


   water_UM_list = []
   file = (here+'/../'+path+'water_amount_OM_15km_400km.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      water_UM_list.append(float(p[1]))
   
   water_UM_run[run] = np.array(water_UM_list)
   water_UM_list.clear()
   print('Run #'+str(run+1)+': Final UM water (OM):      '+str(water_UM_run[run][len(water_UM_run[run])-1]))
   
   
   water_TZ_list = []
   file = (here+'/../'+path+'water_amount_OM_410km_670km.dat')
   f = open(file,'r')
   lines = f.readlines()
   f.close
   for line in lines:
      p = line.split()
      water_TZ_list.append(float(p[1]))
   
   water_TZ_run[run] = np.array(water_TZ_list)
   water_UM_list.clear()
   print('Run #'+str(run+1)+': Final TZ water (OM):      '+str(water_TZ_run[run][len(water_TZ_run[run])-1]))

   run = run + 1


fig, ax = plt.subplots(4, 1, squeeze=False)
ax[0][0].plot(time_step_Myr_run[0],water_ocean_run[0], linestyle=style[0],    linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 20000, mfc='#ffffff')
ax[0][0].plot(time_step_Myr_run[1],water_ocean_run[1], linestyle=style[1],    linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 20000, mfc='#ffffff')
ax[0][0].plot(time_step_Myr_run[2],water_ocean_run[2], linestyle=style[2],    linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 20000, mfc='#ffffff')
ax[0][0].plot(time_step_Myr_run[3],water_ocean_run[3], linestyle=style[3],    linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 20000, mfc='#ffffff')

ax[0][0].set_xlim(0., evolved_time)
ax[0][0].set_xscale('linear')
ax[0][0].set_xticks(time_ticks)
ax[0][0].set_ylim(1e-2,10)
ax[0][0].set_yscale('log')
ax[0][0].set_ylabel('External reservoir (OM)')
ax[0][0].set_box_aspect(1/box_ratio)


ax[1][0].plot(time_Myr_run[0],water_surface_run[0], linestyle=style[0],    linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[1][0].plot(time_Myr_run[1],water_surface_run[1], linestyle=style[1],    linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[1][0].plot(time_Myr_run[2],water_surface_run[2], linestyle=style[2],    linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
ax[1][0].plot(time_Myr_run[3],water_surface_run[3], linestyle=style[3],    linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 200, mfc='#ffffff')

ax[1][0].set_xlim(0., evolved_time)
ax[1][0].set_xscale('linear')
ax[1][0].set_xticks(time_ticks)
ax[1][0].set_ylim(0.,1.5)
ax[1][0].set_yscale('linear')
ax[1][0].set_yticks([0,0.5,1,1.5])
ax[1][0].set_ylabel('Surface layer (OM)',labelpad=10)
ax[1][0].set_box_aspect(1/box_ratio)

ax[2][0].plot(time_Myr_run[0],water_UM_run[0], linestyle=style[0],    linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[2][0].plot(time_Myr_run[1],water_UM_run[1], linestyle=style[1],    linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[2][0].plot(time_Myr_run[2],water_UM_run[2], linestyle=style[2],    linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
ax[2][0].plot(time_Myr_run[3],water_UM_run[3], linestyle=style[3],    linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 200, mfc='#ffffff')

ax[2][0].set_xlim(0., evolved_time)
ax[2][0].set_xscale('linear')
ax[2][0].set_xticks(time_ticks)
ax[2][0].set_ylim(0.,5.)
ax[2][0].set_yscale('linear')
ax[2][0].set_yticks([0,1,2,3,4,5])
ax[2][0].set_ylabel('Upper mantle (OM)',labelpad=18)
ax[2][0].set_box_aspect(1/box_ratio)


ax[3][0].plot(time_Myr_run[0],water_TZ_run[0], linestyle=style[0],    linewidth=1, alpha=1, color = palette[0],marker=symbol[0], markersize=4, mec=palette[0], mew='0.5', markevery = 200, mfc='#ffffff')
ax[3][0].plot(time_Myr_run[1],water_TZ_run[1], linestyle=style[1],    linewidth=1, alpha=1, color = palette[1],marker=symbol[1], markersize=4, mec=palette[1], mew='0.5', markevery = 200, mfc='#ffffff')
ax[3][0].plot(time_Myr_run[2],water_TZ_run[2], linestyle=style[2],    linewidth=1, alpha=1, color = palette[2],marker=symbol[2], markersize=4, mec=palette[2], mew='0.5', markevery = 200, mfc='#ffffff')
ax[3][0].plot(time_Myr_run[3],water_TZ_run[3], linestyle=style[3],    linewidth=1, alpha=1, color = palette[3],marker=symbol[3], markersize=4, mec=palette[3], mew='0.5', markevery = 200, mfc='#ffffff')

ax[3][0].set_xlim(0., evolved_time)
ax[3][0].set_xscale('linear')
ax[3][0].set_xticks(time_ticks)
ax[3][0].set_xlabel('Time (Myr)')
if batch == 9:
   ax[3][0].set_ylim(0.,10.)
   ax[3][0].set_yscale('linear')
   ax[3][0].set_yticks([0,2,4,6,8,10])
   ax[3][0].set_ylabel('Transition zone (OM)',labelpad=13)
else:
   ax[3][0].set_ylim(0.,5.)
   ax[3][0].set_yscale('linear')
   ax[3][0].set_yticks([0,1,2,3,4,5])
   ax[3][0].set_ylabel('Transition zone (OM)',labelpad=19)

ax[3][0].set_box_aspect(1/box_ratio)

labels = labels_a
ax[legend_panel][0].legend(ncol = legend_columns, labels = (labels), loc=location, labelspacing = 0.5, facecolor = 'none', edgecolor = 'none', fancybox = False)

ax[0][0].set_title(figure_title)

fig.text(0.01, 0.97, 'a', fontsize=15)
fig.text(0.01, 0.73, 'b', fontsize=15)
fig.text(0.01, 0.50, 'c', fontsize=15)
fig.text(0.01, 0.26, 'd', fontsize=15)


fig.tight_layout()
file_name='fig'+str('{:01d}'.format(batch))+'.pdf'
plt.savefig(here+'/../figures/'+file_name)

