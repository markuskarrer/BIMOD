'''
plots the fluxes and the absolute values of the number density and the mass density
'''

from IPython.core.debugger import Tracer ; debug=Tracer() #insert this line somewhere to debug

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import math
import os
import subprocess
from netCDF4 import Dataset
import traceback
import sys
#self-written functions
import __plotting_functions
import __postprocess_McSnow

#read variables passed by shell script
tstep = int(os.environ["tstep"])
tstep_end = int(os.environ["tstep_end"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MCexp"]
skipMC = False #(os.environ["skipMC"]=="True") #allows to run the scripts also if no McSnow data is there (only 1D-SB runs)
skipSB = True #(os.environ["skipSB"]=="True") #allows to hide plots for SB

if "McSnow_geom_specifier_onestring" in  os.environ.keys():
    McSnow_geom_specifier_onestring = os.environ["McSnow_geom_specifier_onestring"]
else:
    McSnow_geom_specifier_onestring = ''
if "separated_by_sensruns" in os.environ.keys():
    separated_by_sensruns = (os.environ["separated_by_sensruns"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
    separated_by_sensruns_onestring=""
if "separated_by_fallspeedsens" in os.environ.keys():
    separated_by_fallspeedsens = (os.environ["separated_by_fallspeedsens"]=="True") #plot a line for different sensitivity runs
else:
    print "separated_by_fallspeedsens not found in environ: set to False"
    separated_by_fallspeedsens = False
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''
#if "McSnow_geom_list" in os.environ.keys():
#    McSnow_geom_list_str = os.environ["McSnow_geom_list"]
#else:
#    McSnow_geom_list_str = ''
#directory of experiments
directory = MC_dir + "/experiments/"

number_of_plots = 5

figsize_height = 6.0/2.0*(number_of_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))

#get timestep
i_timestep=(tstep/10)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)

#plot different lines for different sensitivity runs
if separated_by_sensruns:
    sensrun_list = separated_by_sensruns_onestring

    sensrun_list = separated_by_sensruns_onestring.split('_')
else:
    sensrun_list = [""]
#plot different lines for different McSnow geometries
McSnow_geom_list = McSnow_geom_specifier_onestring.split('_')

axnum = plt.subplot2grid((number_of_plots, 1), (0, 0))
axmass = plt.subplot2grid((number_of_plots, 1), (1, 0))
axmean = plt.subplot2grid((number_of_plots, 1), (2, 0))
axnumflux = plt.subplot2grid((number_of_plots, 1), (3, 0))
axmassflux = plt.subplot2grid((number_of_plots, 1), (4, 0))

for i_sensrun, sensrun_now in enumerate(sensrun_list): #loop over different SB sensitivity runs
    for i_sensMC, sensrun_now_MC in  enumerate(McSnow_geom_list): #loop over different McSnow geometry sensitivity runs
        if separated_by_sensruns:#modify the experiment string
            #replace the geometry string
            experiment_splitted = experiment.split('_',2) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
            experiment_splitted[1] = sensrun_now #replace experiment part
            experiment = "_".join(experiment_splitted)
            

        if len(McSnow_geom_list)>1:#modify the experiment string to read the file with the write geometry
            experiment_splitted = experiment.split('_',3) #this results e.g. in ['1d', 'powerlawJplate', 'xi10000000_
            experiment_splitted[2] = sensrun_now_MC #replace experiment part
            experiment = "_".join(experiment_splitted)
            McSnow_geom_list_str = McSnow_geom_specifier_onestring
        else:
            McSnow_geom_list_str = str(McSnow_geom_list[0])
        #ATTENTION: quick fix
        if "stick1__" in experiment:
            experiment="__".join(experiment.split("stick1__")) #remove the string
        if "stick2__" in experiment:
            experiment="__".join(experiment.split("stick2__")) #remove the string
        if "stick3__" in experiment:
            experiment="__".join(experiment.split("stick3__")) #remove the string

            #debug()
        #read hei2massdens to get average /summed up values per height
        filestring_hei2massdens = directory + experiment + "/hei2massdens.dat"
        timestep = tstep/10 #TODO: do not hardcode the .. minute output interval here
        timestep_end = tstep_end/10 #TODO: do not hardcode the .. minute output interval here
        if skipMC:
            hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,timestep_end=timestep_end,empty_flag=True) #get empty arrays with the same dimension as the true arrays in order to not plot the McSnow data
        else:
            hei2massdens = __postprocess_McSnow.read_hei2massdens(filestring_hei2massdens,timestep=timestep,timestep_end=timestep_end)
            
        #######
        #read twomoment-data
        #######
        #define filestring
        #filestring = directory + experiment + "/twomom_d.dat"
        #load netCDF4 file
        twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
        #create dictionary for all variables from twomom_d.ncdf
        twomom = dict()

        #if necessary change name of variables
        varlist = twomom_file.variables
        #read 2mom variables to twomom dictionary
        for var in varlist:#read files and write it with different names in Data
            #do a skipSB similar to skipMC above
            if skipSB and var!="heights":
                #create empty variables so we dont see them in the plot
                twomom[var] = np.nan*np.ones_like(twomom_file.variables[var])
            else:
                twomom[var] = np.squeeze(twomom_file.variables[var])

        #from IPython.core.debugger import Tracer ; Tracer()()
        ####################################
        #plot mixing ratio + number density
        ####################################
        #number density
        mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

        
        if mass_num_flag==2:
            ax2 = axnum.twiny()
        else: #in case there is no need for a second axis, just pass the first ax twice
            ax2 = axnum
            
        axnum = __plotting_functions.plot_moments(axnum,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle=['-','--','-.',':'][i_sensrun],forced_markerMC=['','x','o','d'][i_sensMC])

        #mass density
        mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

        if mass_num_flag==2:
            ax2 = axmass.twiny()
        else: #in case there is no need for a second axis, just pass the first ax twice
            ax2 = axmass
            
        axmass = __plotting_functions.plot_moments(axmass,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle=['-','--','-.',':'][i_sensrun],forced_markerMC=['','x','o','d'][i_sensMC])

        ##
        #plot normalized mixing ration (q/qn) labelling normq...
        ##
        #calculate q/qn respectively Md/Nd

        ax2 = axmean
        
        axmean = __plotting_functions.plot_normmix(axmean,ax2,twomom,hei2massdens,i_timestep,forced_linestyle=['-','--','-.',':'][i_sensrun],forced_markerMC=['','x','o','d'][i_sensMC])


        ############
        #plot fluxes
        ############
        #number flux
        mass_num_flag = 0 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

        if mass_num_flag==2:
            ax2 = axnumflux.twiny()
        else: #in case there is no need for a second axis, just pass the first ax twice
            ax2 = axnumflux
            
        axnumflux = __plotting_functions.plot_fluxes(axnumflux,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle=['-','--','-.',':'][i_sensrun],forced_markerMC=['','x','o','d'][i_sensMC])

        #mass flux
        mass_num_flag = 1 #0-> plot only number flux; 1-> plot only mass flux; 2-> plot both 

        if mass_num_flag==2:
            ax2 = axmassflux.twiny()
        else: #in case there is no need for a second axis, just pass the first ax twice
            ax2 = axmassflux
            
        axmassflux = __plotting_functions.plot_fluxes(axmassflux,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=mass_num_flag,forced_linestyle=['-','--','-.',':'][i_sensrun],forced_markerMC=['','x','o','d'][i_sensMC])

if len(sensrun_list)>1:#add labels for the different sensruns
    for i_sensrun,sens_run in enumerate(sensrun_list):
        for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
            ax.plot(np.nan,np.nan,color='k',linestyle=['-','--','-.',':'][i_sensrun],label=sens_run)
            ax.legend()
if len(McSnow_geom_list)>1:#add labels for the different sensruns
    for i_sensMC,sensrun_now_MC in enumerate(McSnow_geom_list):
        for ax in [axnum,axmass,axmean,axnumflux,axmassflux]:
            ax.plot(np.nan,np.nan,color='k',linestyle='--',marker=['','x','o','d'][i_sensMC],label=['default','binary','N_mono \ndependent'][i_sensMC])
            ax.legend()

#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create directory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/fluxes_" + switch_off_processes_str + '_' + McSnow_geom_list_str + separated_by_sensruns_onestring + '_' + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
