'''
create an overview panel with properties of McSnow and PAMTRA output
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
import __postprocess_PAMTRA

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"]
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MCexp"]
adapt_version = int(os.environ["adapt_version"]) #reading the files of the appropriate adaption version
if "separated_by_categories" in os.environ.keys():
    separated_by_categories = (os.environ["separated_by_categories"]=="True") #plot a line for individual categories also
else:
    separated_by_categories = False 
if "separated_by_sensruns" in os.environ.keys():
    separated_by_sensruns = (os.environ["separated_by_sensruns"]=="True") #plot a line for different sensitivity runs
    separated_by_sensruns_onestring= os.environ["model_setup_specifier_onestring"]
else:
    print "separated_by_sensruns not found in environ: set to False"
    separated_by_sensruns = False
if "switch_off_processes" in os.environ.keys():
    switch_off_processes_str = os.environ["switch_off_processes"]
else:
    switch_off_processes_str = ''
#directory of experiments
directory = MC_dir + "/experiments/"

num_plots = 8 #number of subplots from pamtra variables and twomoment output
figsize_height = 6.0/2.0*(num_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))


if "semiidealized" in os.environ:
    height_bounds = __postprocess_McSnow.read_init_vals(MC_dir) #read the heights from the init_vals.txt file

else: #no height limits if it is not run by the McSnow_Pamtra_ICONinit script
    height_bounds = [3850,0] #set some default heigh-bounds
    
##############################
#plot McSnow pamtra output
##############################

#define file to read
print "adaptv:",adapt_version
if adapt_version==1:
    filename = "/adaptv1_t" + str(tstep) #+ ".nc"
elif adapt_version==2:
    filename = "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
elif adapt_version==3:
    filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
    
pam_filestring = directory + experiment + filename +".nc"

#read pamtra output to pamData dictionary
pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
axwaterfall = plt.subplot2grid((num_plots, 1), (4, 0), rowspan=4)
#plot spectrogram
axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=9.6,color='b',linestyle='--',z_lim=height_bounds[::-1])
axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=35.5,color='r',linestyle='--',z_lim=height_bounds[::-1])
axwaterfall = __plotting_functions.plot_waterfall(axwaterfall,pamData,freq=95.0,color='g',linestyle='--',z_lim=height_bounds[::-1])
#plot reflectivities
axrefl = plt.subplot2grid((num_plots, 1), (0, 0))

axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,pamData,linestyle='--')
#plot other radar-moments
axvDoppler = plt.subplot2grid((num_plots, 1), (1, 0))
axvDoppler = __plotting_functions.plot_pamtra_highermoments(axvDoppler,pamData,linestyle='--',moment='vDoppler')
axswidth = plt.subplot2grid((num_plots, 1), (2, 0))
axswidth = __plotting_functions.plot_pamtra_highermoments(axswidth,pamData,linestyle='--',moment='swidth')
axskewn = plt.subplot2grid((num_plots, 1), (3, 0))
axskewn = __plotting_functions.plot_pamtra_highermoments(axskewn,pamData,linestyle='--',moment='skewn')


#save figure
plt.tight_layout()
if not os.path.exists('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment)
out_filestring = "/spectra_" + switch_off_processes_str + "_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])

#plot spectrogram in new figure

#plot spectrograms
plt.clf() #clear figure


'''
for i_scheme,scheme in enumerate(["SB","MC"]):
    if scheme=="SB":
        #define file to read
        pam_filestring = directory + experiment + "/PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min.nc'
        #read pamtra output to pamData dictionary
        pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
    elif scheme=="MC":
        try:
            if adapt_version==1:
                filename = "/adaptv1_t" + str(tstep) #+ ".nc"
            elif adapt_version==2:
                filename = "/adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
            elif adapt_version==3:
                filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
                
            pam_filestring = directory + experiment + filename +".nc"
            pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)

        except:
            print "no spectral data found in:", pam_filestring, ". Is radar_mode set to simple?"
            continue
    #read pamtra output to pamData dictionary
    number_of_spectrograms = 3*2 #3 frequency and two schemes (McSnow and SB)
    for i_freq,freq in enumerate([9.6,35.5,95.0]):
        ax = plt.subplot2grid((number_of_spectrograms, 1), (i_scheme*3+i_freq, 0))
        ax = __plotting_functions.plot_pamtra_spectrogram(ax,pamData,freq=freq)
        ax.set_ylim([height_bounds[1],height_bounds[0]]) #this might be the range from the initialization height to the maximum of the masses
#plt.tight_layout()
out_filestring = "/spectrogram_" + testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/overview_panel/' + experiment + out_filestring + '.pdf'])
'''
