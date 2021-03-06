'''
create an overview panel with properties of McSnow and PAMTRA output
'''

from IPython.core.debugger import Tracer ; debug=Tracer() #insert this line somewhere to debug

import numpy as np
import xarray as xr
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
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]
#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #ATTENTION: this is overwritten
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MCexp"]

#directory of experiments
directory = MC_dir + "/experiments/"

num_plots = 3 #number of subplots from pamtra variables and twomoment output
figsize_height = 6.0/2.0*(num_plots)
fig	=	plt.figure(figsize=(8.0,figsize_height))#figsize=(4, 4))


##############################
#plot McSnow pamtra output
##############################

#define file to read
filename = "/adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'

pam_filestring = directory + experiment + filename +".nc"

#read pamtra output to pamData dictionary
pamData = __postprocess_PAMTRA.read_pamtra(pam_filestring)
#plot reflectivities
axrefl = plt.subplot2grid((num_plots, 1), (0, 0))
axrefl = __plotting_functions.plot_pamtra_Ze(axrefl,pamData,linestyle='--')
#plot DWRs 
axDWR = plt.subplot2grid((num_plots, 1), (1, 0))
axDWR.plot(pamData["Ze"][:,0]-pamData["Ze"][:,1],pamData["height"],label="DWR_X_Ka") 

axDWR.plot(pamData["Ze"][:,1]-pamData["Ze"][:,2],pamData["height"],label="DWR_Ka_W") 
axDWR.legend()
axDWR.set_xlabel("DWR [dB]")
axDWR.set_ylabel("height [m]")


axSpec = plt.subplot2grid((num_plots, 1), (2, 0))
i_h_plot,h_plot = find_nearest(pamData["height"],3000)
#debug()
#i_eight = 0 #around 5km-> -15degC #ATTENTION: change profile
i_freq = 1 #Ka-Band
#plot the spertrum
axSpec.plot(-pamData["Radar_Velocity"][i_freq,:],pamData["Radar_Spectrum"][i_h_plot,i_freq,:],label="forw op")
axSpec.text(1,-45,"height=" + str(h_plot))
#add observed spectrum
selObsSpec = xr.open_dataset('/home/jdias/Projects/radarCode/tripexPolPro/selSpecDT.nc')
axSpec.plot(selObsSpec.kaSpec.vel,selObsSpec.kaSpec,label="obs Ka")


axSpec.set_xlabel("DV [m/s]")
axSpec.set_ylabel("power")
axSpec.set_xlim([-2,2])
axSpec.set_ylim([-50,-5])
axSpec.legend()
plt.grid(which="both")


#save figure
plt.tight_layout()
'''
if not os.path.exists('/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment): #create direktory if it does not exists
    os.makedirs('/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment)
out_filestring = "/" +testcase + "_av_" + str(av_tstep) + "_t" + str(tstep).zfill(4) + 'min'
plt.savefig('/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment + out_filestring + '.pdf', dpi=400)
plt.savefig('/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment + out_filestring + '.png', dpi=400)
print 'The pdf is at: ' + '/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment + out_filestring + '.pdf'
subprocess.Popen(['evince','/home/mkarrer/Dokumente/plots/4Jose_bimodality/' + experiment + out_filestring + '.pdf'])
'''
