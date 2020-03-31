'''
this script creates ncdf data from McSnows mass2fr..dat file with the original and some postprocessed variables (like temporal averages)
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import numpy as np
import sys #for some debugging
import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
import os #import environment variables
from netCDF4 import Dataset 
import time

#self written #these files are linked here from the pythoncode/functions directory
import __postprocess_McSnow

#read variables passed by shell script
tstep = int(os.environ["tstep"]) #string with 4 numbers and 'min'
tstep_end = int(os.environ["tstep_end"]) #string with 4 numbers and 'min'
if "tstep_step" in os.environ.keys():
    tstep_step = int(os.environ["tstep_step"]) 
else:
    tstep_step=10
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_expdir = os.environ["MCexp"]

#directory of experiments
directory = MC_expdir + "/experiments/"

#read timestep from string
dtc = int(re.search(r'dtc(.*?)_nrp', experiment).group(1)) #this line gets the values after dtc and before _nrp -> timestep of collision in s

#perform temporal averaging
t_after_full = 0 #time after the timestep with the full output in seconds
#while t_after_full<av_tstep: #1. merge all timesteps <av_tstep in one tuple
for i_tstep,tstep_now in enumerate(range(tstep,tstep_end,tstep_step)):
    filestring_mass2fr = directory + experiment + "/mass2fr_" + str(tstep_now).zfill(4) + 'min_' + str(t_after_full).zfill(2) + 's' + ".dat"
    print 'reading: ' + filestring_mass2fr
    if i_tstep==0:
        SP_nonaveraged = (__postprocess_McSnow.read_mass2frdat(experiment,filestring_mass2fr),) #save first timestep to SP_nonaveraged
    else:
        SP_nonaveraged = SP_nonaveraged + (__postprocess_McSnow.read_mass2frdat(experiment,filestring_mass2fr),) #add SP of other timestep to SP_nonaveraged
        
#from IPython.core.debugger import Tracer ; Tracer()()    
if bool(SP_nonaveraged[0]): #skip the rest if no SP are there (f.e. when only analyzing the SB output) #the trick here is that empty dicts evaluate to False
    #2. compute the 'average' of the SP list, which results in a longer list (approx SP_pertimestep*averaged timesteps) but approx the same number of RP
    SP = __postprocess_McSnow.average_SPlists(SP_nonaveraged)
    #get number of timesteps which can be averaged averaged
    num_timesteps = len(SP_nonaveraged)

    #create netCDF4 file
    output_file=directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + '-' + str(tstep_end).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf'
    dataset = Dataset(output_file, 'w',  format='NETCDF4_CLASSIC')
    #ATTENTION: this has been previously:
    #dataset = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf', 'w',  format='NETCDF4_CLASSIC')

    #Global attributes
    dataset.description = 'list with SP and their properties'
    dataset.history = 'Created ' + time.ctime(time.time())

    #create Dimensions
    dim_SP_all = dataset.createDimension('dim_SP_all_av' +str((num_timesteps)*dtc),len(SP['xi']))

    #create Variables and fill them
    unit = ['kg','1','kg','1','kg m-3','kg','m','m','1','m2','m s-1','m2','kg','m']
    for i,key in enumerate(SP.keys()):
        dataset.createVariable(key,np.float32,('dim_SP_all_av' +str((num_timesteps)*dtc),)); dataset.variables[key].units = unit[i]
        dataset.variables[key][:] = SP[key]

    print output_file,dataset.history
    dataset.close()

#################################################
#also save SB-twomoment scheme output to netCDF4
#################################################

twomom_bool = True #TODO: do not hardcode this here?
if not twomom_bool:
    sys.exit(0)
#read number of vertical levels from string
nz = int(re.search(r'nz(.*?)_', experiment).group(1)) #this line gets the values after dtc and before _nrp -> timestep of collision in s

#read the twomom_d.dat file to the dictionary twomom
filestring_twomom_d = MC_expdir + '/experiments/' + experiment + '/twomom_d.dat' 
twomom = __postprocess_McSnow.read_twomom_d(filestring_twomom_d,nz) #returns the twomom variables in a dictionary

#create height axis from dz (layer thickness) array
heights = np.cumsum(twomom['dz'][0,::-1])[::-1] #fix index 0 is allowed unless level thickness change with height

#create netCDF4 file
dataset = Dataset(directory + experiment + '/twomom_d.ncdf', 'w',  format='NETCDF4_CLASSIC')
print directory + experiment + '/twomom_d.ncdf'

#Global attributes
dataset.description = 'Output of SB twomoment-scheme (embedded in McSnow) converted to netCDF4'
dataset.history = 'Created ' + time.ctime(time.time())

#create Dimensions
dim_height = dataset.createDimension('height',nz)
dim_time = dataset.createDimension('time',twomom[twomom.keys()[0]].shape[0])

#create Variables and fill them
#1. coordinate variables
dim_times = dataset.createVariable('times',np.float64,('time',))
dim_heights = dataset.createVariable('heights',np.int32,('height'))
timestep_out = 30 #TODO: do not hardcode this
dataset.variables['times'][:] = np.arange(0,twomom[twomom.keys()[0]].shape[0])*timestep_out; dataset.variables['times'].units = 'min'
dataset.variables['heights'][:] = heights; dataset.variables['heights'].units = 'm' 
#2. real variables
#unit = ['kg','1','kg','1','kg m-3','kg','m','m','1','m2','m s-1','m2','kg','m'] #TODO: define units
for i,key in enumerate(twomom.keys()):
    dataset.createVariable(key,np.float32,('time','height')); #dataset.variables[key].units = unit[i]
    dataset.variables[key][:] = twomom[key]

print "twomom_d file created at: ",dataset.history,directory + experiment + '/twomom_d.ncdf'
dataset.close()
