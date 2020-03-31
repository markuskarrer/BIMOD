import numpy as np


def read_pamtra(pam_filestring):
    '''
    read and plot pamtra output
    INPUT: pam_filestring: full string of the path of the pamtra_output
    '''
    #import modules
    from netCDF4 import Dataset #for reading NetCDF

    #load file
    pam_file = Dataset(pam_filestring,mode='r')
    #create dictionary for all variables from PAMTRA
    pamData = dict()
    #print pam_file.variables.keys() 

    #if necessary change name of variables
    varlistin = pam_file.variables
    varlistout = pam_file.variables
    #read PAMTRA variables to pamData dictionary
    for varin,varout in zip(varlistin,varlistout):#read files and write it with different names in Data
        pamData[varout] = np.squeeze(pam_file.variables[varin])
        
    return pamData