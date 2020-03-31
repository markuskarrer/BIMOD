#import modulus
import numpy as np
import sys

from IPython.core.debugger import Tracer ; debug=Tracer()


# define functions to postprocess McSnow

def read_mass2frncdf(experiment,filestring):
    '''
    read SP-properties from the mass2fr....ncdf file into the SP-dictionary
    INPUT experiment: descriptor (also folder name) of the experiment
    '''
    from netCDF4 import Dataset
    #load file
    SP_file = Dataset(experiment+ '/' + filestring,mode='r')
    #create dictionary for all superparticle variables
    SP = dict()

    #if necessary change name of variables
    varlist = SP_file.variables
    #read PAMTRA variables to pamData dictionary
    for var in varlist:#read files and write it with different names in Data
        SP[var] = np.squeeze(SP_file.variables[var])


    return SP

    
def read_mass2frdat(experiment,filestring):
    '''
    read SP-properties from the mass2fr.dat file into the SP-dictionary
    INPUT experiment: descriptor (also folder name) of the experiment
    '''
    
    #load file from .dat
    #SP_fullinfo = np.loadtxt(experiment + '/'+ filestring)
    SP_fullinfo = np.loadtxt(filestring)
    
    #create dictionary
    SP = dict()
    #names of variables in mass2fr
    varnames =             ["m_tot","Frim","height","d_rime","vt","xi",    "rhor","a_rime","mr_crit","diam",    "proj_A",   "mm",            "m_i",   "m_wat"]
    #                       total   rime             ??    fall  multi-    rime     ?     crit mass  diameter  projected   monomer          mass of   mass of 
    #                       mass    fraction         ??    speed plicity  density   ?  compl. infilling          area     multiplicity       ice       water
    if SP_fullinfo.shape[0]==0 or len(SP_fullinfo.shape)==1:
        #print "error: check if", filestring, "is empty or just one SP is there (in postprocess_McSNow read_mass2frdat() )"
        #print "if you want to continue anyway (f.e. for runs with only the SB-scheme) just press enter"
        #raw_input()
        pass
    else:
        if not SP_fullinfo[0,:].shape[0]==len(varnames):
            print "error: compare number of variables in varnames (", len(varnames), ") with variables in MCsnows write_icefraction() (", not SP_fullinfo[0,:].shape[0], ") from mo_output.f90"
            sys.exit(1)
        #fill SP-dictionary with different variables
        for i,key in enumerate(varnames):
            SP[key] = SP_fullinfo[:,i]

    return SP

def read_hei2massdens(filestring,timestep=0,timestep_end=0,empty_flag=False):
    '''
    read height-profiles from the mass2fr.dat file into the SP-dictionary
    INPUT   experiment: descriptor (also folder name) of the experiment
            timestep: integer defining which output step to choose
            timestep_end: integer defining which output step to choose as the end of average (we just take all available)
            empty_flag: return arrays with zeros (but the same shape); this can be done to not plot McSnow data
    '''
    if timestep_end==0:
        timestep_end=timestep #set end timestep to start if not defined
    #load file from .dat
    heightprofiles_fullinfo = np.loadtxt(filestring)
    #create dictionary
    hei2massdens = dict()
    #names of variables in hei2massdens
    varnames = ["z", "Ns",
                "Nd",           "Md",          "Fn",           "Fm",            "Fmono", 
                "Nd_mm1",   "Md_mm1",           "Fn_mm1",   "Fm_mm1",       "Fmono_mm1", 
                "Nd_unr",   "Md_unr",           "Fn_unr",   "Fm_unr",       "Fmono_unr", 
                "Nd_grp",   "Md_grp",           "Fn_grp",   "Fm_grp",       "Fmono_grp", 
                "Nd_liq",   "Md_liq",           "Fn_liq",   "Fm_liq",       "Fmono_liq"]
    #        height  number of SP 
    #         number density    mass density    number flux   mass flux     monomer flux
    #   same values as above but only for: pristine
    #                                      unrimed
    #                                      graupel
    #                                      liquid

    if not heightprofiles_fullinfo[0,:].shape[0]==len(varnames):
        print "error: compare number of variables in varnames (",len(varnames),") with variables in MCsnows write_hprofiles() (", heightprofiles_fullinfo[0,:].shape[0],")from mo_output.f90"
        sys.exit(1)
    #fill SP-dictionary with different variables
    #read z-array first to crop hei2massdens array to the relevant timesteps
    z_for_crop = heightprofiles_fullinfo[:,0]
    index_start_of_timesteps = np.where(z_for_crop==0)[0] #this is an array with indices which indicate at which line a new timestep starts
    try:
        N_heights = index_start_of_timesteps[1] - index_start_of_timesteps[0]
    except:
        print "error in __postprocess_McSnow.read_hei2massdens: is the file: " + filestring + " not complete?"
        

    i_start = index_start_of_timesteps[timestep]; 

    #going to the beginning of next timestep and then one line back fails if the last timestep is read in -> treat this case seperately
    if (index_start_of_timesteps.shape[0]-1)==timestep_end: #means: last timestep
        i_end = heightprofiles_fullinfo.shape[0]
    else: #not last timestep
        if timestep_end==0: #means: no average interval selected
            i_end = index_start_of_timesteps[timestep+1]-1
        else:
            i_end = index_start_of_timesteps[timestep_end+1]
    for i,key in enumerate(varnames):
        if empty_flag:
            hei2massdens[key] = np.zeros_like(heightprofiles_fullinfo[i_start:i_end,i])
        else:
            if key=="z": #z is structured differently
                tmp_var_reshaped =  heightprofiles_fullinfo[i_start:i_end,i].reshape((timestep_end-timestep+1,N_heights)) #extract array with (n_heights,n_timesteps_to_average=timestep_end-timestep)
                hei2massdens[key] = np.mean(tmp_var_reshaped,axis=0) #average over all timesteps
            else:
                tmp_var_reshaped =  heightprofiles_fullinfo[i_start:i_end,i].reshape((timestep_end-timestep+1,N_heights)) #extract array with (n_heights,n_timesteps_to_average=timestep_end-timestep)
                hei2massdens[key] = np.mean(tmp_var_reshaped,axis=0) #average over all timesteps

    #calculate unrimed properties as residuum
    for prop in ["Nd","Md","Fn","Fm","Fmono"]: #loop over all properties
        prop_tmp = hei2massdens[prop][:] #save total
        for categ in ["mm1","unr","grp","liq"]:
            prop_tmp = prop_tmp - hei2massdens[prop + '_' + categ][:] #subtract quantity of each category
        hei2massdens[prop + "_rimed"] = prop_tmp[:] #save residuum for _rimed of this property
    for prop in ["Nd","Md","Fn","Fm","Fmono"]: #loop over all properties
        hei2massdens[prop + "_rimed"] = np.where(hei2massdens["Md_rimed"]>2e-6,hei2massdens[prop + "_rimed"],0) #remove noisy height-levels

    return hei2massdens

def average_SPlists(SP_nonaveraged):
    '''
    average the SP output file 'massfr2_....dat' over some timesteps
    INPUT: SP_nonaveraged: tuple containing dictionaries with the SP-list
    '''
    #get number of timesteps which are averaged
    num_timesteps = len(SP_nonaveraged)
    SP_averaged = dict()

    for key in SP_nonaveraged[0].keys():
        SP_averaged[key]=[]
        for i in range(0,num_timesteps):
            SP_averaged[key] = np.concatenate((SP_averaged[key],SP_nonaveraged[i][key]))
        
    #divide multiplicity of merged SP-dictionary by the number of timesteps in order to receive snapshot and not a sum of the superparticle over the average time
    roundbool = False #true-> round to integer; false -> do not round
    if roundbool:
        total_xi_before_round = np.sum(SP_averaged['xi'])
        SP_averaged['xi'] = np.rint(SP_averaged['xi']/num_timesteps) #round to integer
        total_xi_after_round =  np.sum(SP_averaged['xi'])*num_timesteps
        print 'error made by rounding after merging of timestep:'
        print 'total_xi_before: ',total_xi_before_round ,'total_xi_after-total_xi_before: ',total_xi_after_round-total_xi_before_round,' ratio: ',(total_xi_after_round-total_xi_before_round)/total_xi_before_round
    else: # and 'xi' in SP_averaged.keys:
        SP_averaged['xi'] = SP_averaged['xi']/num_timesteps #non-integer multiplicity is allowed
    #elif roundbool and (not 'xi' in SP_averaged.keys):
    #    pass #this is a workaround for just SB runs


    return SP_averaged

def read_twomom_d(filestring,nz):
    '''
    read output from twomoment-scheme (embedded in McSnow) from the twomom_d.dat file into the twomom-dictionary
    INPUT   experiment: descriptor (also folder name) of the experiment #ATTENTION: removed, because it was not used
            filestring: full path of file which has the twomom_d data
            nz: number of vertical levels in output file
    '''
    
    #load file from .dat
    twomom_fullinfo = np.loadtxt(filestring)
    
    #create dictionary
    twomom = dict()
    #names of variables in mass2fr
    varnames_twomom = ["dz",        "rho",    "pres",          "qv",              "qc","qnc",                   "qr","qnr", "qi","qni", "qs","qns", "qg","qng", "qh","qnh" ,"ninact", "t",      "w" ,"fr","fnr","fi","fni","fs","fns","fg","fng","fh","fnh"]
    #           (Schichtdichte, Gesamtdichte, Druck, Wasserdampfmassendichte, Wolkenwasser massen und anzahldichte, regen, Eis,       Schnee, Graupel ,    Hagel ,            IN, Temperatur, Vertikale Geschwindigkeit w, fluxes of moments)
    if not twomom_fullinfo.shape[1]/len(varnames_twomom)==nz:
        print "error in postprocess_McSnow.read_twomom_d: compare number of variables in varnames_twomom(", len(varnames_twomom) ,") with variables in MCsnows t_2mom_indices in mo_2mom_mcrph_driver.f90(", twomom_fullinfo.shape[1] ,"); the ratio must be a multiple of nz "
        sys.exit(1)
    
    #fill twomom-dictionary with different variables
    for i,key in enumerate(varnames_twomom):
        twomom[key] = np.zeros([twomom_fullinfo.shape[0],nz])
        for timestep in range(0,twomom_fullinfo.shape[0]):
            twomom[key][timestep,:] = twomom_fullinfo[timestep,(i*nz):(i+1)*nz] #dimension of each key is now [time,height]
    '''
    print "qi",twomom["qi"][timestep,:]
    print "qni",twomom["qni"][timestep,:]
    print "fi",twomom["fi"][timestep,:]
    print "fni",twomom["fni"][timestep,:]
    raw_input()
    '''
            
    return twomom
    
def read_atmo(experiment,filestring_atmo):
    '''
    read atmospheric variables from McSnow experiment directory
    INPUT:  experiment: specifier of experiment
    '''

    #load file from .dat
    atmo_fullinfo = np.loadtxt(filestring_atmo)
    
    #create dictionary
    atmo = dict()
    #names of variables in atmo
    varnames_atmo = ["z", "rho", "T", "p", "eta", "ssat", "rh", "psatw", "psati", "lwc"]
    
    if not atmo_fullinfo[0,:].shape[0]==len(varnames_atmo):
        print "error: compare number of variables in varnames_atmo with variables in MCsnows plot_atmo() from mo_check.f90"
        sys.exit(1)
    for i,key in enumerate(varnames_atmo):
        atmo[key] = atmo_fullinfo[:,i]
    
    return atmo

def interpolate_atmo(atmo,heightvec):
    '''
    interpolate atmospheric variables given in arbitrary vertical levels to heightvec
    INPUT:  atmo: dictionary with atmospheric variables
            heightvec: vector on which the interpolation should be fitted
    '''
    
    save_atmoz=atmo["z"] #save atmo["z"] so that it can also be interpolated as a key in atmo
    reverse=1
    if save_atmoz[1]<save_atmoz[0]: #must be increasing: if not reverse them and all variables
        save_atmoz=save_atmoz[::-1]
        reverse=-1 #save_atmoz=save_atmoz[::-1]
    
    for i,key in enumerate(atmo.keys()):
        atmo[key] = np.interp(heightvec,save_atmoz,atmo[key][::reverse])
        
    return atmo

def interpolate_2height(dictio,target_heightvec,curr_heightvec):
    '''
    interpolate atmospheric variables given in arbitrary vertical levels to heightvec
    INPUT:  atmo: dictionary with atmospheric variables
            heightvec: vector on which the interpolation should be fitted
    '''
    
    if curr_heightvec[1]<curr_heightvec[0]: #must be increasing: if not reverse them and all variables
        curr_heightvec=curr_heightvec[::-1]
        reverse=-1 #save_atmoz=save_atmoz[::-1]
    else:
        reverse=1
        
    for i,key in enumerate(dictio.keys()):
        try:
            dictio[key] = np.interp(target_heightvec,curr_heightvec,dictio[key][::reverse])
        except:
            print 'variable: (' + key + ') could not be interpolated to heightvec'
    
        
    return dictio

def read_MCdistribution(filestring,n_tsteps):
    '''
    read the distribution of the McSnow run
    ARGUMENTS:
        filestring: name of file to read
        n_tsteps: number of timesteps in file
    '''

    #load file from .dat
    distributions_fullinfo = np.loadtxt(filestring) #distributions_fullinfo.shape -> (...,26) (200radii*timesteps,26variables)
    distributions_fullinfo_reshaped = np.reshape(distributions_fullinfo,(n_tsteps,200,26))
    
    #distributions_fullinfo_reshaped = [1,:,2] #timestep 1 var 2 all sizes
    distribution_dic=dict()

    #names of variables in hei2massdens
    varnames = ["radius", 
                "Nd",           "Md",          "Fn",           "Fm",            "Fmono", 
                "Nd_liq",   "Md_liq",           "Fn_liq",   "Fm_liq",       "Fmono_liq",
                "Nd_mm1",   "Md_mm1",           "Fn_mm1",   "Fm_mm1",       "Fmono_mm1", 
                "Nd_unr",   "Md_unr",           "Fn_unr",   "Fm_unr",       "Fmono_unr", 
                "Nd_grp",   "Md_grp",           "Fn_grp",   "Fm_grp",       "Fmono_grp"] 
    for i_key,key in enumerate(varnames):
        if key!="radius":
            distribution_dic[key] = distributions_fullinfo_reshaped[:,:,i_key] #dimensions (diam,timestep)
        else:
            distribution_dic["radius"] = distributions_fullinfo[0:200,0]
     
        #print key,distribution_dic[key][:,-1]
    return distribution_dic

def separate_by_height_and_diam(SP,nbins=100,nheights=51,model_top=500,diamrange=[-9,0],calconly="None"):
    '''
    count number of real particles (RP) in height-diameter bins
    INPUT:  SP:dictionary with properties of all SP of one timestep 
            nbins: number of bins in diameter array
            nheights: number of heights in height array
            model_top: top of model (also top of height area) [m]
            diamrange: specifies range of diameter array from 10**firstvalue to 10**secondvalue
            calconly: 'None': calc all height bins, else calculate only height bins which starts with [<calconly-array>]
    '''

    #create height vector
    zres = model_top/(nheights-1) #vertical resolution
    heightvec = np.linspace(0,model_top,nheights) #start with 0+zres and go n_heigts step up to model_top
    #heightvec[-1] = np.inf #be sure to have all particles in heighest bin
    #define arrays with sp_diam
    d_bound_ds = np.logspace(diamrange[0],diamrange[1],nbins+1) #array from 1nm to 1m
    d_ds = d_bound_ds[:-1] + 0.5*np.diff(d_bound_ds) #diameter at center of bins
    #create dictionary for binned values and fill it with zeros
    binned_val = dict()
    binned_val["d_counts"] = np.zeros([nheights-1,nbins]) #number of RP at each h-D bin
    binned_val["d_counts_no_mult"] = np.zeros([nheights-1,nbins]) #number of SP at each h-D bin
    binned_val["av_Frim"] = np.zeros([nheights-1,nbins]) #mass averaged rime fraction for each h-D bin
    binned_val["av_rhor"] = np.zeros([nheights-1,nbins]) #mass averaged rime density for each h-D bin
    binned_val["av_mm"] = np.zeros([nheights-1,nbins]) #multiplicity weighted monomer number
    binned_val["nav_vt"] = np.zeros([nheights-1,nbins]) #multiplicity weighted fall speed
    binned_val["mass_in_bin"] = np.zeros([nheights-1,nbins]) #total mass within bin
    
    #get the number of particles (SP*sp_multiplicity) at each bin (binned by height and sp_diameter)
    for i in range(0,nheights-1):
        if isinstance(calconly,basestring) or (heightvec[i] in calconly): #skip heights which are not needed for plotting
            for j in range(0,nbins):
                condition_in_bin = np.logical_and(
                                np.logical_and(d_bound_ds[j]<=SP["diam"],SP["diam"]<d_bound_ds[j+1]),
                                np.logical_and(heightvec[i]<=SP["height"],SP["height"]<heightvec[i+1]),
                                )
                
                binned_val["d_counts"][i,j] = np.sum(np.where(condition_in_bin,SP["xi"],0))
                binned_val["d_counts_no_mult"][i,j] = np.sum(np.where(condition_in_bin,1,0))
                binned_val["mass_in_bin"][i,j] = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["xi"],0))
                #get total number of RP per bin
                multipl_bin = np.sum(np.where(condition_in_bin,SP["xi"],0))            
                #get sum of qirim per bin
                qirim_bin = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["Frim"]*SP["xi"],0))
                #get sum of qitot
                qitot_bin = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["xi"],0))
                binned_val["av_Frim"][i,j] = qirim_bin/qitot_bin #calc. mass averaged rime fraction
                #calc. rime mass averaged rime density
                binned_val["av_rhor"][i,j] = np.sum(np.where(condition_in_bin,SP["m_tot"]*SP["Frim"]*SP["xi"]*SP["rhor"],0))/qitot_bin
                #calc. multipl averaged
                binned_val["av_mm"][i,j] = np.sum(np.where(condition_in_bin,SP["mm"]*SP["xi"],0))/multipl_bin
                #calc number averaged fall speed
                binned_val["nav_vt"][i,j] = np.sum(np.where(condition_in_bin,-SP["vt"]*SP["xi"],0))/multipl_bin
                #print i,j,binned_val["mass_in_bin"][i,j]
        else:
            pass

    binned_val["RPpSP"] = binned_val["d_counts"]/binned_val["d_counts_no_mult"]
    return binned_val,heightvec,d_bound_ds,d_ds,zres


def calculate_Vbox(experiment,zres):
    '''
    calculate the volume of the self (considering h) and the McSnow defined (considering A) box volume
    INPUT:  experiment: experiment-descriptor which includes information on runscripts
            delta_z: height steps in self-defined z-vector
    '''
    import re #for searching and selecting in strings
    #read box area from string
    box_area = float(re.search('ba(.*)', experiment).group(1)) *1./100. #this line gets the values after ba and before /mass2fr -> should always be the last specification in the experiment-string; 1/100m^2 is the specification of the unit in the McSnow runscripts
    Vbox=box_area*zres
    return Vbox

def conv_num2numpm3(var,Vbox):
    #convert number in height bin [#] to numper per volume [#/m3]
    return var/Vbox
    
def return_parameter_mD_AD_rel(experiment):
    '''
    select the m-D and A-D relationship based on the beginnig of the experiment name
    '''
    import re
    rhoi = 919.#from McSnows mo_atmo_types.f90
    rhol = 1e3 #from McSnows mo_atmo_types.f90
    #read string which specifies the m-D and A-D relationship
    #try:
    ADmD_spec_string = re.search('1d_(.*)_xi', experiment).group(1)
    #except: #this was the standard (deviations from the standard configurations should have been between 1d_ and _xi
    #    try:
    #        ADmD_spec_string = re.search('1d_(.*)_xi', experiment).group(1)
    #    except: #ultima ratio: take the default one
    #        ADmD_spec_string=''
    if True:
        unr_bet = 2.1 ; unr_alf = 0.0028*10.**(2.*unr_bet-3.) #from McSnows mo_mass2diam.f90 #m-D
        unr_gam = 1.88 ; unr_sig = 0.2285*10**(2.*unr_gam-4.) #A-D
        print "\n\n\n\n MAKE parameter in McSnow consistent with SB setup ( in __postprocess_McSnow)"
    else:
        if ADmD_spec_string=='':
            unr_bet = 2.1 ; unr_alf = 0.0028*10.**(2.*unr_bet-3.) #from McSnows mo_mass2diam.f90 #m-D
            unr_gam = 1.88 ; unr_sig = 0.2285*10**(2.*unr_gam-4.) #A-D
        elif ADmD_spec_string=='mDADJaggdent':
            #   Jussi's aggregates of dendrites
            unr_alf = 0.01243; unr_bet = 2.00000 #found commented in McSnows mo_mass_diam.f90 #m-D
            unr_sig = 0.05625; unr_gam = 1.81000 #A-D
        elif ADmD_spec_string=='SBice':
            unr_bet = 1.0/0.39; unr_alf = (1.0/0.835)**unr_bet
            unr_gam = 1.88 ; unr_sig = 0.2285*10**(2.*unr_gam-4.) #A-D#unr_sig = 0.0; unr_gam = 0.0 #not needed for the powerlaw fall speed

    
    
    sph_gam = 2.0
    sph_sig = np.pi/4.
    print "ADmD_spec_string",ADmD_spec_string
    Dth = (rhoi * np.pi/6 * 1./unr_alf)**(1./(unr_bet-3))
    mth = np.pi/6. * rhoi * Dth**3 #from McSnows mo_mass2diam.f90
    return mth,unr_alf,unr_bet,rhoi,rhol,Dth,unr_sig,unr_gam,sph_sig,sph_gam

def read_init_vals(MC_dir):
    '''
    read heights of initialization and maximum qi+qs
    '''
    import csv
    import os
    with open(MC_dir + "/input/init_vals.txt","rb") as txtfile: #http://effbot.org/zone/python-with-statement.htm explains what if is doing; open is a python build in
        initvals_reader = csv.reader(txtfile, delimiter=' ', quoting=csv.QUOTE_NONE, lineterminator=os.linesep) #quoting avoids '' for formatted string; lineterminator avoids problems with system dependend lineending format https://unix.stackexchange.com/questions/309154/strings-are-missing-after-concatenating-two-or-more-variable-string-in-bash?answertab=active#tab-top
        height_bounds = [0,10000] #set some default heigh-bounds
        for i_row,row_content in enumerate(initvals_reader):
            height_bounds[i_row] = int(row_content[0])
            
        return height_bounds
def kernel_estimate(D_SP_list,Dgrid,sigmai,weight="None",space='loge'): #taken and adapted from mo_output.f90
    '''
    calculate the kernel estimate (adapted from McSnow's mo_output routines (f.e. write_distributions_meltdegree)
    INPUT:  D_SP_list: list of diameters of SP
            Dgrid: array of diameter on which the kde is calculated
            sigmai: bandwidth of kde
            weight: weight applied during integration
    '''
    number_sp = len(D_SP_list)
    
    #radius from diameter
    D_vec = Dgrid #radius array for results
    D_SP_vec = D_SP_list #radius of SP
    N_D = np.zeros_like(D_vec)
    
    if not isinstance(weight,basestring): #if there is a real weight
        weight_total = sum(weight)
    expdiff_prefactor = 1./np.sqrt(2.*np.pi)/sigmai #calculate expdiff here to save computational time
    for i_rad,rad in enumerate(D_vec):
        for i_SP,r in enumerate(D_SP_vec): 
            #print i_SP,i_rad,rad,r
            #calculate weight
            if space=='loge':
                expdiff = expdiff_prefactor * np.exp(-0.5*((np.log(rad)-np.log(r))/sigmai)**2) #r and rad is already log-scaled
            elif space=='lin':
                expdiff = expdiff_prefactor * np.exp(-0.5*(((rad)-(r))/sigmai)**2) #r and rad is already log-scaled 
            elif space=='D2':
                expdiff = expdiff_prefactor * np.exp(-0.5*(((rad)**2-(r)**2)/sigmai)**2) #r and rad is already log-scaled
            #integrate over each SP
            if isinstance(weight,basestring): #if there is no weight
                N_D[i_rad] += expdiff/number_sp #ATTENTION: add sp%xi to this summation as in mo_output
            else:
                #print 'here',weight[i_SP],weight[i_SP]*expdiff/number_sp
                N_D[i_rad] += weight[i_SP]*expdiff/weight_total #ATTENTION: add sp%xi to this summation as in mo_output
            #print r,rad,expdiff

            
    return N_D

def kde_statsmodels_m(x, x_grid, bandwidth=0.2, **kwargs):
    
    from statsmodels.nonparametric.kernel_density import KDEMultivariate #for multivariate KDE
    """Multivariate Kernel Density Estimation with Statsmodels"""
    kde = KDEMultivariate(x, bw=np.array(bandwidth * np.ones_like(x)),
                          var_type='c', **kwargs)

    return kde.pdf(x_grid) #return the pdf evaluated at the entries of x_grid

def kde_bandwidth_estimated(x,x_grid,guessed_bandwidths):
    
    #the next 2 packages are to estimate the bandwith by a machine learning approach    
    from sklearn.neighbors import KernelDensity
    from sklearn.grid_search import GridSearchCV
    
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': guessed_bandwidths},
                    cv=5) # ...-fold cross-validation
    #evaluate the grid to find the best parameter
    grid.fit(x[:, None])
    #copy result of parameter estimate to "kde"
    kde = grid.best_estimator_
    

    #evaluate the grid to find the best parameter
    grid.fit(x[:, None])
    #copy result of parameter estimate to "kde"
    #grid.best_estimator_
    
    return kde.bandwidth,np.exp(kde.score_samples(x_grid[:, None])) #score_samples returns the log of the pdf -> we need to exp() it
