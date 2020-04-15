'''
this script includes some general utilities (not specific for a tool like McSnow, PAMTRA,...)
'''

#from IPython.core.debugger import Tracer ; Tracer()()


import numpy as np
import sys

def gen_shortest_strings_without_duplicate(count_len=4):
    '''
    generate all combinations which can be done with 3-4 string consistent off small,big letters and numbers
    '''
    import itertools #used to generate all permutations

    #allocate list of strings
    if count_len==3: #for a string of 3 we have 238328 possible entries
        n_entries=238328
        count_str = ['___']*n_entries
    elif count_len==4: #for a string of 4 we have 14776336 possible entries
        n_entries=14776336
        count_str = ['___']*n_entries
    elif 14776336<number_ofSP:
        print "check number of SP (exceeding anticipated number in this adaption; in __general_utilities.gen_shortest_strings_without_duplicate)"
        sys.exit(0)

    #generate list of strings with itertools
    #if __name__ == "__main__":
    chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

    i=0
    for item in itertools.product(chars, repeat=count_len):
        count_str[i] ="".join(item)
        i=i+1
    return count_str


def find_nearest(array, value):
    '''
    find element in ARRAY with smallest difference to VALUE
    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

################
# some functions to convert atmospheric variables into  each other
################

def get_SB_constants(constlist): #kwargs 
    '''
    returns some constants from the SB-scheme
    INPUT: **kwargs:  arbitrary number of keyword arguments depending on which the necessary constants from the SB scheme are returned 
    '''
    constants = dict()
    
    if "T_3" in constlist:
        constants["T_3"] =  273.15 #! [K]     melting temperature of ice/snow #from mo_atmotypes
    if "A_w" in constlist:
        constants["A_w"]  = 1.72693882e1 #, &  !..constant in saturation pressure - liquid #from mo_atmotypes
    if "B_w" in constlist:
        constants["B_w"]  = 3.58600000e1  #!..constant in saturation pressure - liquid #from mo_atmotypes
    if "A_i" in constlist:
        constants["A_i"]  = 2.18745584e1 #, &  !..constant in saturation pressure - ice #from mo_atmotypes
    if "B_i" in constlist:
        constants["B_i"]  = 7.66000000e0  #!..constant in saturation pressure - ice #from mo_atmotypes
    if "e_3" in constlist:
        constants["e_3"]  = 6.10780000e2  #!..saturation pressure at T = T_3#from mo_atmotypes
        
    return constants

def get_Pamtra_constants(constlist):
    '''
    returns some constants from PAMTRA
    INPUT: **kwargs:  arbitrary number of keyword arguments depending on which the necessary constants from PAMTRA are returned 
    '''
    constants = dict()
    
    if "r_d" in constlist:
        constants["r_d"] = 287.0596736665907    #! gas constant of dry air
    if "r_v" in constlist:
        constants["r_v"] = 461.5249933083879    #! gas constant of water vapor
        
    return constants
#taken from SB
def calc_rhi(atmo): #calculate rh over ice from rh over water and saturation pressure equation in SB
    '''
    calculate the relative humidity with respect to ice from the relative humidity with respect to water and temperature
    INPUT:  atmo: dictionary which must contain numpy arrays including the keys rh and T
    '''

    psat_w = psatw(atmo)
    psat_i = psati(atmo)
    rhi = psat_w/psat_i*atmo["rh"]
    
    return rhi

def calc_rh(atmo): #calculate rh over water from rh over ice(rhi) and saturation pressure equation in SB
    '''
    calculate the relative humidity with respect to ice from the relative humidity with respect to water and temperature
    INPUT:  atmo: dictionary which must contain numpy arrays including the keys rh and T
    '''

    psat_w = psatw(atmo)
    psat_i = psati(atmo)
    rh = psat_i/psat_w*atmo["rhi"]

    return rh

def psatw(atmo): #PURE ELEMENTAL FUNCTION psat_water( T ) result ( psatw ) in mo_atmo.f90
    '''
    calculate relative humidity with respect to water from temperature
    INPUT: atmo: dictionary which must contain numpy arrays including the key T
    '''
    constants = get_SB_constants(constlist=["T_3","A_w","B_w","e_3"])
    return constants["e_3"] * np.exp(constants["A_w"] * (atmo["T"]-constants["T_3"]) / (atmo["T"]-constants["B_w"]))

def psati(atmo): #PURE ELEMENTAL FUNCTION psat_water( T ) result ( psatw ) in mo_atmo.f90
    '''
    calculate relative humidity with respect to ice from temperature
    INPUT: atmo: dictionary which must contain numpy arrays including the key T
    '''
    constants = get_SB_constants(constlist=["T_3","A_i","B_i","e_3"])
    print 'used constants for conversions of atmospheric variables', constants
    return constants["e_3"] * np.exp(constants["A_i"] * (atmo["T"]-constants["T_3"]) / (atmo["T"]-constants["B_i"]))

#taken from McSNow

def psat_water(T):
    
    constants = get_SB_constants(constlist=["T_3","A_w","B_w","e_3"]) #these are the same in ICON-SB and McSnow

    psatw = constants["e_3"] * np.exp( constants["A_w"] * (T - constants["T_3"]) / (T - constants["B_w"]) )

    return psatw

def calc_qv_from_rh(rh,p,T):
    
    constants = get_Pamtra_constants(constlist=["r_d","r_v"]) #these are the same in PAMTRA and McSnow
    qvsat = constants["r_d"]/constants["r_v"] * psat_water(T) / (p - (1.-constants["r_d"]/constants["r_v"])*psat_water(T))
    qv = rh/100.*qvsat
    
    return qv

#taken from PAMTRA
def q2abs(spec_var,qv,t,p,q_all_hydro="NaN"):#function q2abs(spec_var,t,p,qv,q_all_hydro) in src/conversions.f90
    '''
    convert "spec var" from [kg kg-1] to [kg m-3]
    INPUT:  spec_var: variable which is converted
            t: temperature in K
            p: pressure in hPa
    '''
    #        use kinds
    #        use constants, only: r_d, r_v
    constants = get_Pamtra_constants(constlist=["r_d","r_v"])


    #real(kind=dbl), intent(in) :: spec_var,& ! specific variable to convert [kg/kg]
    #qv,& 
    #q_all_hydro

    if q_all_hydro=="NaN":
        q_all_hydro = 0 #neglect 
        print "WARNING: neglect hydrometeor weight in conversion from specific to absolute quantities, because not given as input (in __general_utilities.q2abs() )"
    p_Pa = p*100.
    q2abs = spec_var*p/(constants["r_d"]*(1.+(constants["r_v"]/constants["r_d"]-1.)*qv-q_all_hydro)*t)
    return q2abs

def rh2vap(temp_p,pres_p,rh):
    '''
    calculate vapour mixing ratio "hum_massmix" from INPUT:
        INPUT:  temperature in K
                pressure in Pa
                relative humidity rh in %
    '''
    #reverse vap2rh from pamtras conversions.f90

    tpt = 273.16               # triple point temperature
    estpt = 611.14             # saturation vapor triple point temperature
    r_d = 287.0596736665907    #gas constant of dry air
    r_v = 461.5249933083879    #gas constant of water vapor
    mmv = 18.0153e-3           # molar mass of vapor
    mmd = 28.9644e-3           # molar mass of dry air
    vapor_hc  = 2.5008e+6      # vaporization heat constant
    sublim_hc  = 2.8345e+6     # sublimation heat constant

    '''
    !     XPABSM air pressure in Pa
    !     ZTEMP air temperature in K
    !     XRM water vapor mass mixing ratio kg/kg


    REAL(kind=dbl) :: XCPV               ! Cpv (vapor)
    REAL(kind=dbl) :: XCL,XCI            ! Cl (liquid), Ci (ice)
    REAL(kind=dbl) :: XLMTT              ! Melting heat constant
    REAL(kind=dbl) :: XALPW,XBETAW,XGAMW ! Const saturation vapor pressure  (liquid)
    REAL(kind=dbl) :: XALPI,XBETAI,XGAMI ! Consts saturation vapor pressure  (solid ice)
    real(kind=dbl) :: ztemp,xrm,xpabsm,zwork31,zwork32
    '''

    zwork32 = 0
    XCPV   = 4. * r_v
    XCL    = 4.218e+3
    XCI    = 2.106e+3
    XLMTT  = sublim_hc - vapor_hc
    XGAMW  = (XCL - XCPV) / r_v
    XBETAW = (vapor_hc/r_v) + (XGAMW * tpt)
    XALPW  = np.log(sublim_hc) + (XBETAW /tpt) + (XGAMW *np.log(tpt))
    XGAMI  = (XCI - XCPV) / r_v
    XBETAI = (sublim_hc/r_v) + (XGAMI * tpt)
    XALPI  = np.log(estpt) + (XBETAI /tpt) + (XGAMI *np.log(tpt))

    ZTEMP = temp_p
    #XRM = hum_massmix
    XPABSM = pres_p

    #!     * humidit351 relative par rapport 340 l'eau liquide
    ''' #original code: but we want to get qv from rh, so the other way round 
    if (ZTEMP >= tpt):
        ZWORK31=np.exp(XALPW - XBETAW/ZTEMP - XGAMW*np.log(ZTEMP))
        ZWORK31=(mmv/mmd)*ZWORK31/(XPABSM-ZWORK31)
        ZWORK32=100.*XRM/ZWORK31
    elif (ZTEMP < tpt):
        #!     * humidit351 relative par rapport 340 la glace
        ZWORK31=np.exp(XALPI - XBETAI/ZTEMP - XGAMI*np.log(ZTEMP))
        ZWORK31=(mmv/mmd)*ZWORK31/(XPABSM-ZWORK31)
        ZWORK32=100.*XRM/ZWORK31
            #if you skip the second line
            #ZWORK32=100.*XRM/((mmv/mmd)*ZWORK31/(XPABSM-ZWORK31)) (Eq 1)

    vapor2rh=ZWORK32
    '''
    #if (ZTEMP >= tpt): #naming as here: https://github.com/PyAOS/aoslib/blob/master/aoslib/src/spechum.f
    #	#calculate vapor pressure
    #	ZWORK31 = np.exp(XALPW - XBETAW/ZTEMP - XGAMW*np.log(ZTEMP))
    #elif (ZTEMP < tpt):
    #	ZWORK31 = np.exp(XALPI- XBETAI/ZTEMP - XGAMI*np.log(ZTEMP))

    ZWORK31 = np.where(ZTEMP >= tpt,np.exp(XALPW - XBETAW/ZTEMP - XGAMW*np.log(ZTEMP)),np.exp(XALPI- XBETAI/ZTEMP - XGAMI*np.log(ZTEMP))) #translated the if clause above into pythonic handling of arrays
    #reverse Eq 1
    hum_massmix = 1./100.* rh * ((mmv/mmd)*ZWORK31/(XPABSM-ZWORK31))
    #from IPython.core.debugger import Tracer ; Tracer()()
    return hum_massmix

def rh2mixr(RH,p,T):
    '''purpose: conversion relative humidity (unitless) to mixing ratio [kg/kg]'''
    Mw=18.0160 # molecular weight of water
    Md=28.9660 # molecular weight of dry air
    R =  8.31432E3 # gas constant
    Rd = R/Md # specific gas constant for dry air
    Rv = R/Mw # specific gas constant for vapour
    Lv = 2.5e6 # heat release for condensation of water vapour [J kg-1]
    eps = Mw/Md
    es = esat(T)
    print "es",es,"RH",RH,"p",p
    return Mw/Md*RH*es/(p-RH*es)
def esat(T):
    ''' get sateration pressure (units [Pa]) for a given air temperature (units [K])'''
    from numpy import log10
    TK = 273.15
    e1 = 101325.0
    logTTK = log10(T/TK)
    esat =  e1*10**(10.79586*(1-TK/T)-5.02808*logTTK+ 1.50474*1e-4*(1.-10**(-8.29692*(T/TK-1)))+ 0.42873*1e-3*(10**(4.76955*(1-TK/T))-1)-2.2195983) 
    return esat
################
#end:  some functions to convert atmospheric variables into  each other
################

def mask_and_interpolate(x,y,mask_value,mask_smaller=True):
    '''
    this functions adds additional data points to the x and y array at exactly the "mask_value" and masks all the smaller (bigger) value for mask_smaller=True (mask_smaller=False)
    '''
    
    #get a boolean vector where the points of the sign change (in this case: bigger or smaller than mask_value) is marked by True (1 element shorter than the original vector)
    sign_changes = (np.diff(np.sign(y-mask_value)) != 0)
    #get indices of the sign change
    sign_changes_indices = np.where(sign_changes)[0]
    sign_changes_posneg = np.diff(np.sign(y-mask_value))/2 #-1 for a negative +1 for a positive change  
    x_at_maskedvalue = np.zeros(len(sign_changes_indices))
    for i in range(0,len(sign_changes_indices)):
        #for each sign change take the surrounding x- and y-value
        if sign_changes_posneg[sign_changes_indices[i]]>0: #sign changes from negative to positive: insert before
            x_around = x[sign_changes_indices[i]:(sign_changes_indices[i]+2)]
            y_around = y[sign_changes_indices[i]:(sign_changes_indices[i]+2)]
        elif sign_changes_posneg[sign_changes_indices[i]]<0: #sign changes from positive to negative: insert after
            x_around = x[(sign_changes_indices[i]+1):(sign_changes_indices[i]+3)]
            y_around = y[(sign_changes_indices[i]+1):(sign_changes_indices[i]+3)]
        else: #no sign change at all (f.e. the spectral power is always below min_shown (here: mask_value)
            x_around = [np.nan,np.nan]
            y_around = [np.nan,np.nan]
        #calculate interpolated values for y=mask_value
        x_at_maskedvalue[i] = (x_around[1]-x_around[0])/(y_around[1]-y_around[0])*(mask_value-y_around[0])+x_around[0]
        #insert the interpolated values
        if y_around[1]-y_around[0]>0: #sign changes from negative to positive: insert before
            x = np.insert(x,sign_changes_indices[i]+1,x_at_maskedvalue[i])
            y = np.insert(y,sign_changes_indices[i]+1,mask_value)
        elif y_around[1]-y_around[0]<0: #sign changes from positive to negative: insert after
            x = np.insert(x,sign_changes_indices[i]+2,x_at_maskedvalue[i])
            y = np.insert(y,sign_changes_indices[i]+2,mask_value)
    #mask the arrays
    y_masked = np.ma.array(y,mask=y<mask_value)
    x_masked = np.ma.array(x,mask=y<mask_value)

    return x_masked,y_masked
