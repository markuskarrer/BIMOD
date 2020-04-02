#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import __fallspeed_relations

from IPython.core.debugger import Tracer ; debug=Tracer()
def massPlawCgs2SI(aCgs, bCgs):
    """
        conversion os the power law
        parameters from CGS to SI 
        
        m(D) = aD^b
    """
    
    aSi = (100**(bCgs)) * (aCgs/1000.) 
    bSi = bCgs
    
    return aSi, bSi
def areaPlawCgs2SI(aCgs, bCgs):
    """
        conversion os the power law
        parameters from CGS to SI 
        
        A(D) = aD^b
    """
    
    aSi = (100**(bCgs)) * (aCgs/10000.) 
    bSi = bCgs
    
    return aSi, bSi



def gammaFunc(particle, massAreaParam, psdParamTM):
    """ 
        calculation of the gamma psd for 
        a selected particle
        
    """
    
    
    N0 = psdParamTM[particle]['N0']
    mu = psdParamTM[particle]['mu']
    Lambda = psdParamTM[particle]['Lambda']
    binSize = psdParamTM[particle]['binSize'] #[microns]
    
    minSize = massAreaParam[particle]['D_range'][0] * 1e-6 #[microns]
    maxSize = massAreaParam[particle]['D_range'][1] * 1e-6#[microns]
    
     #[microns]
    diameters = np.arange(minSize, maxSize, binSize) #[microns]
    diameterBins = diameters[:-1] + 0.5 * np.diff(diameters) #[microns]

    D = diameterBins #* 1e-6#[meters]
    
    N_D = N0 * D**(mu) * np.exp(-Lambda * D)
    
    return N_D, D

def getParticleSettings(particle, massAreaParam, psdParamTM):
    
    """
        Full settings of a chosen particle
    
    """
    
    am, bm = massPlawCgs2SI(massAreaParam[particle]['alpha'],
                            massAreaParam[particle]['beta'])
    aA, bA = areaPlawCgs2SI(massAreaParam[particle]['gamma'],
                            massAreaParam[particle]['sigma'])

    psd, D = gammaFunc(particle, massAreaParam, psdParamTM)
    qn = np.sum(psd * psdParamTM[particle]['binSize']) # [#/m^4]
    mass = am * D**bm #[kg]
    qm = np.sum(psd * mass * psdParamTM[particle]['binSize']) # total mass [kg]
    
    settingsDic  = {'am':am, # [SI]
                    'bm':bm, # [SI]
                    'aA':aA,
                    'bA':bA,
                    'gam':psdParamTM[particle]['gam'],
                    'mu':psdParamTM[particle]['mu'],
                    'LAMBDA':psdParamTM[particle]['Lambda'],
                    'N0':psdParamTM[particle]['N0'], 
                    'psd':psd, # [#/m^4]
                    'mass':mass, # mass per bins [Kg]
                    'D':D, # diameters [m]
                    'dD':psdParamTM[particle]['binSize'], # bin size [m]
                    'qn':qn, # total concentration [# /m^3]
                    'qm':qm, # total mass [Kg / m^3]
                    }
    
    return settingsDic
    
def getIconParam(pyTmParam):
    """
       convertion from distribution in size 
       to a distribution in mass
      
       psd parameters for mass gamma
       f(x) = A x^nu exp(-lambda x^muI )
    """
    

    bgeo = 1/ pyTmParam['bm']
    ageo = 1/(pyTmParam['am']**bgeo)
    print pyTmParam["gam"]
    lambd = pyTmParam['LAMBDA'] * ageo**pyTmParam['gam']
    muI = pyTmParam['gam'] * bgeo

    nu = pyTmParam['mu'] * bgeo + bgeo - 1
    A = pyTmParam['N0'] * bgeo * ageo ** (pyTmParam['mu']+1)
    
    paramDic = {'ageo':ageo,
                'bgeo':bgeo,
                'lambda':lambd,
                'muI':muI,
                'nu':nu,
                'A':A}
    
    return paramDic
def main(display_params=False):
    #(mitchell 1996)
    #[D_range: mu m, mass: g (alpha, beta), area: cm^2 (gamma, sigma)]
    massAreaParam = {'Hexagonal_plates2':{'alpha':0.00739,'beta':2.45, 'gamma':0.65, 'sigma':2.00, 'D_range':(50,3000)},
                 'Crystal_with_sector_like_branches2':{'alpha':0.00142,'beta':2.02, 'gamma':0.55, 'sigma':1.97, 'D_range':(50,3000)},                 
                 'Crystal_with_sector_like_branches2_expND':{'alpha':0.00142,'beta':2.02, 'gamma':0.55, 'sigma':1.97, 'D_range':(100,3000)},                 
                 'Aggregates_of_side_planes':{'alpha':0.0033,'beta':2.2, 'gamma':0.2285, 'sigma':1.88, 'D_range':(300,4100)},
                }

    psdParamTM = {'Hexagonal_plates2':{'N0':1.50641526924e+20, 'mu':3,'gam':1.0, 'Lambda':1.5e4, 'binSize':10e-6},
              'Crystal_with_sector_like_branches2':{'N0':1316944120.15, 'mu':0.3,'gam':1.0, 'Lambda':6e3, 'binSize':10e-6},                 
              'Crystal_with_sector_like_branches2_expND':{'N0':1.4e7, 'mu':0,'gam':2.02, 'Lambda':3.8e6, 'binSize':10e-6},                 
              'Aggregates_of_side_planes':{'N0':1.5e37, 'mu':9,'gam':1.0, 'Lambda':9e3, 'binSize':10e-6},
             }

    #particle='Crystal_with_sector_like_branches2' # mode 2
    #particle='Hexagonal_plates2' # mode 2
    # particle='Aggregates_of_side_planes' # mode 1

    particleSettings = dict()
    particleInconParam = dict()
    for particle in ['Crystal_with_sector_like_branches2','Aggregates_of_side_planes']:
        particleSettings[particle] = getParticleSettings(particle, massAreaParam, psdParamTM)

        particleInconParam[particle] = getIconParam(particleSettings[particle])
        if display_params:
            print '\n'
            print '##',particle 
            print 'size dist mass: N(m)=A*m**nu*exp(-lambda*m**muI) qm qn'
            print "A",particleInconParam[particle]["A"],"nu",particleInconParam[particle]["nu"],"lambda",particleInconParam[particle]["lambda"],"muI",particleInconParam[particle]["muI"],"qm",particleSettings[particle]["qm"],"qn",particleSettings[particle]["qn"],"qm/qn",particleSettings[particle]["qm"]/particleSettings[particle]["qn"],"Dmean",(particleSettings[particle]["qm"]/particleSettings[particle]["qn"]/particleSettings[particle]["am"])**(1./particleSettings[particle]["bm"])
            print 'particle prop m(D)=am*D**bm A(D)=aA*D**bA'
            print "am",particleSettings[particle]["am"],"bm",particleSettings[particle]["bm"],"aA",particleSettings[particle]["aA"],"bA",particleSettings[particle]["bA"]
            diam=np.logspace(-6,-2,100)
            plt.loglog(diam,particleSettings[particle]["am"]*diam**particleSettings[particle]["bm"])
            plt.xlabel("diam [m]")
            plt.ylabel("mass [kg]")
            #plt.show()
            print '\n'
            particleSettings[particle]["A"] = particleSettings[particle]["aA"]*particleSettings[particle]["D"]**particleSettings[particle]["bA"]
            particleSettings[particle]["m"] = particleSettings[particle]["am"]*particleSettings[particle]["D"]**particleSettings[particle]["bm"]
            vterm = __fallspeed_relations.calc_vterm("HW10",particleSettings[particle]["m"],particleSettings[particle]["D"],particleSettings[particle]["A"])
            fig,axes = plt.subplots(1,2)
            #m-D
            axes[0].loglog(particleSettings[particle]["D"],particleSettings[particle]["m"])
            axes[0].set_xlabel("D [m]")
            axes[0].set_ylabel("mass [kg]")
            #v-D
            axes[1].semilogx(particleSettings[particle]["D"],vterm)
            axes[1].set_xlabel("D [m]")
            axes[1].set_ylabel("vterm [m/s]")
            axes[1].text(1e-3,0.2,particle)

            #plt.show()

            #--- plotting PSDs ---#
            #plt.semilogy(particleSettings['D'],particleSettings['psd'], 
            #             label = 'qn: {0}'.format(particleSettings['qn']))
            # plt.xlim(0,3e-3*1e3)
            # plt.ylim(1e-10,1e9)
            #plt.ylabel('# [m-4]')
            #plt.xlabel('diameter [m]')
            #plt.legend()
            #plt.show()

            #plt.plot(particleSettings['D'], particleSettings['mass'], 
            #             label = 'qm: {0}'.format(particleSettings['qm']))
            # plt.xlim(0,3e-3*1e3)
            # plt.ylim(1e-10,1e9)
            #plt.ylabel('# [kg / m-4]')
            #plt.xlabel('diameter [m]')
            #plt.legend()
            #plt.show()
    return massAreaParam,particleSettings,particleInconParam

if __name__=="__main__":

    main(display_params=True)
