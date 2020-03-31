
# coding: utf-8
#import packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#import other self-defined functions
import __postprocess_McSnow
import __general_utilities
from IPython.core.debugger import Tracer ; debug=Tracer()


# define some general functions which are useful for nice formatting of the plots

def proper_font_and_fig_size(number_of_plots,aspect_ratio=1./3.,legend_fontsize='large'):
    '''
    optimize the appearance of the plot (figure size, fonts)
    '''
    
    import matplotlib.pylab as pylab

    #increase font sizes
    params = {'legend.fontsize': legend_fontsize,
        'figure.figsize': (15, 5),
        'axes.labelsize': 'x-large', #size: Either a relative value of 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large' or an absolute font size, e.g., 12
        'axes.titlesize':'x-large',
        'xtick.labelsize':'x-large',
        'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)
    #define figure
    figsize_height = 1./aspect_ratio*number_of_plots
    fig, axes = plt.subplots(nrows=number_of_plots, ncols=1, figsize=(8.0,figsize_height))
    
    return fig,axes

# define some plotting functions

def plot1Dhist(x,bins,normed=1,facecolor='blue',alpha=1, 
               xlabel='',ylabel='counts',title='',logflag=0,axlims=[0]):
    '''
    input: x->main variable of the histogram, bins, normed-> normalize flag, number_of_bins, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    if logflag==1 or logflag==3:
        plt.xscale('log')
    if logflag==2 or logflag==3:
            logy=True
    else: logy=False
    
    # the histogram of the data
    n, bins, patches = plt.hist(x, bins=bins, normed=normed, facecolor=facecolor, alpha=alpha, log=logy)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.axis([10**-5, 1, 1, 10**6])
    if len(axlims)==4:
        plt.axis(axlims)
    plt.grid(True)
    
def plot1Dhistline(ax,x,counts,xlabel='',ylabel='counts',title='',linelabel='',logflag=0,xlims='None',ylims='None',linestyle='-',color='r'):
    '''
    input: x-> x-axes of bar-plot of the histogram, counts-> corresponding counts to x-axis, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    
    # the histogram of the data
    if logflag==0:
        ax.plot(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==1:
        ax.semilogx(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==2:
        ax.semilogy(x,counts,label=linelabel,linestyle=linestyle,color=color)
    elif logflag==3:
        ax.loglog(x, counts,label=linelabel,linestyle=linestyle,color=color) #, np.diff(x) , facecolor=facecolor, alpha=alpha, log=logy)
    
    if not ylims=='None':
        ax.set_ylim([ylims[0],ylims[1]])
    if not xlims=='None':
        ax.set_xlim([xlims[0],xlims[1]])    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax

def plotbar(x,counts,facecolor='blue',alpha=1, 
               xlabel='',ylabel='counts',title='',logflag=0,axlims=[0]):
    '''
    input: x-> x-axes of bar-plot of the histogram, counts-> corresponding counts to x-axis, facecolor, 
    alpha-> transparency value, xlabel,ylabel,title, logflag-> 0: no 1: logx 2:logy 3:loglog
    '''
    if logflag==1 or logflag==3:
        plt.xscale('log')
    if logflag==2 or logflag==3:
            logy=True
    else: logy=False
    
    # the histogram of the data
    plt.bar(x[:-1], counts[:-1], np.diff(x) , facecolor=facecolor, alpha=alpha, log=logy)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if len(axlims)==4:
        plt.axis(axlims)
    plt.grid(True)
    plt.rcParams.update({'font.size': 20})

def pcol_height_diam(ax,diam,heights,var,mincol=-999,maxcol=-999,logcol=0,cmap='viridis_r',collabel='',ticks=np.array([-999])):
    '''
    create a 2D-histogram with axis diameter & height : min&max, ticks and label of colorbar can be specified
    INPUT:  ax: axis handle (subplots or similar should be organized by calling script)
            diam: array of diameters
            height: array of heights
            var: variable which should be plotted now
            
            mincol,maxcol: minimum and maximum colors of the colorbar  (can be done automatically)
            logcol: 0-> linear colorbar 1-> logarithmic colorbar
            cmap: define colormap  (can be done automatically)
            collabel: label for colorbar
            ticks: specifiy ticks of colorbar (can be done automatically)
    '''

    if mincol==-999: #set range of colorbar automatically
        pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap) #norm=norm)
        if logcol==0:
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        elif logcol==1:
            if mincol<=0: print "Error: logarithmic colorbar and non-positive and mincol<= 0 is not possible"
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,norm=colors.LogNorm())
    else: #set range of colorbar via mincol, maxcol
        if logcol==0:
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,vmin=mincol,vmax=maxcol) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        elif logcol==1:
            if mincol<=0: print "Error: logarithmic colorbar and non-positive and mincol<= 0 is not possible"
            pcol = ax.pcolormesh(diam,heights,var,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=mincol,vmax=maxcol)) #,norm=colors.LogNorm(vmin=mincol,vmax=maxcol))
        
    #modify axis
    ax.set_xscale("log")
    #TODO: give this as an input
    ax.set_xlim([10**-5,10**-1])
    ax.set_ylim([0,heights[-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #plot colorbar
    
    
    if ticks[0]==-999: #automatic ticks
        col = plt.colorbar(pcol)
    else:#predefined ticks
        col = plt.colorbar(pcol,ticks=ticks)

    col.set_label(collabel, rotation=90)
    #set labels
    ax.set_xlabel('diameter / m ')
    ax.set_ylabel('height / m ')
    return ax
def plot_pamtra_Ze(ax,pamData,linestyle='-',marker=' ',nolabel=False,forcedcolor='auto',forcedlabel='auto',onlyfreq='None'):
    '''
    plot pamtra output
    INPUT:  pamData: dictionary with PAMTRA variables
            linestyle: force a linestyle (default is '-')
            marker: put optionally a marker to the lines
            nolabel: if True then dont add a legend entry for this line
            forcedcolor: if defined, select a color instead of the automatically selected one (which is defined by the frequency)
            forcedlabel: if defined, select a label instead of the automatically selected one (which is the frequency)

            onlyfreq: select one frequency (by its indice in the pamData["frequency"] list)
    '''

    #create figure for plotting of pamtra output

    #plot X-,Ka-,W-Band in three different colors
    #from IPython.core.debugger import Tracer ; Tracer()()
    for i in range(0,len(pamData["frequency"])): #all available frequencies in the pamtra files are plotted
        if not onlyfreq=='None':
            if not onlyfreq==i: #skip all non-selected frequencies
                continue
        if forcedcolor=='auto':
            color=np.array(['b','r','g'])[i]
        else:
            color=forcedcolor

        if linestyle=='-': #just create labels for the first series of frequencies (which should have '-' as a linestyle)
            if nolabel:
                label='__None'
            else:
                if forcedlabel=='auto':
                    label=label='{:5.1f}GHz'.format(pamData["frequency"][i])
                else:
                    label=forcedlabel

            ax.plot(pamData["Ze"][:,i],pamData["height"],color=color,linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=20,label=label)
        else:
            ax.plot(pamData["Ze"][:,i],pamData["height"],color=color,linestyle=linestyle,marker=marker,markerfacecolor='None',markevery=5)
    #set range
    ax.set_xlim([-40,55]) #range of Ze
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_xlabel("reflectivity / dBz")
    ax.set_ylabel("height / m")

    return ax

def plot_pamtra_highermoments(ax,pamData,linestyle='-',marker=' ',moment='all',forcedcolor='auto',forcedlabel='auto',onlyfreq='None'):
    ''' #TODO: there might not be a reason to distinguish between the reflectivity and other moments
    plot pamtra output
    INPUT: pamData: dictionary with PAMTRA variables (or the same naming)
            linestyle: force a linestyle (default is '-')
            marker: put optionally a marker to the lines
            nolabel: if True then dont add a legend entry for this line
            forcedcolor: if defined, select a color instead of the automatically selected one (which is defined by the frequency)
            forcedlabel: if defined, select a label instead of the automatically selected one (which is the frequency) #TODO: not used here

            onlyfreq: select one frequency (by its indice in the pamData["frequency"] list)
    '''

    if pamData["height"].ndim>1: #the observational data have more heights (for each velocity) which are all the same
        pamData["height"] = pamData["height"][:,0]
    #create figure for plotting of pamtra output
    #plot X-,Ka-,W-Band in three different colors
    for i in range(0,len(pamData["frequency"])): #all available frequencies in the pamtra files are plotted
        if not onlyfreq=='None':
            if not onlyfreq==i: #skip all non-selected frequencies
                continue
        #take given color or set automatically
        if forcedcolor=='auto':
            color=np.array(['b','r','g'])[i]
        else:
            color=forcedcolor
        markersymbol = '' #deactivate marker by default
        if linestyle=='-': #just create labels for the first series of frequencies (which should have '-' as a linestyle)
            
            if moment=='all' or moment=='swidth':
                ax.plot(pamData["Radar_SpectrumWidth"][:,i] ,pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')
            if moment=='all' or moment=='vDoppler':
                #ax.plot(pamData["Radar_MeanDopplerVel"][:,i],pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')
                ax.plot(pamData["Radar_MeanDopplerVel"][:,i],pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')

            if moment=='all' or moment=='skewn':
                if moment=='all': markersymbol='o'
                ax.plot(pamData["Radar_Skewness"][:,i]      ,pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')

        else:
            if moment=='all' or moment=='swidth':
                ax.plot(pamData["Radar_SpectrumWidth"][:,i] ,pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')
            if moment=='all' or moment=='vDoppler':
                if moment=='all': markersymbol='*'
                ax.plot(pamData["Radar_MeanDopplerVel"][:,i],pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')
            if moment=='all' or moment=='skewn':
                if moment=='all': markersymbol='o'
                ax.plot(pamData["Radar_Skewness"][:,i]      ,pamData["height"],color=color,linestyle=linestyle,marker=markersymbol,markerfacecolor='None',markevery=5,label='__None')

    #set range
    #from IPython.core.debugger import Tracer ; Tracer()()
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    if moment=='all':
        ax.set_xlabel("skewness (circles) & spectrum width & mean vDoppler (stars) / m s-1")
    elif moment=='swidth':
        ax.set_xlim([0,0.4]) 
        ax.set_xlabel("spectral width / m s-1")
    elif moment=='vDoppler':
        ax.set_xlim([0,2]) 
        ax.set_xlabel("mean vDoppler / m s-1")
    elif moment=='skewn':
        ax.set_xlim([-2,2]) 
        ax.set_xlabel("skewness / 1")
    ax.set_ylabel("height / m")
    #plt.grid()

    return ax

def plot_pamtra_spectrogram(ax,pamData,freq=35.5,cmap='viridis_r'):
    '''
    plot sprectogram from pamtra output
    INPUT: pamData: dictionary with PAMTRA variables
    freq: frequency, from which the spectrum should be used
    cmap: name of colormap
    '''
    import sys
    #choose spectrum of right frequency
    freqindex = np.where(abs(pamData["frequency"]-freq)<0.01)[0] #we need some tolerance because the frequency is not always saved correctly
    if freqindex.size==0: #give error if frequency is not in file
        print "frequency: ", freq ,"is not in the PAMTRA output file"
        print "choose between: ",pamData["frequency"]," or change frequencies in run_pamtra during runtime"
        sys.exit(0)

    #get axis and spectrogram data from pamData
    min_shown = -45 #minimum value shown in the diagram
    Radar_Spectrum = np.squeeze(pamData["Radar_Spectrum"][:,freqindex,:]) #dimensions [height, nfft]
    Radar_Spectrum_masked = np.ma.array(Radar_Spectrum,mask=Radar_Spectrum<min_shown) #do not show small values
    Radar_Velocity = np.squeeze(pamData["Radar_Velocity"][freqindex,:])  #dimension [nfft]
    height = np.squeeze(pamData["height"])
    
    #plot
    pcol = ax.pcolormesh(-Radar_Velocity,height,Radar_Spectrum_masked,rasterized=True,cmap=cmap,vmin=min_shown) #norm=norm)

    #set range
    ax.set_xlim([-3,1]) #range of Ze
    ax.set_ylim([0,pamData["height"][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #plot labels and create colorbar
    ax.set_xlabel("Doppler velocity / m s-1")
    ax.set_ylabel("height / m")
    col = plt.colorbar(pcol)
    col.set_label("spectral power @" + str(freq) + "GHz / dB", rotation=90)
    #plot label of used frequency
    #plt.text(0.4, 4500, str(freq) + 'GHz', fontsize=12)
    
    return ax

def plot_waterfall(ax,pamData,freq=35.5,color='b',linestyle='-',vel_lim=[0,3],z_lim=[0,5000]):
    '''
    plot spectra from pamtra output as a waterfall plot
    INPUT: pamData: dictionary with PAMTRA variables #it does not have to be Pamtra (but the same nomenclatura of the keys)
    freq: frequency, from which the spectrum should be used
    color, linestyle: arguments of the plot() functions
    vel_lim: limits to the velocity
    z_lim: limits to the height-axis
    '''

    #choose spectrum of right frequency
    freqindex = np.where(abs(pamData["frequency"]-freq)<0.01)[0] #we need some tolerance because the frequency is not always saved correctly
    if freqindex.size==0: #give error if frequency is not in file
        print "frequency: ", freq ,"is not in the PAMTRA output file"
        print "choose between: ",pamData["frequency"]," or change frequencies in run_pamtra during runtime"
        sys.exit(0)
    #define number of heights in the waterfall plot
    dz_heights = 500. #interval of heights at which the spectra are plotted #n_heights = 10
    #get axis and spectrogram data from pamData
    min_shown = -25. #minimum value shown in the diagram
    max_shown = 30. #minimum value shown in the diagram
    Radar_Spectrum = np.squeeze(pamData["Radar_Spectrum"][:,freqindex,:]) #dimensions [height, nfft]
   
    Radar_Velocity = np.squeeze(pamData["Radar_Velocity"][freqindex,:])  #dimension [nfft]
    #normalize the Radar Spectrum
    Radar_Spectrum_lin = 10.**(Radar_Spectrum/10.)
    Radar_Spectrum_lin_normalized = Radar_Spectrum_lin[:,:-1]/abs(np.diff(Radar_Velocity))
    Radar_Spectrum = 10.*np.log10(Radar_Spectrum_lin_normalized)
    #crop the Radar_Velocity vector to be compatible with the Radar Spectrum shape
    if Radar_Velocity.ndim==1: 
        Radar_Velocity = Radar_Velocity[:-1]
    else: #some observational data have a varying Radar Velocity with height
        Radar_Velocity = Radar_Velocity[:,:-1]
    #from IPython.core.debugger import Tracer ; Tracer()()
    
    height = np.squeeze(pamData["height"])
    heightmax = 10000. #max(height)
    if height.ndim>1: #the observational data have the heightvector for each velocity bin saved
        height = height[:,0]   #the heights should be the same for each velocity bin, so we take just 0 here
    for index_plot_heights,heights_dz in enumerate(range(int(heightmax),int(dz_heights),-int(dz_heights))):
        plot_height, index_radar_height = __general_utilities.find_nearest(height, heights_dz)

        if Radar_Velocity.ndim>1: #this is the case if the Radar_Velocity bins are changing with height
            Radar_Velocity_now = Radar_Velocity[index_radar_height,:] 
        else: 
            Radar_Velocity_now = Radar_Velocity

        #save the base line for the spectrogram of each height
        baseline = heights_dz #ATTENTION: here plotted is now not exactly the height from the radar data file#plot_height #height[i_height]

        topline = plot_height+dz_heights #height[i_height+height.shape[0]/n_heights]

        #define a scaling facor for dBz into the height y-axis
        mult_dBz = dz_heights/(max_shown-min_shown) #max((np.max(Radar_Spectrum[i_height])-min_shown),0.1) #2*n_heights #multiplicator for dBz on the y (height) axis to make it visible and not too large
        #mask out everything below mask_value and interpolate additional values to mask_value
        [Radar_Velocity_masked_and_interpolated,Radar_Spectrum_masked_and_interpolated] = __general_utilities.mask_and_interpolate(Radar_Velocity_now,Radar_Spectrum[index_radar_height],mask_value=min_shown,mask_smaller=True) #do not show small values

        #plot the spectrum
        ax.plot(Radar_Velocity_masked_and_interpolated,baseline + (Radar_Spectrum_masked_and_interpolated-min_shown)*mult_dBz,color=color,linestyle=linestyle)
        #plot the limits of the heights and label them with the corresponding dBz value
        ax.axhline(y=baseline,color='k',linestyle='--',linewidth=0.5)
    dBz_scale_xpos = 0.05
    ax.text(dBz_scale_xpos,baseline+0.05*(topline-baseline),str(min_shown) + " dBz")
    ax.text(dBz_scale_xpos,topline-0.15*(topline-baseline),str(dz_heights/mult_dBz+min_shown) + " dBz")
    
    ax.annotate(s='', xy=(dBz_scale_xpos,baseline), xytext=(dBz_scale_xpos,topline), arrowprops=dict(arrowstyle='<->'))
    #set range
    ax.set_xlim([vel_lim[0],vel_lim[1]]) #range of Ze
    ax.set_ylim([z_lim[0],z_lim[1]])
    #plot labels and create colorbar
    ax.set_xlabel("Doppler velocity / m s-1")
    ax.set_ylabel("layer top height / m")

    
    return ax

def plot_McSnows_vt_in_spectrogram(ax,pamData,SP,experiment,cmap='viridis_r'):
    '''
    plot number of RP in spectrogram bins
    INPUT: pamData: dictionary with PAMTRA variables
    SP dictionary with McSnows variables for each SP
    experiment: string describing the experiment (needed for calculating Vbox)
    cmap: name of colormap
    '''

    #get axis and spectrogram data from pamData
    Radar_Velocity = np.squeeze(pamData["Radar_Velocity"][0,:])  #dimension [nfft]
    height = np.squeeze(pamData["height"])
    binned_val=dict()
    binned_val["counts_spectro"] = np.zeros([len(height)-1,len(Radar_Velocity)]) #multiplicity weighted fall speed
    #get the number of particles (SP*sp_multiplicity) at each bin (binned by height and sp_diameter)
    for i in range(0,len(height)-1):
        for j in range(0,len(Radar_Velocity)-1):
                    condition_in_bin = np.logical_and(
                                    np.logical_and(Radar_Velocity[j]<=SP["vt"],SP["vt"]<Radar_Velocity[j+1]),
                                    np.logical_and(height[i]<=SP["height"],SP["height"]<height[i+1]),
                                    )
                    binned_val["counts_spectro"][i,j] = np.sum(np.where(condition_in_bin,SP["xi"],0))
    
    #calculate volume of box to derive number density later
    zres = height[1]-height[0]
    Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
    numdens_spectro=binned_val["counts_spectro"]/Vbox
    #plot
    masked_counts_spectr = np.ma.array(numdens_spectro,mask=numdens_spectro<=0)
    pcol = ax.pcolormesh(-Radar_Velocity,height,masked_counts_spectr,rasterized=True,cmap=cmap,norm=colors.LogNorm(vmin=10,vmax=10**5)) #norm=norm)
    #from IPython.core.debugger import Tracer ; Tracer()()

    #set range
    ax.set_xlim([-3,1]) #range of Ze
    
    #change height limits and ticks
    ax.set_ylim([0,height[-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))

    #plot labels and create colorbar
    ax.set_xlabel("fall speed / m s-1")
    ax.set_ylabel("height / m")
    col = plt.colorbar(pcol)
    col.set_label("number density / m-3", rotation=90)
    return ax

def plot_twomom_moments(ax,ax2,twomom,i_timestep,add_Dmean=True):
    '''
    plot moments (number density and mass density) over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            add_Dmean (boolean): add also Dmean into plot
    '''

    import matplotlib.ticker
    #plot mixing ratios
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    markerlist = ["o","^","s","*","x","d"]
    
    axisq = [];axisDmean = []
    for i,cat in enumerate(['c','i','r','s','g','h']):
        axisq += ax.semilogx(twomom['q' + cat][i_timestep,:],twomom['heights'],color='r',marker=markerlist[i],markevery=20,markerfacecolor='None')
        if add_Dmean:
            axisDmean += ax2.semilogx(twomom['D_mean_' + cat][i_timestep,:],twomom['heights'],color='k',linestyle='',marker=markerlist[i],markevery=20,markerfacecolor='None')

    #plot number concentrations
    axisqn = []
    for i,cat in enumerate(['c','i','r','s','g','h']):
        axisqn += ax2.semilogx(twomom['qn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None')

    #change height limits
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_xlabel("mixing ratio / kg m-3", color="r")
    ax.tick_params(axis='x', colors='r')    
    ax2.set_xlabel("mean diameter (symbols only) / m     number concentration / m-3",color="k")
    ax.set_ylabel("height / m")
    #set a fix range of the x-axis
    xmin=-9;xmax=1
    x2min=-5;x2max=9
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
    #set xticks
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1))
    ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    for label in ax2.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    
    #add legend
    ax.legend(axisqn,["cloud w.","cloud ice","rain","snow","graupel","hail"],loc='center right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return ax

def plot_fluxes(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=2,forced_linestyle='None',forced_markerMC='',nocatdist=False):
    '''
    plot number and mass fluxes of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            mass_num_flag: boolean which determines if the mass or the number flux should be plotted
            forced_linestyle: force a specific linestyle (otherwise the linestyle is selected automatically by either the dictname (twomom or hei2massdens) or the key_name (std); 
            nocatdist: dont distinguish between categories

    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    #colors[-1] = 'red'
    colors=["blue","green","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all Sb categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']);    twomom['q_all'] = np.zeros_like(twomom['qni']);     twomom['fn_all'] = np.zeros_like(twomom['qni']);     twomom['f_all'] = np.zeros_like(twomom['qni'])
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
        twomom['fn_all']+= twomom['fn' + cat]
        twomom['f_all']+= twomom['f' + cat]
    McSnow_plot_only_every = 4 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    if mass_num_flag==0:  #plot number flux
        axisqn = []
        print twomom.keys()
        for i,cat in enumerate(SB_specs):
            if forced_linestyle=='None' or forced_linestyle=='-':
                if any(twomom['fn' + cat][i_timestep,:]>0):
                    axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i],linestyle='-')
            else:
                if any(twomom['fn' + cat][i_timestep,:]>0):
                    axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color=colors[i],linestyle=forced_linestyle)
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if nocatdist and cat!="": continue #continue if we dont want a classification 
            if any(hei2massdens["Fn" + cat][:])>0:
                if forced_markerMC=='': #only add a label in case of no marker, which is the first MCSnow sensrun
                    axisqn += ax2.semilogx(hei2massdens["Fn" + cat][:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
                
                else:
                    axisqn += ax2.semilogx(hei2massdens["Fn" + cat][:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle='--',marker=forced_markerMC,markevery=4)
    if mass_num_flag==1:     #plot mass flux
        axisq = [];axisDmean = []
        if mass_num_flag==1:
            color='r'
        elif mass_num_flag==2:
            color='r'
        for i,cat in enumerate(SB_specs):        
            if forced_linestyle=='None' or forced_linestyle=='-':
                if any(twomom['fn' + cat][i_timestep,:]>0):
                    axisq += ax2.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i],linestyle='-')
            else:
                if any(twomom['fn' + cat][i_timestep,:]>0):
                    axisq += ax2.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=colors[i],linestyle=forced_linestyle)
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if nocatdist and cat!="": continue #continue if we dont want a classification 
            if any(hei2massdens["Fm" + cat][:])>0:
                if forced_markerMC=='': #only add a label in case of no marker, which is the first MCSnow sensrun
                    axisq += ax.semilogx(hei2massdens["Fm" + cat][:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
                else:
                    axisq += ax.semilogx(hei2massdens["Fm" + cat][:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle='--',marker=forced_markerMC,markevery=4)
    elif mass_num_flag==2:
        for i,cat in enumerate(SB_specs):
            if any(twomom['f' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=color,marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])

        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if any(hei2massdens["Fm" + cat][:])>0:
                axisq += ax2.semilogx(hei2massdens["Md" + cat][:],hei2massdens['z'],color='b',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])

        for i,cat in enumerate(SB_specs):
            if any(twomom['fn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if any(hei2massdens["Fn" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Fn" + cat][:],hei2massdens['z'],color='g',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    mass_min=-8; mass_max=0
    num_min=2; num_max=6
    if mass_num_flag==2:
        ax.set_xlabel("mass flux density / kg m-2 s-1", color="r")
        ax.tick_params(axis='x', colors='r')    
        ax2.set_xlabel("number flux density / m-2 s-1",color="k")
        #set a fix range of the x-axis
        xmin=mass_min;xmax=mass_max
        x2min=num_min;x2max=num_max
    elif mass_num_flag==0:
        ax.set_xlabel("$F_N$ [m-2 s-1]",color="k")
        xmin=num_min;xmax=num_max
    elif mass_num_flag==1:
        ax.set_xlabel("$F_m$ [kg m-2 s-1]",color="k")
        xmin=mass_min;xmax=mass_max
    else:
        print "error: mass_num_flag in plot_fluxes must be in [0,1,2]"


    #set a fix range of the x-axis
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1)) #set xticks
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    
    #add legend
    if mass_num_flag in [0,2]:
        ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    else: #if there is no qn in the plot label q
        ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return ax

def plot_MC_profiles(ax,hei2massdens,i_timestep,var_flag=0,forced_linestyle='-',forced_markerMC='',top_height=5000,plot_totalice=True,catonly=None,show_label=True,force_gray=False):
    '''
    plot number and mass fluxes of McSnow over height
    INPUT:  ax: first x-axis
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            var_flag: which variable should be plotted 0:number 1: mass 2: number flux 3: mass flux
            forced_linestyle: force a specific linestyle (otherwise the linestyle is selected automatically by either the dictname (twomom or hei2massdens) or the key_name (std); 
            plot_totalice: plot a line for the sum of the categories
            catonly: plot a specific category (if None plot all avalable)
    '''
    import matplotlib.ticker
    colors=["blue","green","orange","brown","y","red"]
    markerlist = ["^","s","*","x","d",""]
    MC_specs = ["_mm1","_unr","_grp","_rimed","_liq",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    if show_label:
        labellistMC = ["Nmono=1","Nmono>1","graupel","rimed","liquid","all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    else:
        labellistMC = ["__","__","__","__","__","__"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all Sb categories

    McSnow_plot_only_every = 4 #plot not every height point of the noisy McSnow fluxes (for better visibility)
        

    for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
        if force_gray:
            colors[j]="gray"
        if (cat=="" and not plot_totalice): #plot all categories but not the total
            continue
        elif (catonly!=None and catonly!=cat): #plot only a specific category
            continue
	if var_flag==0:  #plot number flux
            plot_var=hei2massdens["Nd" + cat]
        if var_flag==1:  #plot mass flux
            plot_var=hei2massdens["Md" + cat]
        if var_flag==2:  #plot number flux
            plot_var=hei2massdens["Fn" + cat]
        if var_flag==3:  #plot mass flux
            plot_var=hei2massdens["Fm" + cat]
        if var_flag==4:  #mean mass
            plot_var=hei2massdens["Md" + cat]/hei2massdens["Nd" + cat]    

        if any(hei2massdens["Md" + cat]>2e-6): #there might be some numerical rime therefore here is not >0
            if forced_markerMC=='' and forced_linestyle=='-': #only add a label in case of no marker, which is the first MCSnow sensrun
                ax.semilogx(plot_var[:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='-',marker=forced_markerMC)
            else:
                ax.semilogx(plot_var[:-1:McSnow_plot_only_every],hei2massdens['z'][:-1:McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle=forced_linestyle,marker=forced_markerMC,markevery=4)
                
                
    #from IPython.core.debugger import Tracer ; Tracer()()
    #change height limits and set y-label
    ax.set_ylim([0,top_height])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height [m]")
    
    
    #add labels and legend
    num_min=-1; num_max=7
    mass_min=-9; mass_max=-2
    nmass_min=-10; nmass_max=-4
    fnum_min=-2; fnum_max=5
    fmass_min=-9; fmass_max=-3
    #fmass_min=-4; fmass_max=-3
    '''
    fmass_min=-8; fmass_max=0
    fnum_min=2; fnum_max=6
    mass_min=-8; mass_max=-1
    num_min=0; num_max=6
    nmass_min=-10; nmass_max=-4
    '''
    if var_flag==0:
        ax.set_xlabel("number density / m-3",color="k")
        xmin=num_min;xmax=num_max
    elif var_flag==1:
        ax.set_xlabel("mass density / kg m-3",color="k")
        xmin=mass_min;xmax=mass_max
    elif var_flag==2:
        ax.set_xlabel("$F_N$ [m-2 s-1]",color="k")
        xmin=fnum_min;xmax=fnum_max
    elif var_flag==3:
        ax.set_xlabel("$F_m$ [kg m-2 s-1]",color="k")
        xmin=fmass_min;xmax=fmass_max
    elif var_flag==4:
        ax.set_xlabel("$m_{mean}$ [kg]",color="k")
        xmin=nmass_min;xmax=nmass_max
    else:
        print "error: var_flag in plot_fluxes must be in [0,1,2,3,4]"


    #set a fix range of the x-axis
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1)) #set xticks
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    
    #add legend
    if var_flag in (1,3):
        ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    else: #if there is no qn in the plot label q
        ax.legend(loc='lower right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return ax

def plot_moments(ax,ax2,twomom,hei2massdens,i_timestep,mass_num_flag=2,forced_linestyle='None',forced_markerMC='',nocatdist=False):
    '''
    plot number and mass concentration of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            mass_num_flag: boolean which determines if the mass or the number flux should be plotted
            forced_linestyle: force a specific linestyle (otherwise the linestyle is selected automatically by either the dictname (twomom or hei2massdens) or the key_name (std); 
            nocatdist: dont distinguish between categories
    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    colors=["blue","green","brown","orange","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','g','h','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr","_grp","_rimed",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","graupel","hail","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","graupel","rimed","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all Sb categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']);    twomom['q_all'] = np.zeros_like(twomom['qni']);     twomom['fn_all'] = np.zeros_like(twomom['qni']);     twomom['f_all'] = np.zeros_like(twomom['qni'])
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
        
    McSnow_plot_only_every = 4 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    if mass_num_flag==0:  #plot number concentrations
        axisqn = []
        for i,cat in enumerate(SB_specs):
            if ('qn' + cat) in twomom.keys() and any(twomom['qn' + cat][i_timestep,:]>0):
                if forced_linestyle=='None' or forced_linestyle=='-':
                    axisqn += ax2.semilogx(twomom['qn' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
                else:
                    axisqn += ax2.semilogx(twomom['qn' + cat][i_timestep,:],twomom['heights'],color=colors[i],label='__None',linestyle=forced_linestyle)
                if ('qn' + cat + '_std') in twomom.keys(): #add +- standard deviation in plot if available
                    ax2.semilogx(twomom['qn' + cat][i_timestep,:]-twomom['qn' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    ax2.semilogx(twomom['qn' + cat][i_timestep,:]+twomom['qn' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if nocatdist and cat!="": continue #continue if we dont want a classification 
            if ("Nd" + cat) in hei2massdens.keys() and any(hei2massdens["Nd" + cat][:])>0:
                if forced_markerMC=='': #only add a label in case of no marker, which is the first MCSnow sensrun
                    #from IPython.core.debugger import Tracer ; Tracer()()
                    axisqn += ax2.semilogx(hei2massdens["Nd" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--',marker=forced_markerMC)
                else:
                    axisqn += ax2.semilogx(hei2massdens["Nd" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle='--',marker=forced_markerMC,markevery=4)
    if mass_num_flag==1:     #plot mixing ratios
        axisq = [];axisDmean = []
        #if mass_num_flag==1:
        #    color='r'
        #elif mass_num_flag==2:
        #    color='r'
        for i,cat in enumerate(SB_specs):
            if ("q" + cat) in twomom.keys() and any(twomom['q' + cat][i_timestep,:]>0):
                if forced_linestyle=='None' or forced_linestyle=='-':
                    axisq += ax.semilogx(twomom['q' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
                else:
                    axisq += ax.semilogx(twomom['q' + cat][i_timestep,:],twomom['heights'],color=colors[i],label='__None',linestyle=forced_linestyle)

                if ('q' + cat + '_std') in twomom.keys(): #add +- standard deviation in plot if available
                    ax2.semilogx(twomom['q' + cat][i_timestep,:]-twomom['q' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    ax2.semilogx(twomom['q' + cat][i_timestep,:]+twomom['q' + cat + '_std'][i_timestep,:],twomom['heights'],color=colors[i],linestyle='--',label='_dontshowinlegend')
                    #print "q" + cat,twomom['q' + cat]
        for j,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if nocatdist and cat!="": continue #continue if we dont want a classification 
            if ("Md" + cat) in hei2massdens.keys() and any(hei2massdens["Md" + cat][:])>0:
                if forced_markerMC=='': #only add a label in case of no marker, which is the first MCSnow sensrun
                    axisq += ax.semilogx(hei2massdens["Md" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
                else:
                    axisq += ax.semilogx(hei2massdens["Md" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle='--',marker=forced_markerMC,markevery=4)

    elif mass_num_flag==2:
        for i,cat in enumerate(SB_specs):
            if any(twomom['f' + cat][i_timestep,:])>0:
                axisq += ax.semilogx(twomom['f' + cat][i_timestep,:],twomom['heights'],color=color,marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up            
            if any(hei2massdens["Md" + cat][:])>0:
                axisq += ax2.semilogx(hei2massdens["Md" + cat][:],hei2massdens['z'],color='b',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
        for i,cat in enumerate(SB_specs):
            if any(twomom['fn' + cat][i_timestep,:])>0:
                axisqn += ax2.semilogx(twomom['fn' + cat][i_timestep,:],twomom['heights'],color='k',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistSB[i])
        for i,cat in enumerate(MC_specs): #choose "" as an entry to get all summed up
            if any(hei2massdens["Nd" + cat][:])>0:
                axisqn += ax2.semilogx(hei2massdens["Nd" + cat][:],hei2massdens['z'],color='g',linestyle='--',marker=markerlist[i],markevery=20,markerfacecolor='None',label=labellistMC[i])
    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    mass_min=-8; mass_max=-1
    num_min=0; num_max=6
    if mass_num_flag==2:
        ax.set_xlabel("mass density / kg m-3", color="r")
        ax.tick_params(axis='x', colors='r')    
        ax2.set_xlabel("number density / m-3",color="k")
        #set a fix range of the x-axis
        xmin=mass_min;xmax=mass_max
        x2min=num_min;x2max=num_max
    elif mass_num_flag==0:
        ax.set_xlabel("number density / m-3",color="k")
        xmin=num_min;xmax=num_max
    elif mass_num_flag==1:
        ax.set_xlabel("mass density / kg m-3",color="k")
        xmin=mass_min;xmax=mass_max
    else:
        print "error: mass_num_flag in plot_fluxes must be in [0,1,2]"



    #set a fix range of the x-axis
    ax.set_xlim([10**xmin,10**xmax]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(xmin,xmax,xmax-xmin+1)) #set xticks
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    
    labels = ax2.get_xticklabels()
    #add legend
    ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)
    
    return ax

def plot_normmix(ax,ax2,twomom,hei2massdens,i_timestep,forced_linestyle='None',forced_markerMC='',nocatdist=False):
    '''
    plot normalized mixing ratio/ mean particle mass (ratio of mass density and number density) of SB and McSnow over height
    INPUT:  ax: first x-axis
            ax2: second x-axis
            twomom: dictionary with (McSnows) two-moment variables
            hei2massdens: dictionary with McSnow variables
            i_timestep: timestep used (the dicts contain all timesteps of the simulation output; i_timestep is defined by governing scripts)
            forced_linestyle: force a specific linestyle (otherwise the linestyle is selected automatically by either the dictname (twomom or hei2massdens) or the key_name (std); ATTENTION: forced linestyle comes with no labels
            nocatdist: dont distinguish between categories
    '''
    import matplotlib.ticker
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    #colors[-1] = 'red'
    colors=["blue","green","red"]
    markerlist = ["^","s","*","x","d",""]
    SB_specs = ['i','s','_all'] #specify here what categories should be used; all would be ['i','r','s','g','h']
    MC_specs = ["_mm1","_unr",""] #specify here what categories should be used; all would be ["_mm1","_unr","_grp","_liq","_rimed",""]#
    labellistSB = ["cloud ice","snow","SB all"] #check consistency with SB_specs; all would be ["cloud ice","rain","snow","SB graupel","hail"]
    labellistMC = ["pristine","unrimed agg.","MC all"] #check consistency with MC_specs; all would be ["pristine","unrimed agg.","MC graupel","liquid","rimed","MC all"]
    #sum up all SB categories
    twomom['qn_all'] = np.zeros_like(twomom['qni']); twomom['q_all'] = np.zeros_like(twomom['qni']);
    for i,cat in enumerate(['i','s']): #['i','s','g','h','r']):#TODO: include other categories when analyzing rimed cases 
        twomom['qn_all']+= twomom['qn' + cat]
        twomom['q_all']+= twomom['q' + cat]
    #calculate normalized mixing ratios
    #first for the sum of all categories
    twomom['nq_all'] = twomom['q_all']/twomom['qn_all']
    hei2massdens['nq_all'] = hei2massdens['Md']/hei2massdens['Nd']
    #second for all categories individually
    for i,cat in enumerate(SB_specs):
        twomom['nq' + cat] = twomom['q' + cat]/twomom['qn' + cat]
    for j,cat in enumerate(MC_specs):
        if any(hei2massdens["Md" + cat]>1e-7):
            hei2massdens['nq' + cat] = hei2massdens["Md" + cat]/hei2massdens["Nd" + cat]
    
    axisnq = [];axisDmean = [] #initialize handles
    McSnow_plot_only_every = 4 #plot not every height point of the noisy McSnow fluxes (for better visibility)
    for i,cat in enumerate(SB_specs): #['i','s','g','h','r','_all']):#TODO: include other categories when analyzing rimed cases 
        if ('nq' + cat) in twomom.keys() and any(twomom['nq' + cat][i_timestep,:]>0):
            if forced_linestyle=='None' or forced_linestyle=='-':
                axisnq += ax2.semilogx(twomom['nq' + cat][i_timestep,:],twomom['heights'],color=colors[i],label=labellistSB[i])
            else:
                axisnq += ax2.semilogx(twomom['nq' + cat][i_timestep,:],twomom['heights'],color=colors[i],label='__None',linestyle=forced_linestyle)
                
    for j,cat in enumerate(MC_specs):  #["_mm1","_unr","_grp","_liq","_rimed",""]): #choose "" as an entry to get all summed up #TODO: include other categories when analyzing rimed cases 
        if nocatdist and cat!="": continue #continue if we dont want a classification 
        if ("nq" + cat) in hei2massdens.keys() and any(hei2massdens["nq" + cat][:]>0):
            if forced_markerMC=='': #only add a label in case of no marker, which is the first MCSnow sensrun
                axisnq += ax2.semilogx(hei2massdens["nq" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label=labellistMC[j],linestyle='--')
            else:
                axisnq += ax2.semilogx(hei2massdens["nq" + cat][::McSnow_plot_only_every],hei2massdens['z'][::McSnow_plot_only_every],color=colors[j],label='_dontshowinlegend',linestyle='--',marker=forced_markerMC,markevery=4)
    #change height limits and set y-label
    ax.set_ylim([0,twomom["heights"][0]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    ax.set_ylabel("height / m")
    
    
    #add labels and legend
    nmass_min=-10; nmass_max=-4
    ax.set_xlabel("mean particle mass / kg") #"normalized mixing ratio / kg",color="k")


    #set a fix range of the x-axis
    ax.set_xlim([10**nmass_min,10**nmass_max]) #[ax.get_xlim()[0],ax.get_xlim()[1]*100])
    ax.set_xticks(np.logspace(nmass_min,nmass_max,nmass_max-nmass_min+1)) #set xticks
    '''
    #remove every second ticklabel
    for label in ax.xaxis.get_ticklabels()[::2]:
        label.set_visible(False)
    if mass_num_flag==2: #do the same for the second x-axis if existent
        ax2.set_xlim([10**x2min,10**x2max]) #[ax2.get_xlim()[0],ax2.get_xlim()[1]*1000])
        ax2.set_xticks(np.logspace(x2min,x2max,x2max-x2min+1))
        for label in ax2.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
    '''
    #labels = ax2.get_xticklabels()
    #add legend
    ax.legend(loc='upper right') #,loc='center left', bbox_to_anchor=(1, 0.5)) #position: the "center left" of the box is at (x=1,y=0.5, in relative coordinates of the whole plot)

    return ax

def plot_atmo(ax,ax2,atmo):
    '''
    plot atmospheric variables in one panel
    INPUT:  ax: first x-axis
            ax2: second x-axis
            atmo: dictionary with atmospheric variables
    '''
    #plot T, RHw and RHi in one panel
    handles = []
    if 'T' in atmo.keys():
        handles += ax.plot(atmo['T']-273.15,atmo['z'],color='r',label='T')
        ax.set_xlabel("temperature / $^\circ$C")
        if 'T_std' in atmo.keys():
            handles += ax.plot(atmo['T']-273.15-atmo['T_std'],atmo['z'],color='r',linestyle='--',label='__none')
            handles += ax.plot(atmo['T']-273.15-atmo['T_std'],atmo['z'],color='r',linestyle='--',label='__none')
    if 'rh' in atmo.keys():
        handles += ax2.plot(atmo["rh"],atmo['z'],color='b',label='RHw')
        ax2.plot(np.ones(atmo['z'].shape)*100,atmo['z'],color='grey',linestyle='--',label='_dontshowinlegend')
        ax2.set_xlabel("relative humidity / %")
        if 'rh_std' in atmo.keys():
            handles += ax2.plot(atmo["rh"]+atmo["rh_std"],atmo['z'],color='b',linestyle='--',label='__none')
            handles += ax2.plot(atmo["rh"]-atmo["rh_std"],atmo['z'],color='b',linestyle='--',label='__none')
    #plot rhi (if rh,psatw and psati present calculate it from them) elif rhi is there take it directly
    #if all (key in atmo.keys() for key in ("rh","psatw","psati")):
    #    handles += ax2.plot(atmo["rh"]*atmo["psatw"]/atmo["psati"],atmo['z'],color='b',linestyle='--',label='RHi')
    if 'rhi' in atmo.keys():
        handles += ax2.plot(atmo["rhi"],atmo['z'],color='g',linestyle='-',label='RHi')
        if 'rhi_std' in atmo.keys():
            handles += ax2.plot(atmo["rhi"]+atmo["rhi_std"],atmo['z'],color='g',linestyle='--',label='__none')
            handles += ax2.plot(atmo["rhi"]-atmo["rhi_std"],atmo['z'],color='g',linestyle='--',label='__none')
    #change height limits
    ax.set_ylim([0,np.nanmax(atmo['z'])])#atmo['z'][-1]])
    ax.set_yticks(np.arange(0,ax.get_ylim()[1]+1,1000))
    #add labels and legend
    ax.set_ylabel("height / m")
    #create space for legend
    ax.set_xlim([ax.get_xlim()[0],ax.get_xlim()[1]+0])
    low_RH=50 #define lower end of RH by hand
    ax2.set_xlim([low_RH,ax2.get_xlim()[1]+15])
    # added these three lines
    labs = [h.get_label() for h in handles]
    ax.legend(handles,labs,loc='upper right')
    #from IPython.core.debugger import Tracer ; Tracer()()
    return ax
 
def plot_vD_scatter(ax,v_SPlist,D_SPlist,xlog=False,ylog=False):
    '''
    plot velocity over diameter
    INPUT:  v_SPlist: fall speed for a list of SP
            D_SPlist: diameter for a list of SP
    '''
    #plot the lists against each other
    ax.scatter(D_SPlist,v_SPlist,s=1,rasterized=True)
    if xlog:
        ax.set_xscale("log")
        ax.set_xlim([1e-5,3e-2])
    if ylog:
        ax.set_yscale("log")
    #define labels
    ax.set_xlabel("diameter of SP / m") #"normalized mixing ratio / kg",color="k")
    ax.set_ylabel("fall speed of SP / m s-1") #"normalized mixing ratio / kg",color="k")
    #set xlimits
    if not xlog:
        ax.set_xlim(left=0)
    return ax
#below plots for 3D-data 
def pcolor_timeseries(fig,ax,time,height,var,varname='',time_formatted='None',lin_log=0,unit=''):
    '''
    plot timeseries (written for plotting ICON-meteogram)
    INPUT:  time-> timevector
            height-> heightvector
            var -> variable which should be plotted
    '''
    import matplotlib.colors as colors
    
    #####
    #variable specific settings
    #####
    #default settings #these are overwritten after "specifications for variable groups or individual variables"
    lin_log=0 #lin is default
    cmap = "viridis" #take viridis as default colorbar
    #if lin_log==0:
    varmin=var.min()
    varmax=var.max()
    #elif lin_log==1:
    #    varmin=max(var.min(),1e-100) #ATTENTION:max(..,1e-100) avoids only the crash of the script; define lower limit for variables for which zeros occur below
    #    varmax=var.max()
        
    ###
    #specifications for variable groups or individual variables
    ###
    if varname in ("qc","qr","qi","qs","qg","qh"): #these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e-8
        varmax = 1e-2
    if varname in ("qnc","qnr"): #these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e0
        varmax = 1e12
    if varname in ("qni","qns","qng","qnh"):#these are logscaled
        lin_log = 1 #set y-axis to logscale
        #set limits
        varmin = 1e0
        varmax = 1e6
    if varname=="rhi":
        #set limits
        varmin = 80
        varmax = 120
        cmap = "bwr"
    #normalize according to lin_log
    if lin_log==0:
        norm = colors.Normalize(vmin=varmin,vmax=varmax)
    elif lin_log==1:
        norm = colors.LogNorm(vmin=varmin,vmax=varmax)
    
    #plot
    ax_now = ax.pcolormesh(time,height,var,rasterized=True,norm=norm,cmap=cmap)#,vmin=mincol,vmax=maxcol,rasterized=True,norm=norm,cmap=new_cmap)
    #create colorbar
    cbar = plt.colorbar(ax_now,ax=ax) #, format='%d')
    cbar.ax.set_ylabel(varname + " / " + unit)
    #ticklabels
    #derive location of ticks dynamically depending on timeranges
    if np.sum((time_formatted.minute==0) & (time_formatted.second==0))>15: #label every 2 hours for long timeranges (f.e. 15 hours); specify hours by the last number
        ticklocations = np.where((time_formatted.minute==0) & (time_formatted.second==0) & ((time_formatted.hour%2)==0))[0]
    elif 2<np.sum((time_formatted.minute==0) & (time_formatted.second==0))<15: #label every hour for delta_t_hour = [2,15]
        ticklocations = np.where((time_formatted.minute==0) & (time_formatted.second==0))[0]
    else:
        print "not implemented x-axis labelling for short timeranges"
        sys.exit(1)
    #set location of ticks
    ax.set_xticks(time[ticklocations])
    #generate tick labels and set them
    timestrings = [ "{:02d}:{:02d}".format(time_formatted[i].hour,time_formatted[i].second) for i in ticklocations]
    ax.set_xticklabels(timestrings) #"{:2.0f}:{:2.0f}".format(time[ticklocations].hour.values,time[ticklocations].minute.values))
    #labelling
    ax.set_xlabel("time / UTC")
    ax.set_ylabel("height / m")

    return ax
