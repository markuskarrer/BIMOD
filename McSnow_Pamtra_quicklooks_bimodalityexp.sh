#!/bin/bash

# THIS SCRIPT DOES ALL THE WORK: 
#- compiling and running McSnow 
#- run PAMTRA on the output
#- plot the output of McSNow and PAMTRA
# 

#############
# you might not always need to do all together (McSnow+PAMTRA+ plotting)
#, so you can give some input arguments to deactivate some of the procedures: e.g. "bash McSnow_Pamtra_quicklooks_bimodalityexp.sh 1 0 0 0 0" will only compile McSnow
#INPUT ARGUMENTS: Boolean (0) dont do (1) do
#$1 RECOMPILE
#$2 RUN MCSNOW
#$3 CREATE NCDF FROM MCSNOW RUN
#$4 RUN PAMTRA
#$5 DO PLOTTING
#
#
# choose carefully:     - is recompiling necessary
#                       - which script to run (line below #select which script to run) and the namelist of this script (the script will create a new directory for different namelist parameters (but not for all?)
#                       -select timestep of output with tstep 
#############
if [ "$2" == "1" ]; then #McSnow can somehow not be run with open evince instances in certain folders
    killall evince #otherwise there is an "directory not empty" error
fi

set -e #exit if error occurs

#if cheops is shut down use old folder for tests
export MC=/home/mkarrer/Dokumente/McSnow_bimodal/MCSNOW #TODO: REPLACE BY <MC_PATH>
export MCexp=/data/optimice/McSnowoutput/BIMOD/ #TODO: THINK ABOUT THIS
export BIMOD=$(pwd -P) #current directory

export depos_from0mto=0 # [100 m] 0 for SDA 10 for SDAdeposition and SDARdeposition
export riming_fromdeposto=0 # [100 m]  0 for SDA 10 for SDAdeposition 20 for SDARdeposition
#define what this script should do
if [ -z "$1" ]; then
    execwhat="r0Mc0ncdf0Pam0pl1" #recompile, McSnow, create ncdf with corresponding averages, PAMTRA, overview_panel #set 1 to activate and 0 to deactivate one of these steps
else #the following allows the user to define what the script should do by (currently 5) booleans
    execwhat="r"$1"Mc"$2"ncdf"$3"Pam"$4"pl"$5
    echo $execwhat
fi
######
# choose setup (testcase)
######
###
#use next line to plot McSnow runs with different geometries
###
export McSnow_testcase=("bimodal2mode") 

#switch of processes (for McSnow this is so far not used)
export switch_off_processes_list=("_") # "_noicesnow_nosnowself" "_nosnowself" "_noicesnow") # _noiceself, _noicesnow, _nosnowself

#something just needed for SB
export model_setup="_" #dont mess with SB settings if you are not running SB
export model_setup_specifier_onestring="_"
for testcase in "${McSnow_testcase[@]}"
do
    #loop over different namelist settings (relevant for running McSnow not for plotting)
    ssat_array=(0) #(0 10 20 30 40 50)  #supersaturation over ice: [1/1000] -> 1 is 0.001
    stick_array=(2) #sticking efficiency: 0: E_s=1, 1: PK97, 2: Connolly12, 3: 0.5*Connolly12
    ncl_array=(10) #(10 20 50) #nucleation rate [3/100 SP/sm3] #setting a high numer gets expensive (CPU-time and memory!); compensate this by a high multiplicity (xi0 in runscript: McSnow_runscripts/1d_bimodal2mode)
    nclmass_array=(1000) #(100 1000 5000) #
    McSnow_geom=2 # 1.constant 2.binary 3.monodep #we need binary here
    len_namelistcomb=${#ssat_array[@]}
    for (( i_ssat = 0 ; i_ssat < ${#ssat_array[@]} ; i_ssat++ )); do
        ssat=${ssat_array[$i_ssat]}
        for (( i_stick = 0 ; i_stick < ${#stick_array[@]} ; i_stick++ )); do
            stick=${stick_array[$i_stick]}
            for (( i_ncl = 0 ; i_ncl < ${#ncl_array[@]} ; i_ncl++ )); do
                ncl=${ncl_array[$i_ncl]}
                for (( i_nclmass = 0 ; i_nclmass < ${#nclmass_array[@]} ; i_nclmass++ )); do
                    nclmass=${nclmass_array[$i_nclmass]}
                #bash doesnt care about indentation; for looping over other variables just yank the for loop and the following line

        
    echo "#################"
    echo "analyzing testcase: "$testcase "with ssat=" $ssat ", stick=" $stick ",ncl=," $ncl ",nclmass," $nclmass
    echo "#################"
    export testcase; export McSnow_geom; export model_setup; export switch_off_processes #make them visible for called routines
            
    export tstep="300" #start of timestep analyzed by PAMTRA and plotting routine in minutes
    export tstep_end="600" #end timestep analyzed by PAMTRA and plotting routine in minutes (defines end of averaging period)
    export av_tstep="5" #average the McSnow output in the PAMTRA adaption by ... s #do not set this lower than dtc (f.e. "5"s)

    if [[ "$execwhat" == *r1* ]] ; then #recompile for a change in the model
        echo "############"
        echo "recompiling"
        echo "############"

        cd $MC

        #alternative if cheops is shut down
        make release_omp #alternative if cheops is shut down
    fi
    #get foldername from current runscript (even when it is not executed
    cd $BIMOD
    cp McSnow_runscripts/1d_${testcase%_*} $MC/run/1d
    export experiment=$($MC/run/1d "onlyname" $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto  $ncl $nclmass)
    echo "analyze experiment:" $experiment
    
    if [[ "$execwhat" == *Mc1* ]] ; then #run McSnow (f.e 1d-model)
        echo "############"
        echo "run McSnow"
        echo "############"

        cd $MC/run #alternative if cheops is shut down
        echo "SETTINGS: " $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto $testcase $habit $ncl $nclmass
        ./1d "fullrun" $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto  $ncl $nclmass
    fi

    if [[ "$execwhat" == *ncdf1* ]] ; then #run McSnow_dat2ncdf.py i.e. produce a netCDF in the right time averaging step
    #convert the mass2fr...dat to ncdf and derive temporal averages
        cd $BIMOD
        echo "converting McSnowoutput from .dat to .ncdf"
        python McSnow_dat2ncdf.py
    fi

    #switch between different versions of adaption #this is not in the next if clause because it is also needed in the plotting routines
    export adapt_version=3 #1-> histogram as input 2-> own category for each SP 3-> KDE as input

    if [[ "$execwhat" == *Pam1* ]] ; then #let PAMTRA run on output
        echo "############"
        echo "run Pamtra"
        echo "############"

        cd /home/mkarrer/Dokumente/pamtra/ #TODO: REPLACE by <PAM_path>
        #run PAMTRA
        python McSnow_adapt_bimodal.py
    fi
     if [[ "$execwhat" == *pl1* ]] ; then #produce quicklook of McSnow and the corresponding PAMTRA run
        echo "############"
        echo "plotting"
        echo "############"

        cd $BIMOD

        #the next two lines I used for plotting McSnow model output  (if you need something like this just ask me, they are not clean, but we might get them running)
        #python overview_panel.py 
        #python plot_fluxes_and_densities.py
        
        export plot_totalice="False"
        python bimodalitystudy_DWRs.py
        python plot_spectra_and_moments.py

    fi                 


    done #for ssat
    done #for stick
    done #for ncl
    done #for nclmass
done #for testcase


cd $BIMOD #go back to start directory

