#!/bin/bash

LOGFILE=/user/mweber/systematics.log

#Suffix determines which jobs to start: 
# "" (empty)    : Default job start
# "_singlefake" : Singlefake job start
# "_doublefake" : Doublefake job start
SUFFIX="$1"

# Second argument determines if jobs are to be submit or collected
# "-s" means submit jobs
# "-c" means collect jobs
ACTION="$2"

newroot()
{
    source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh 
    source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.01/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh
}

oldroot()
{
    cd ~/CMSSW_4_2_8_patch7/src
    eval `scramv1 runtime -sh`
}

mycollect()
{
   # skip configuration file
   if [[ "$1" == "-a" ]]
   then
       shift ; shift
   fi
   newroot 
   echo collect.py -f $1 $2 
   collect.py -f $1 $2 ;
   oldroot
   collect.py -r $1 $2 
   sleep 600 
   newroot 
   collect.py $1 $2
}

action()
{
    echo "ACTION: ${ACTION}"
    if [[ "${ACTION}" == "-s" ]] ; then
	submit.py $*
    elif [[ "${ACTION}" == "-c" ]] ; then
	mycollect $*
    else
	echo "You must specify -s (submit) or -c (collect) as second argument"
    fi
}

## original default submission
#for name in default13
#do
#    TYPE=${name}${SUFFIX}
#    echo "Start submitting ${TYPE}" >> ${LOGFILE}
#    action -a ${TYPE} ${TYPE} 
#done

## systematics for relative isolation
#for reliso in 03 05 06 08
#do
#    TYPE=reliso_${reliso}${SUFFIX}
#    echo "Start submitting ${TYPE}" >> ${LOGFILE}
#    action -a ${TYPE} ${TYPE} 
#done

# systematics for a possible jet cut bias
for jetptmin in 50 60 80
do
    TYPE=jetptmin_${jetptmin}${SUFFIX}
    echo "Start submitting ${TYPE}" >> ${LOGFILE}
    action -a ${TYPE} ${TYPE}
done

# systematics for a possible bias due to Z peak scaling / Z modelling
# One does not need to rerun the analysis for this.
# Needs to be done offline ->> take care in singlefake / doublefake
#for zpeak in 86_96 powheg
#do 
#    TYPE=zpeak_${zpeak}
#    echo "Start submitting ${TYPE}" >> ${LOGFILE}
#    action -a ${TYPE} ${TYPE}
#done

# systematics for values outside histogram range
for fakeratemethod in zero
do 
    TYPE=fakeratemethod_${fakeratemethod}${SUFFIX}
    echo "Start submitting ${TYPE}" >> ${LOGFILE}
    action -a ${TYPE} ${TYPE}
done

## systematics for a possible trigger bias
#for trigger in singlemu mu8_jet40
#do
#    TYPE=triggerbias_${trigger}${SUFFIX}
#    echo "Start submitting ${TYPE}" >> ${LOGFILE}
#    action -a ${TYPE} ${TYPE} ${trigger}
#done
