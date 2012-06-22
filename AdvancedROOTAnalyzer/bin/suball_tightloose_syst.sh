#!/bin/bash

LOGFILE=/user/mweber/systematics.log

# original default submission
echo "Start submitting default13" > ${LOGFILE}
submit.py default13

# systematics for relative isolation
for reliso in 03 05 06 08
do
    TYPE=reliso_${reliso}
    echo "Start submitting ${TYPE}" > ${LOGFILE}
    submit.py -a ${TYPE} ${TYPE} 
done

# systematics for a possible jet cut bias
for jetptmin in 50 60 80
do
    TYPE=jetptmin_${jetcut}    
    echo "Start submitting ${TYPE}" > ${LOGFILE}
    submit.py -a ${TYPE} ${TYPE}
done

# systematics for a possible bias due to Z peak scaling / Z modelling
# One does not need to rerun the analysis for this.
# Needs to be done offline -> take care in singlefake / doublefake
#for zpeak in 86_96 powheg
#do 
#    TYPE=zpeak_${zpeak}
#    echo "Start submitting ${TYPE}" > ${LOGFILE}
#    submit.py -a ${TYPE} ${TYPE}
#done

# systematics for values outside histogram range
for fakeratemethod in zero
do 
    TYPE=fakeratemethod_${fakeratemethod}
    echo "Start submitting ${TYPE}" > ${LOGFILE}
    submit.py -a ${TYPE} ${TYPE}
done

# systematics for a possible trigger bias
for trigger in singlemu mu8_jet40
do
    TYPE=triggerbias_${trigger}
    echo "Start submitting ${TYPE}" > ${LOGFILE}
    submit.py -a ${TYPE} ${TYPE} ${trigger}
done
