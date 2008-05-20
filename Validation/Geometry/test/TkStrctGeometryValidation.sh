#! /bin/bash

# the CMSSW_x_y_z to compare with actual one:
export Version=$1

export here=$PWD
export reportFile=TkStrctGeometryValidation.log

cd $here
echo "------------------------------------------------------" | tee    $reportFile
echo "   VALIDATION OF THE TRACKER GEOMETRY"                  | tee -a $reportFile
echo "     NEW = $CMSSW_VERSION vs OLD = $Version"            | tee -a $reportFile
echo "------------------------------------------------------" | tee -a $reportFile

echo "Working area:" $here | tee -a $reportFile
eval `scramv1 runtime -sh`

# OLD: from CMSSW_1_6_* no more updated:
#export referenceDir=/afs/cern.ch/cms/data/CMSSW/Validation/Geometry/reference/Tracker
# NEW: the reference files can be chosen
#export referenceDir=/afs/cern.ch/cms/performance/tracker/activities/validation/ReferenceFiles/$Version/Geometry
#echo "Reference area:" $referenceDir | tee -a $reportFile
#
## Create Images/ directory if it does not exist
#if [ ! -d Images ]; then
#    echo "Creating directory Images/" | tee -a $reportFile
#    mkdir Images
#    echo "...done" | tee -a $reportFile
#fi
##
#
# Download the source file
if [ ! -e single_neutrino.random.dat ]; then
    echo "Download the Monte Carlo source file..." | tee -a $reportFile
    wget `cat $CMSSW_RELEASE_BASE/src/Validation/Geometry/data/download.url`
    echo "...done" | tee -a $reportFile
fi
#

## Download the reference files and rename them to 'old'
#echo "Download the reference 'old' files..." | tee -a $reportFile
#cp $referenceDir/matbdg_TkStrct.root       matbdg_TkStrct_old.root 
#echo "...done" | tee -a $reportFile
##

# Run all the Tracker scripts and rename files as 'new'
echo "Run all the scripts to produce the 'new' files..." | tee -a $reportFile
#
echo "Running Tracker Structure..." | tee -a $reportFile
rm -rf TkStrct.txt
cmsRun $CMSSW_BASE/src/Validation/Geometry/test/runP_TkStrct.cfg       > TkStrct.txt
echo "...done" | tee -a $reportFile
#
#cp matbdg_TkStrct.root       matbdg_TkStrct_new.root 
echo "...done" | tee -a $reportFile
#

## Produce the 'new' plots
#echo "Run the Tracker macro MaterialBudget.C to produce the 'new' plots..." | tee -a $reportFile
#root -b -q 'MaterialBudget.C("TkStrct")'
#echo "...done" | tee -a $reportFile
##
#
## Compare 'old' and 'new' plots
#echo "Run the Tracker macro TrackerMaterialBudgetComparison.C to compare 'old and 'new' plots..." | tee -a $reportFile
#root -b -q 'TrackerMaterialBudgetComparison.C("TkStrct")'
#echo "...done" | tee -a $reportFile
##
#
### Run the Tracker ModuleInfo analyzer (to compare position/orientation of Tracker Modules)
##echo "Run the Tracker ModuleInfo analyzer to print Tracker Module info (position/orientation)..." | tee -a $reportFile
##cmsRun $CMSSW_RELEASE_BASE/src/Geometry/TrackerGeometryBuilder/test/trackerModuleInfo.cfg
##echo "...done" | tee -a $reportFile
###
##
### Compare the ModuleInfo.log file with the reference one
##echo "Compare the ModuleInfo.log (Tracker Module position/orientation) file with the reference one..." | tee -a $reportFile
##if [ -e diff_info.temp ]; then
##    rm -rf diff_info.temp
##fi
###
##diff ModuleInfo.log $referenceDir/ModuleInfo.log > diff_info.temp
##if [ -s diff_info.temp ]; then
##    echo "WARNING: the module position/orientation is changed, check diff_info.temp file for details" | tee -a $reportFile
##else
##    echo "Tracker Module position/orientation OK" | tee -a $reportFile
##fi
##echo "...done" | tee -a $reportFile
###
##
### Run the Module Numbering (only Microstrip) check algorithm and print the tail
##echo "Run the Tracker ModuleNumbering analyzer to print Tracker Numbering check..." | tee -a $reportFile
##cmsRun $CMSSW_RELEASE_BASE/src/Geometry/TrackerNumberingBuilder/test/trackerModuleNumbering.cfg
##echo "TRACKER MICROSTRIP NUMBERING... LOOK AT THE RESULTS" | tee -a $reportFile
##tail -7 ModuleNumbering.log | tee -a $reportFile
##if [ -e num.log ]; then
##    rm -rf num.log
##fi
##tail -7 ModuleNumbering.log > num.log
##echo "...done" | tee -a $reportFile
###
##
### Compare the ModuleNumbering.dat file with the reference one
##echo "Compare the ModuleNumbering.dat (Tracker Module position/orientation) file with the reference one..." | tee -a $reportFile
##if [ -e diff_num.temp ]; then
##    rm -rf diff_num.temp
##fi
###
##diff ModuleNumbering.dat $referenceDir/ModuleNumbering.dat > diff_num.temp
##if [ -s diff_num.temp ]; then
##    echo "WARNING: the module numbering is changed, check diff_num.temp file for details" | tee -a $reportFile
##else
##    echo "Tracker Module numbering OK" | tee -a $reportFile
##fi
##echo "...done" | tee -a $reportFile
###
##
### Compare the TrackerNumberingComparison.C, to compare the ModuleNumbering.dat file with the reference, element-by-element mapping both files 
##echo "Run the TrackerNumberingComparison.C macro" | tee -a $reportFile
##cp $referenceDir/ModuleNumbering.dat ModuleNumbering_reference.dat
##root -b -q 'TrackerNumberingComparison.C("ModuleNumbering.dat","ModuleNumbering_reference.dat","NumberingInfo.log")'
##if [ -s NumberingInfo.log ]; then
##    echo "ERROR: a failure in the numbering scheme, see NumberingInfo.log" | tee -a $reportFile
##else
##    echo "Tracker Numbering Scheme OK" | tee -a $reportFile
##fi
##echo "...done" | tee -a $reportFile
###
##
### New test: Check Overlap
##echo "Run the Tracker Check Overlap test" | tee -a $reportFile
##if [ -e trackerOverlap.log ]; then
##    rm -rf trackerOverlap.log
##fi
##cmsRun $CMSSW_RELEASE_BASE/src/Validation/CheckOverlap/test/data/runTracker.cfg > trackerOverlap.log
##grep -A4 'WARNING - ' trackerOverlap.log | tee -a $reportFile
##echo "...done" | tee -a $reportFile
##
###
##echo "Run dddreport..." | tee -a $reportFile
##$CMSSW_RELEASE_BASE/src/DetectorDescription/RegressionTest/scripts/dddreport.sh > dddreport.log
##echo "...done" | tee -a $reportFile
###
##echo "Run domcount..." | tee -a $reportFile
##$CMSSW_RELEASE_BASE/src/DetectorDescription/RegressionTest/scripts/domcount.sh > domcount.log
##echo "...done" | tee -a $reportFile
###
#
echo "TRACKER GEOMETRY VALIDATION ENDED... LOOK AT THE RESULTS" | tee -a $reportFile
