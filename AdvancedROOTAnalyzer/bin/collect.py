#!/bin/env python
import os
import sys
import stat
import optparse
import ara
from ROOT import TEnv

######################################################################
# functions

def check(job):
    """Check log files for errors and warnings"""
    global options
    global myEnv

    print "Checking job", job.name
    good = True
    basedir = myEnv.GetValue("basedir", ".");
    myDir = basedir + '/' + options.selection + '/' + job.period

    errors = [ "exception", "segmentation", "error", "err:" ]
    warnings = [ "warning", "wrn" ]
    requirements = [ "End executing" ]
    
    # List all log files
    for filename in os.listdir(myDir):
        if job.name in filename:
            if "stdout.log" in filename:
                if not ara.check_log(myDir+'/'+filename, errors, warnings, requirements):
                    return False
            if "stderr.log" in filename:
                if not ara.check_log(myDir+'/'+filename, errors, warnings, None):
                    return False                
    return True

def join(job):
    global options
    global myEnv

    basedir = myEnv.GetValue("basedir", ".");
    myDir = basedir + '/' + options.selection + '/' + job.period
    os.chdir(myDir)
    # input file
    try:
        filelist = ara.getCommandOutput2('ls ' + myDir + '/' + job.name + '_*.root').splitlines()
    except RuntimeError:
        print "ERROR: No result files found"
        return False
    # get file dates
    lastDate = 0
    for inputFile in filelist:
        file_stats = os.stat(inputFile)
        lastDate = max(file_stats[stat.ST_MTIME], file_stats[stat.ST_CTIME], lastDate)
    # output file
    joinFile = myDir + '/' + job.name + ".root"
    joinDate = 0
    try:
        file_stats = os.stat(joinFile)
        joinDate = max(file_stats[stat.ST_MTIME], file_stats[stat.ST_CTIME], joinDate)
    except OSError:
        pass
    # Do not join if joinFile already newer than input files
    if (joinDate > lastDate):
        return True
    # now join
    if isinstance (filelist, str):
        os.system("hadd -f " + job.name + ".root " + " " + filelist)
    else:
        os.system("hadd -f " + job.name + ".root " + " ".join(filelist))
    return True

def rm(job):
    global options
    global myEnv

    basedir = myEnv.GetValue("basedir", ".");
    myDir = basedir + '/' + options.selection + '/' + job.period
    print "Removing output file", job.name + ".root"
    os.system('rm -f ' + myDir + '/' + job.name + ".root")

######################################################################
# main
def main():
    global options
    global myEnv

    usage = "usage: %prog [options] selection job | all"
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.configFileName)

    (options, args) = optParser.parse_args()
    if len(args) != 2:
        optParser.print_help()
        return 1

    # setup local arguments
    options.selection = args[0]
    options.job = args[1]

    # read global configuration file
    ara.config.read(options.cfgfile)

    # read ROOT configuration file 
    configPath = os.environ['ARASYS'] + '/config'
    myEnv = TEnv(configPath + '/' + ara.config.get('Analysis', 'configfile'))

    print "Checking and joining job(s):", options.job

    errors = False
    for job in ara.jobs:
        if job.name == options.job or options.job == "all":
            if check(job):
                rc = join(job)
                if not rc:
                    errors = True
            else:
                errors = True
                rm(job)

    if errors:
        print "*******************************************************************************"
        print "* There were errors, check script output                                      *" 
        print "*******************************************************************************"

if __name__ == "__main__":
    sys.exit(main())
