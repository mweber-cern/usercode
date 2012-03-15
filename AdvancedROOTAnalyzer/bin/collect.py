#!/bin/env python
import os
import sys
import ConfigParser
import stat
import optparse
import ara

######################################################################
# functions

def check(job, period):
    """Check log files for errors and warnings"""
    global options

    success = True
    print "Checking job", job
    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period

    errors = [ "exception", "segmentation", "error", "err:" ]
    warnings = [ "warning", "wrn" ]
    requirements = [ "End executing" ]
    
    # List all log files
    for filename in os.listdir(myDir):
        if job in filename:
            if "stdout.log" in filename:
                success = ara.check_log(myDir+'/'+filename, errors, warnings, requirements)
                if not success:
                    return False
            if "stderr.log" in filename:
                success = ara.check_log(myDir+'/'+filename, errors, warnings, None)
                if not success:
                    return False
    return True

def join(job, period):
    global options

    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period
    os.chdir(myDir)
    # input file
    try:
        filelist = ara.getCommandOutput2('ls ' + myDir + '/' + job + '_*.root').splitlines()
    except RuntimeError:
        print "ERROR: No result files found"
        return False
    # get file dates
    lastDate = 0
    for inputFile in filelist:
        file_stats = os.stat(inputFile)
        lastDate = max(file_stats[stat.ST_MTIME], file_stats[stat.ST_CTIME], lastDate)
    # output file
    joinFile = myDir + '/' + job + ".root"
    joinDate = 0
    try:
        file_stats = os.stat(joinFile)
        joinDate = max(file_stats[stat.ST_MTIME], file_stats[stat.ST_CTIME], joinDate)
    except OSError:
        pass
    # Do not join if joinFile already newer than input files
    if (joinDate > lastDate and not options.force):
        return True
    # now join
    if isinstance (filelist, str):
        os.system("hadd -f " + job + ".root " + " " + filelist)
    else:
        os.system("hadd -f " + job + ".root " + " ".join(filelist))
    return True

def rm(job, period):
    global options

    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period
    print "Removing output file", job + ".root"
    os.system('rm -f ' + myDir + '/' + job + ".root")

def check_join(items, all=True):
    errors = False
    for (period, jobs) in items:
        for job in jobs.split(','):
            if all or options.job in job:
                if check(job, period):
                    if not join(job, period):
                        errors = True
                        rm(job, period)
                else:
                    errors = True
                    rm(job, period)
    return errors

######################################################################
# main
def main():
    global options

    usage = "usage: %prog [options] selection [job]"
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-f", "--force-overwrite", dest="force",
                         help="Force overwriting of existing .root files",
                         action="store_true",
                         default=False)

    (options, args) = optParser.parse_args()
    if len(args) < 1 or len(args) > 2:
        optParser.print_help()
        return 1

    try:
        ROOTSYS=os.environ['ROOTSYS']
    except KeyError:
        print("You must setup correct ROOTSYS for AdvancedROOTAnalyzer to work")
        return 1

    # setup local arguments
    options.selection = args[0]
    if len(args) == 2:
        options.job = args[1]
    else:
        options.job = "default"

    # read global configuration file
    ara.config.read(options.cfgfile)

    # Check ROOT version
    rootversion = ara.getCommandOutput2("root-config --version");
    first = rootversion.split("/")[0]
    fix   = rootversion.split("/")[1]
    major = int( first.split(".")[0])
    minor = int( first.split(".")[1])
    if major < 5 or (major == 5 and minor < 32):
        print "You are using to old root version", rootversion
        print "You need to use at least ROOT 5.32.00 for hadd to work correctly"
        print "All your histograms filled by text strings might be wrong"

    print "Checking and joining job(s):", options.job

    errors = False
    try:
        # check if this is a job group
        items = ara.config.items('Jobs:'+options.job)
        # OK, is a job group, check all jobs in job group
        errors = check_join(items, all=True)
    except ConfigParser.NoSectionError:
        # OK, is not a job group, take default group
        # and check all jobs which contain options.job inside their name 
        try:
            items = ara.config.items('Jobs:default')
            errors = check_join(items, all=False)
        except ConfigParser.NoSectionError:
            print "No [Jobs:default] section found in configuration file!"
            errors = True

    if errors:
        print "*******************************************************************************"
        print "* There were errors, check script output                                      *" 
        print "*******************************************************************************"

if __name__ == "__main__":
    sys.exit(main())
