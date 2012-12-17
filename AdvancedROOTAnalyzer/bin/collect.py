#!/bin/env python
import os
import sys
import ConfigParser
import stat
import optparse
import ara

######################################################################
# functions

def resubmit(filename, period):
    global options
    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period
    os.chdir(myDir)
    condor_jobfile = myDir+'/'+filename+"_condor.cfg"
    ara.wait_for_jobs(options.njobs)
    print "Resubmitting", condor_jobfile
    rc = ara.getCommandOutput2("condor_submit " + condor_jobfile)
    return rc

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
    
    # Find all job files from directory
    jobfiles = []
    for filename in os.listdir(myDir):
        if filename.startswith(job+'_') and filename.endswith("_condor.cfg"):
            jobfiles.append(filename[:filename.find("_condor.cfg")])

    if len(jobfiles) == 0:
        print 'Error, no job files found for process', job
        return False

    # Check all job files
    good = True
    for filename in jobfiles:
        success = ara.check_log(myDir+'/'+filename+"_stdout.log", 
                                errors, warnings, requirements)
        if not success:
            good = False
            if (options.resubmit):
                resubmit(filename, period)
        else:
            success = ara.check_log(myDir+'/'+filename+"_stderr.log",
                                    errors, warnings, None)
            if not success:
                good = False
                if (options.resubmit):
                    resubmit(filename, period)
    return good

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
    # Change environment for ROOT if necessary (takes some time)
    if options.rootversion_workaround:
        command = "source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh ; source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.01/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh ; hadd -f "
    else:
        command = "hadd -f "
    # now join
    if isinstance (filelist, str):
        os.system(command + job + ".root " + " " + filelist)
    else:
        os.system(command + job + ".root " + " ".join(filelist))
    return True

def rm(job, period):
    global options

    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period
    print "Removing output file", job + ".root"
    os.system('rm -f ' + myDir + '/' + job + ".root")

def check_join(jobgroup, job, period):
    if not check(job,period):
        rm(job, period)
        return False
    if not join(job, period):
        rm(job, period)
        return False
    return True

######################################################################
# main
def main():
    global options

    usage = """usage: %prog [options] selection [jobgroup] [pattern]

Will check output of CONDOR jobs and join files in given jobgroup matching the
specified pattern. "jobgroup" defaults to the default job group named
"default", and "pattern" defaults to "all", i.e. all jobs in the given
jobgroup will be checked and joined."""
    optParser = optparse.OptionParser(usage)
    defaultnjobs=20
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-f", "--force-overwrite", dest="force",
                         help="Force overwriting of existing .root files",
                         action="store_true",
                         default=False)
    optParser.add_option("-r", "--resubmit", dest="resubmit", action="store_true",
                         help="Resubmit failed jobs")
    optParser.add_option("-n", "--njobs", dest="njobs",
                         help="how many jobs to run at once if resubmitting",
                         default=defaultnjobs)

    (options, args) = optParser.parse_args()
    if len(args) < 1 or len(args) > 3:
        optParser.print_help()
        return 1

    try:
        ROOTSYS=os.environ['ROOTSYS']
    except KeyError:
        print("You must setup correct ROOTSYS for AdvancedROOTAnalyzer to work")
        return 1

    # get number of jobs
    try:
        options.njobs = int(options.njobs)
    except ValueError:
        print "Number of jobs has to be a number"
        return 1

    if options.njobs <= 0:
        print "Number of jobs must be greater than zero"
        return 2

    # read global configuration file
    ara.config.read(options.cfgfile)

    # setup local arguments
    options.selection = args[0]
    options.jobgroup = 'default'
    options.pattern = 'all'
    if len(args) >= 2:
        options.jobgroup = args[1]
    if len(args) >= 3:
        options.pattern = args[2]

    # Check ROOT version
    options.rootversion_workaround = False
    rootversion = ara.getCommandOutput2("root-config --version").strip()
    first = rootversion.split("/")[0]
    fix   = rootversion.split("/")[1]
    major = int( first.split(".")[0])
    minor = int( first.split(".")[1])
    if major < 5 or (major == 5 and minor < 32):
        print "You are using an old root version", rootversion
        print "Enabling workaround for calling hadd"
        options.rootversion_workaround = True

    errors = False
    try:
        # check if this is a job group
        section = 'jobgroup:'+options.jobgroup
        items = ara.config.items(section)
        # OK, is a job group
        print "Checking and joining jobs from job group", options.jobgroup
        for (period, jobs) in items:
            print "Job group period", period, "consists of jobs", jobs
            for job in jobs.split(','):
                if options.pattern in job or options.pattern == 'all':
                    print "Collecting job", job
                    errors = check_join(options.jobgroup, job, period) and errors
    except ConfigParser.NoSectionError:
        print section, "is not a job group."
        print "You must configure job groups in file %s" % ( options.cfgfile )
        return 1

    if errors:
        print "*******************************************************************************"
        print "* There were errors, check script output                                      *" 
        print "*******************************************************************************"
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
