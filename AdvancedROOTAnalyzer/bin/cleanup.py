#!/bin/env python
import os
import sys
import ConfigParser
import stat
import optparse
import ara

######################################################################
# functions

def rm(job, period):
    global options

    basedir = ara.config.get('Analysis', 'basedir')
    myDir = basedir + '/' + options.selection + '/' + period
    print "Removing output file", job + ".root"
    os.system('rm -f ' + myDir + '/' + job + ".root")

def cleanup(items, all=True):
    for (period, jobs) in items:
        for job in jobs.split(','):
            if all or options.job in job:
                rm(job, period)

######################################################################
# main
def main():
    global options

    usage = "usage: %prog [options] selection [job]"
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)

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
    print rootversion
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
        errors = cleanup(items, all=True)
    except ConfigParser.NoSectionError:
        # OK, is not a job group, take default group
        # and check all jobs which contain options.job inside their name 
        try:
            items = ara.config.items('Jobs:default')
            errors = cleanup(items, all=False)
        except ConfigParser.NoSectionError:
            print "No [Jobs:default] section found in configuration file!"
            errors = True

    if errors:
        print "*******************************************************************************"
        print "* There were errors, check script output                                      *" 
        print "*******************************************************************************"

if __name__ == "__main__":
    sys.exit(main())
