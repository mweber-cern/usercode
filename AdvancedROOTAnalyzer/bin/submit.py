#!/bin/env python
import os
import sys
import ConfigParser
import optparse
import ara
import shutil
import time

######################################################################
# misc helper functions
def getType(job):
    if "DATA" in job.upper():
        return "data"
    elif "SIGNAL" in job.upper():
        return "signal"
    elif "BACKGROUND" in job.upper():
        return "background"
    else:
        return "mc"

######################################################################
# submitters

def submit_condor_job(executable, arguments, inbox, outbox, jobname):
    global options
    scriptdir = os.environ['ARASYS'] + '/config'
    condor_template_file = open(scriptdir + '/' + options.template)
    condor_template = condor_template_file.read()
    stderr = "_".join((jobname, "stderr.log"))
    stdout = "_".join((jobname, "stdout.log"))
    cfgFile = "_".join((jobname, "condor.cfg"))
    repMap = {}
    repMap["executable"] = executable
    repMap["arguments"] = arguments
    repMap["inputfiles"] = inbox
    repMap["stdout"] = stdout
    repMap["stderr"] = stderr
    repMap["log"] = "condor.log"
    repMap["outputfiles"] = ",".join((stderr, stdout, outbox))
    content = condor_template % repMap
    jobfile = open(cfgFile, "w")
    jobfile.write(content)
    jobfile.close()
    ara.wait_for_jobs(options.njobs)
    rc = ara.getCommandOutput2("condor_submit " + cfgFile)
    return rc

def submit_condor_with_wrapper(executable, arguments, inbox, outbox, jobname):
    scriptdir = os.environ['ARASYS'] + '/bin'
#    job_starter_file = open(scriptdir + "/job_starter_template.sh")
#    job_starter = job_starter_file.read()
#    repMap = {}
#    content = job_starter % repMap
#    job_starter_filename = "job_starter.sh"
#    exefile = open(job_starter_filename, "w")
#    exefile.write(content)
#    exefile.close()
#    os.system("chmod 755 " + job_starter_file)
    job_starter_filename = scriptdir + "/job_starter.sh";
    submit_condor_job(job_starter_filename, 
                      " ".join((executable, arguments)), 
                      inbox, outbox, jobname)

def submit(jobgroup, job, period):
    global options

    # read files from directory
    # call program and record both stderr and stdout
    filelist = []
    section = jobgroup+':'+period
    try:
        for thisdir in ara.config.get(section, job).split(','):
            print 'Adding files from directory', thisdir
            # run over files on dCache?
            if thisdir.startswith("dcap:/"):
                # list files on dcache
                path=thisdir[5:]
                files = ara.uberftpls(path)
                for (size, filename) in files:
                    filelist.append(ara.config.get('Grid', 'dcap') + path + '/' + filename)
            else:
                filespec =  os.path.abspath(os.path.expanduser(thisdir)) + '/*.root'
                files = ara.getCommandOutput2('ls ' + filespec)
                for filename in files.splitlines():
                    filelist.append(filename)

    except ConfigParser.NoSectionError:
        print "Could not find section %s in main configuration file %s" % ( section , options.cfgfile )
        print "You must create this section and add entries to it"
        return False

    # determine file splitting
    print 'Found', len(filelist), 'files'
    list_of_lists = ara.split_in_jobs(filelist, options.nsplit)
    n = 0
    for files in list_of_lists:
        print "Job #", n, ": ", len(files), " files"

        # read template
        cfgPath = os.environ['ARASYS'] + '/config'
        cfgName = 'analyzer_' + options.analysiscfg + '.cfg'
        cfgTemplateFile = open(cfgPath + '/' + cfgName)
        cfg = cfgTemplateFile.read()
        repMap = {}
        repMap["sample"] = job
        repMap["inputfile"] = " ".join(files)
        outputFile = job + '_' + str(n) + '.root'
        repMap["outputfile"] = outputFile
        repMap["type"] = getType(job)
        repMap["configPath"] = cfgPath
        content = cfg % repMap

        # create a job configuration file
        basedir = ara.config.get('Analysis', 'basedir')
        # create directory
        myDir = basedir + '/' + options.selection + '/' + period
        os.system("mkdir -p " + myDir)
        jobName = job + '_' + str(n)
        cfgFileNameBase = job + '_' + str(n) + '.cfg'
        cfgFileName = myDir + '/' + cfgFileNameBase
        cfgFile = open(cfgFileName, "w")
        cfgFile.write(content)
        cfgFile.close()

        executable = os.environ['ARASYS'] + '/bin/analyzer_' + options.selection
        arguments = cfgFileNameBase
        # need to go into this directory - submit job from there
        os.chdir(myDir)
        inbox = cfgFileName #",".join((cfgFileName, "Weight3D.root"))
        outbox = outputFile # ",".join((outputFile, "Weight3D.root"))
        submit_condor_with_wrapper(executable, cfgFileNameBase,
                                   inbox, outbox, jobName)

        time.sleep(3)

        # Next file
        n += 1

    return True

######################################################################
# main
def main():
    usage = """usage: %prog [options] selection [jobgroup] [pattern]

Will create a and submit jobs from the given jobgroup matching the pattern to
CONDOR. The job files and output will be stored in BASEDIR/selection in
subdirectories, one for each configured period (for configuration see default
configuration file). "jobgroup" defaults to the default job group named
"default", and "pattern" defaults to "all", i.e. all jobs in the given jobgroup
will be submitted."""
    optParser = optparse.OptionParser(usage)
    defaulttemplate="condor_template.cfg"
    defaultnjobs=20
    defaultsplit=1000
    defaultanalysiscfg='default'
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global ARA configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-t", "--template", dest="template",
                         help="condor template file",
                         default=defaulttemplate)
    optParser.add_option("-n", "--njobs", dest="njobs",
                         help="how many jobs to run at once",
                         default=defaultnjobs)
    optParser.add_option("-s", "--split", dest="nsplit",
                         help="split each process into this number of jobs",
                         default=defaultsplit)
    optParser.add_option("-a", "--analysis-config", dest="analysiscfg",
                         help="analysis configuration file (default: '"+defaultanalysiscfg+"')",
                         default=defaultanalysiscfg)
                         

    global options
    (options, args) = optParser.parse_args()
    if len(args) < 1 or len(args) > 3:
        optParser.print_help()
        return 1

    if os.environ['CMSSW_BASE'] == '':
        raise "You must setup correct CMSSW version for AdvancedROOTAnalyzer to work"

    # get number of splits
    try:
        options.nsplit = int(options.nsplit)
    except ValueError:
        print "Split has to be a number"
        return 1

    if options.nsplit <= 0:
        print "Split must be greater than zero"
        return 2

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

    # copy executable in place
    src = os.environ['ARASYS'] + '/analyzer/analyzer'
    dest = os.environ['ARASYS'] + '/bin/analyzer_' + options.selection
    shutil.copy(src, dest)

    errors = False
    try:
        # check if this is a job group
        section = 'jobgroup:'+options.jobgroup
        items = ara.config.items(section)
        # OK, is a job group
        print "Starting job group", options.jobgroup
        for (period, jobs) in items:
            print "Job group period", period, "consists of jobs", jobs
            for job in jobs.split(','):
                if options.pattern in job or options.pattern == 'all':
                    print "Submitting job", job
                    errors = submit(options.jobgroup, job, period) and errors

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
