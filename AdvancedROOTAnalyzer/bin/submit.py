#!/bin/env python
import os
import sys
import ConfigParser
import optparse
import ara
import getpass
import time
from ROOT import TEnv

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
    wait = True
    while wait:
        # check how many jobs are running already
        joblist = ara.getCommandOutput2("condor_q");
        jobcount = 0
        for line in joblist.splitlines():
            if getpass.getuser() in line:
                jobcount += 1
        # if the maximal job count is exceeded, sleep for a while
        if jobcount >= options.njobs:
            print "Reached " + str(options.njobs) + " jobs, sleeping for 30 s..."
            time.sleep(30)
        else:
            wait = False
    os.system("condor_submit " + cfgFile)

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

def submit(job):
    global options
    global myEnv

    # read files from directory
    # call program and record both stderr and stdout
#    filelist = subprocess.check_output(['ls', os.path.abspath(job.directory + '/*.root')])
    filespecs = []
    for thisdir in job.directory:
        filespec =  os.path.abspath(os.path.expanduser(thisdir)) + '/*.root'
        filespecs.append(filespec)
    filelist = ara.getCommandOutput2('ls ' + " ".join(filespecs))

    # determine file splitting
    list_of_lists = ara.split_in_jobs(filelist.splitlines(), options.nsplit)
    n = 0
    for files in list_of_lists:
        print "Job #", n, ": ", len(files), " files"

        # read template
        cfgPath = os.environ['ARASYS'] + '/config'
        cfgName = 'analyzer_' + options.selection + '.cfg'
        cfgTemplateFile = open(cfgPath + '/' + cfgName)
        cfg = cfgTemplateFile.read()
        repMap = {}
        repMap["inputfile"] = " ".join(files)
        outputFile = job.name + '_' + str(n) + '.root'
        repMap["outputfile"] = outputFile
        repMap["type"] = job.filetype
        content = cfg % repMap

        # create a job configuration file
        basedir = myEnv.GetValue("basedir", ".");
        myDir = basedir + '/' + options.selection + '/' + job.period
        jobName = job.name + '_' + str(n)
        cfgFileNameBase = job.name + '_' + str(n) + '.cfg'
        cfgFileName = myDir + '/' + cfgFileNameBase
        cfgFile = open(cfgFileName, "w")
        cfgFile.write(content)
        cfgFile.close()

        executable = os.environ['ARASYS'] + '/analyzer/ara'
        arguments = cfgFileNameBase
        # need to go into this directory - submit job from there
        os.chdir(myDir)
        inbox = cfgFileName
        outbox = outputFile # ",".join((outputFile, "Weight3D.root"))
        submit_condor_with_wrapper(executable, cfgFileNameBase,
                                   inbox, outbox, jobName)

        # Next file
        n += 1

def make_dirs():
    global options
    global myEnv

    for i in xrange(myEnv.GetValue("N_period", 1)):
        os.system("mkdir -p " + myEnv.GetValue("basedir", ".") + '/'
                  + options.selection + '/'
                  + myEnv.GetValue("period." + str(i) + ".name", "."))


def main():
    usage = "usage: %prog [options] selection job | all"
    optParser = optparse.OptionParser(usage)
    defaulttemplate="condor_template.cfg"
    defaultnjobs=20
    defaultsplit=100
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.configFileName)
    optParser.add_option("-t", "--template", dest="template",
                         help="condor template file",
                         default=defaulttemplate)
    optParser.add_option("-n", "--njobs", dest="njobs",
                         help="how many jobs to run at once",
                         default=defaultnjobs)
    optParser.add_option("-s", "--split", dest="nsplit",
                         help="split each process into this number of jobs",
                         default=defaultsplit)

    global options
    (options, args) = optParser.parse_args()
    if len(args) != 2:
        optParser.print_help()
        return 1

    # get number of splits
    try:
        nsplit = int(options.nsplit)
    except ValueError:
        print "Split has to be a number"
        return 1

    if nsplit <= 0:
        print "Split must be greater than zero"
        return 2

    # get number of jobs
    try:
        njobs = int(options.njobs)
    except ValueError:
        print "Number of jobs has to be a number"
        return 1

    if njobs <= 0:
        print "Number of jobs must be greater than zero"
        return 2

    options.selection = args[0]
    options.job = args[1]
    options.njobs = njobs
    options.nsplit = nsplit

    # read global configuration file
    ara.config.read(options.cfgfile)

    # read ROOT configuration file
    configPath = os.environ['ARASYS'] + '/config'
    global myEnv
    myEnv = TEnv(configPath + '/' + ara.config.get('Analysis', 'configfile'))

    # create directory if necessary
    make_dirs()

    print "Starting job(s):", options.job

    for job in ara.jobs:
        if job.name == options.job or options.job == "all":
            submit(job)

if __name__ == "__main__":
    sys.exit(main())
