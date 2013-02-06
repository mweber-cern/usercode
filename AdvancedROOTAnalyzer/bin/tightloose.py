#!/bin/env python

import os,sys,optparse
import ara

######################################################################
# main
def main():
    usage = """

%prog [-b] [-t] configfile

Batch submission of default job, creation of FakeRate file, submission of fake rate jobs."""
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global ARA configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-t", "--tight-lose-ratio", action="store_true",
                         dest="tlratio",
                         help="compute T/L ratio")
    optParser.add_option("-b", "--submit-base", action="store_true",
                         dest="submitbase",
                         help="submit the base version, too, not only single+double fakes (implies -t)")


    global options
    (options, args) = optParser.parse_args()
    if len(args) < 1 or len(args) > 1:
        optParser.print_help()
        return 1

    basename = args[0]
    basepath = os.environ["ARASYS"]+"/config/analyzer_"+basename
    templateFileName = basepath+".cfg"

    # Check if template file exists
    if not os.path.isfile(templateFileName):
        print "File "+templateFileName+" not found. Exiting."
        return 2

    # Check if base job needs to be submitted
    if (options.submitbase):
        os.system("submit.py -a %s" % ( basename ))
        # wait until all jobs are terminated
        ara.wait_for_jobs(0)
        # check and resubmit if necessary
        os.system("collect.py -r %s " % (basename))
        # wait until all jobs are terminated
        ara.wait_for_jobs(0)

    if (options.tlratio or options.submitbase):
        # create fakerate file
        opwd = os.getcwd()
        os.chdir(os.environ["ARASYS"]+"/root")
        os.system("root -q -b 'tlratio.C(\"%s\")'" % basename)
        os.chdir(opwd)

    # Create descendent configuration file (singlefake)
    repMap = { }
    repMap["AnalysisType: default"] = "AnalysisType: singlefake"
    repMap["FakeRateFile: FakeRate.root"] = "FakeRateFile: FakeRate_%s.root" % ( basename )
    destination = basepath+"_singlefake.cfg"
    ara.create_config_file(templateFileName, repMap, destination)

    # submit
    os.system("submit.py -a %s_singlefake" % basename)

    # Create descendent configuration file (doublefake)
    repMap["AnalysisType: default"] = "AnalysisType: doublefake"
    repMap["FakeRateFile: FakeRate.root"] = "FakeRateFile: FakeRate_%s.root" % ( basename )
    destination = basepath+"_doublefake.cfg"
    ara.create_config_file(templateFileName, repMap, destination)

    # submit
    os.system("submit.py -a %s_doublefake" % basename)

    # wait for jobs to finish and collect jobs
    ara.wait_for_jobs(0)
    os.system("collect.py -r %s_singlefake" % basename)
    os.system("collect.py -r %s_doublefake" % basename)

######################################################################
# entry point if executed
if __name__ == "__main__":
    sys.exit(main())

