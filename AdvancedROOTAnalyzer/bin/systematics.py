#!/bin/env python

import os,sys,optparse
import ara

# define systematics jobs

mydict = eval("""
{ 
  "LooseMuonRelIso" : 
  { 
    "default" : "2.0", 
    "syst" : ( "0.5", "1.0", "1.5", "2.5" )
  },
  "TL_jetpt_min" :
  {
    "default" : "50.",
    "syst" : ( "40.", "60.", "70." )
  },
  "TL_mt_max" :
  {
    "default" : "40.",
    "syst" : ( "30.", "50." )
  }
}
""")

def create_and_submit_jobs(template, replacements, syst, value, basename):
    """Create configuration file "destination" by replacing the words 
    in the dictionary, submit and collect jobs"""
    suffix = "%s_%s_%s" % ( basename, syst, value)
    destination = os.environ["ARASYS"]+"/config/analyzer_%s.cfg"  % ( suffix )
    ara.create_config_file(template, replacements, destination)
    print "Starting systematics", suffix
    os.system("tightloose.py -b " + suffix)

######################################################################
# main
def main():
    usage = """

%prog configfile

Submit systematic studies jobs for the given configuration file. Creates necessary
configuration files and submits and collects the jobs."""
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global ARA configuration file",
                         default=ara.defaultConfigFileName)

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

    # Loop over all systematics
    for syst in mydict:
        print "Starting systematics", syst, "replacing", mydict[syst]["default"], "with", mydict[syst]["syst"]

        repMap = { }
        # Switch off writing tree contents to file, not necessary for systematics,
        # just taking too much disk space
        repMap["FillTree: true"] = "FillTree: false"
        search = syst+': '+mydict[syst]["default"]
        if isinstance(mydict[syst]["syst"], tuple):
            for value in mydict[syst]["syst"]:
                repMap[search] = syst+': '+value
                create_and_submit_jobs(templateFileName, repMap, syst, value, basename)
        else:
            value = mydict[syst]["syst"]
            repMap[search] = syst+': '+value
            create_and_submit_jobs(templateFileName, repMap, syst, value, basename)
        
######################################################################
# entry point if executed
if __name__ == "__main__":
    sys.exit(main())

