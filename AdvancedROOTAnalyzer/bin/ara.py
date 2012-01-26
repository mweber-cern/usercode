#############################################################################
# ara.py - often used functions etc. for CMSSW, grid and cluster jobs etc.
#
# (C) RWTH Aachen University III. Physikalisches Insitut A
# Author: M. Weber

import os
import sys
import ConfigParser
import subprocess

######################################################################
# sanity check for environment setting
if os.environ["ARASYS"] == '':
    raise "You must have setup your Advanced ROOT Analyzer by sourcing bin/setup.sh"

######################################################################
# global variables keeping options and configuration
global configFileName
configFileName = os.environ['ARASYS'] + ("/config/ara.cfg")
global config

config=ConfigParser.ConfigParser()
global options

######################################################################
# Process configuration

class Process:
    def __init__(self, name, directory, filetype, period):
        self.name = name
        self.directory = directory
        self.filetype = filetype
        self.period = period

global jobs
jobs = [ 
    Process("data",
            ["~/SusyWG/CMSSW428_v77/data/DoubleMu_Run2011A-03Oct2011-v1_AOD",
             "~/SusyWG/CMSSW428_v77/data/DoubleMu_Run2011A-05Aug2011-v1_AOD",
             "~/SusyWG/CMSSW428_v77/data/DoubleMu_Run2011A-May10ReReco-v1",
             "~/SusyWG/CMSSW428_v77/data/DoubleMu_Run2011A-PromptReco-v4_AOD",
             "~/SusyWG/CMSSW428_v77/data/DoubleMu_Run2011B-PromptReco-v1_AOD"],
            "data", "2011"),
    Process("signal_LM1",
            ["~/SusyWG/CMSSW428_v77/RPV-signal/signal_LM1_lp211eq0.01"],
            "signal", "2011"),
    Process("background_LM1",
            ["~/SusyWG/CMSSW428_v77/RPV-signal/signal_LM1_lp211eq0.01"],
            "background", "2011"),
    Process("susy_LM1",
            ["~/SusyWG/CMSSW428_v77/RPV-signal/signal_LM1_lp211eq0.01"],
            "mc", "2011"),
#    Process("tt", 
#            "[/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v2/AODSIM"],
#           "mc", "2011")
    Process("tt2l2n",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_TTTo2L2Nu2B_7TeV-powheg-pythia6_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("dyllpf",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1_AODSIM_Reskim"],
            "mc", "2011"),
#    Process("dyll",
#            ["~/SusyWG/CMSSW428_v77/mc/Fall11_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph_v2"],
#            "mc", "2011"),
    Process("stbars", 
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_Tbar_TuneZ2_s-channel_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("stbart",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_Tbar_TuneZ2_t-channel_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("stbartw",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("sts", 
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_T_TuneZ2_s-channel_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("stt",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_T_TuneZ2_t-channel_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"], 
            "mc", "2011"),
    Process("sttw",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola_Fall11-PU_S6_START42_V14B-v1"],
            "mc", "2011"),
    Process("wjetstolnu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_v2"],
            "mc", "2011"),
    Process("wwjetstolnu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph"],
            "mc", "2011"),
    Process("wzjetsto2l2q",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_WZJetsTo2L2Q_TuneZ2_7TeV-madgraph"],
            "mc", "2011"),
    Process("wzjetsto3lnu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_WZJetsTo3LNu_TuneZ2_7TeV-madgraph"],
            "mc", "2011"),
    Process("zzjetsto2l2nu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_ZZJetsTo2L2Nu_TuneZ2_7TeV-madgraph"],
            "mc", "2011"),
    Process("zzjetsto2l2q",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_ZZJetsTo2L2Q_TuneZ2_7TeV-madgraph"],
            "mc", "2011"),
    Process("zzjetsto4l",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_ZZJetsTo4L_TuneZ2_7TeV-madgraph_v1"],
            "mc", "2011"),
    Process("qcd1520mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-15to20_Mu"],
            "mc", "2011"),
    Process("qcd3050mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-30to50_Mu"],
            "mc", "2011"),
    Process("qcd5080mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-50to80_Mu"],
            "mc", "2011"),
    Process("qcd80120mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-80to120_Mu"],
            "mc", "2011"),
    Process("qcd120150mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-120to150_Mu"],
            "mc", "2011"),
    Process("qcd150mu",
            ["~/SusyWG/CMSSW428_v77/mc/Fall11_QCD_Pt-150_MuPt5Enriched"],
            "mc", "2011")
]

######################################################################
# read options and config
def readOptions():
    return

######################################################################
# Execute a command and get its output
def getCommandOutput2(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError, '%s failed with exit code %d' % (command, err)
    return data

######################################################################
# list the contents of the given directory on the grid storage element
# returns a list of [filesize, filename]
# depends on global configuration object
def srmls(dir, verbose=False):
    global config
    user = config.get('Grid', 'user')
    se = config.get('Grid', 'se')
    sp = config.get('Grid', 'sp')
    path = sp+'/'+user+'/'+dir+'/'
    path = path.replace('//', '/')
    offset = 0
    files = []
    foundall = False
    while not foundall:
        if verbose:
            print 'Listing ', se+path, 'at offset', offset
        # call program and record both stderr and stdout
        p = subprocess.Popen(["srmls", "-offset="+str(offset), "-count=999", se+path],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if len(stderr) > 0:
            print "Errors found:", stderr
            sys.exit(1)
        for line in stdout.splitlines():
            # ignore warnings
            if 'WARNING' in line:
                continue
            # split into file size and file name
            try:
                (size, fname) = line.replace(path, '').split()
                size = int(size)
            except ValueError:
                continue
            # remove dir entry
            if fname == '/':
                continue
            if verbose:
                print size, fname
                files.append([size, fname])
        if len(files) == 999:
            offset = offset + 999;
        else:
            foundall = True
    return files

######################################################################
# divide given file list into a predefined number of equally sized chunks
# 
def split_in_jobs(filenames, njobs):
    # catch the case where a single file name is given, not a list
    list_of_lists = [ ]
    if isinstance(filenames, str):
        files = [ ]
        files.append(filenames)
        list_of_lists.append(files)
        return list_of_lists
    # OK, it is a list, determine splitting
    nfiles = len(filenames)
    files_per_job = nfiles / njobs
    if nfiles % njobs != 0:
        files_per_job += 1

    number_of_jobs = nfiles / files_per_job
    if nfiles % files_per_job != 0:
        number_of_jobs += 1

    # produce lists
    n = 0
    m = 0
    while n*files_per_job < nfiles:
        names = [ ]
        while m < files_per_job and n*files_per_job+m < nfiles:
            names.append(filenames[n*files_per_job+m])
            m += 1
        list_of_lists.append(names)
        m = 0
        n += 1
    return list_of_lists

def grep_file(filename, tokens, verbose=True):
    found = False
    searchfor = "|".join(tokens)
    child = os.popen('/bin/bash -c \'grep -Ei \"' + searchfor + '\" ' + filename + '\'')
    data = child.read().splitlines()
    rc = child.close()
    if len(data) > 0:
        found = True
        if verbose:
            print "Tokens matching %s found in file %s: " % ( searchfor, filename )
            for line in data:
                print line
    if rc == 2:
        raise Exception("Error executing grep \"query\" % for file %s " % (searchfor, filenam))
    return found

def check_log(logfilename, errors, warnings, requirements):
    """Check the log files and return True on success.
       errors is an iterable object containing strings marking errors in log file
       warnings is an iterable object containing strings marking warnings in log file
       requirements is an iterable object containing strings the log file that must appear

       If errors are found, print errors and return False
       If warnings are found, print warnings and return True
       If requirements are not found, print message and return False
    """

    # check for errors
    if grep_file(logfilename, errors):
        return False

    # check for warnings
    grep_file(logfilename, errors)

    # check for requirements
    if requirements != None:
        if not grep_file(logfilename, requirements, verbose=False):
            print "Requirements not found in logfile %s: " % ( logfilename )
            return False

    # inform user
    return True
