#############################################################################
# ara.py - often used functions etc. for CMSSW, grid and cluster jobs etc.
#
# (C) RWTH Aachen University III. Physikalisches Insitut A
# Author: M. Weber

import os
import sys
import ConfigParser
import subprocess
import getpass
import time

######################################################################
# global variables keeping options and configuration
global defaultConfigFileName
global config
global options

######################################################################
# setup

# sanity check for environment setting
if os.environ["ARASYS"] == '':
    raise "You must have setup your Advanced ROOT Analyzer by sourcing bin/setup.sh"
defaultConfigFileName = os.environ['ARASYS'] + ("/config/ara.cfg")
config=ConfigParser.ConfigParser()

######################################################################
# Execute a command and get its output
def getCommandOutput2(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError, '%s failed with exit code %d' % (command, err)
    return data

def getStoragePath(directory):
    global config
    user = config.get('Grid', 'user')
    se = config.get('Grid', 'se')
    sp = config.get('Grid', 'sp')
    path = sp+'/'+user+'/'+directory+'/'
    path = path.replace('//', '/')
    return (se, path)

######################################################################
# list the contents of the given directory on the grid storage element
# returns a list of [filesize, filename]
# depends on global configuration object
def srmls(directory, verbose=False):
    (se, path) = getStoragePath(directory)
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

def srmrm(directory, filenames):
    (se, path) = getStoragePath(directory)
    command = 'srmrm ' 
    for filename in filenames:
        print filename
        command = command + " "  + se + path + filename
    msg = getCommandOutput2(command)
    return msg

def srmrmdir(directory):
    (se, path) = getStoragePath(directory)
    return getCommandOutput2('srmrmdir '+se+path)

def split_in_jobs(filenames, njobs):
    """divide given file list into a predefined number of equally sized chunks"""   
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
    grep_file(logfilename, warnings)

    # check for requirements
    if requirements != None:
        if not grep_file(logfilename, requirements, verbose=False):
            print "Requirements not found in logfile %s: " % ( logfilename )
            return False

    # inform user
    return True

def wait_for_jobs(njobs):
    wait = True
    while wait:
        # check how many jobs are running already
        ntrial = 0
        while (ntrial < 3):
            try:
                joblist = getCommandOutput2("condor_q");
                break
            except RuntimeError as message:
                print message
                print "Waiting for condor to settle..."
                time.sleep(60)
                ntrial += 1
        jobcount = 0
        for line in joblist.splitlines():
            if getpass.getuser() in line:
                jobcount += 1
        # if the maximal job count is exceeded, sleep for a while
        if jobcount >= njobs:
            print "Reached " + str(njobs) + " jobs, sleeping for 30 s..."
            time.sleep(30)
        else:
            wait = False
