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
import tempfile

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

def getUserStoragePath(directory):
    global config
    user = config.get('Grid', 'user')
    sp = config.get('Grid', 'sp')
    path = sp+'/'+user+'/'+directory+'/'
    path = path.replace('//', '/')
    return path

######################################################################
# list the contents of the given directory on the grid storage element
# returns a list of [filesize, filename]
# depends on global configuration object
def srmls(path, verbose=False):
    offset = 0
    files = []
    foundall = False
    fullpath = config.get('Grid', 'se')+path
    while not foundall:
        if verbose:
            print 'Listing', fullpath, 'at offset', offset
        # call program and record both stderr and stdout
        p = subprocess.Popen(["srmls", "-offset="+str(offset), "-count=999", fullpath],
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

def uberftpls(path, verbose=False):
    ftp = config.get('Grid', 'ftp')
    offset = 0
    fullpath = ftp+path
    files = []
    if verbose:
        print 'Listing ', fullpath
    # create temporary file with commands for uberftp
    commands = """cd %s
ls
""" % ( path )
    # call program and record both stderr and stdout
    p = subprocess.Popen(["uberftp", config.get('Grid','ftp')], 
                         stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(commands)
    if len(stderr) > 0:
        print "Errors found:", stderr
        sys.exit(1)
    for line in stdout.splitlines():
        # remove uberftp prompt
        line = line.replace("uberftp> ", "")
        # split into file size and file name
        try:
            size = int(line.split()[4])
            fname = line.split()[8]
        except ValueError:
            continue
        except IndexError:
            continue
        if verbose:
            print size, fname
        files.append([size, fname])
    return files

def srmrm(path, filenames):
    se = config.get('Grid', 'se')
    command = 'srmrm '
    for filename in filenames:
        print filename
        command = command + " "  + se + path + filename
    msg = getCommandOutput2(command)
    return msg

def srmrmdir(path):
    se = config.get('Grid', 'se')
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
        raise Exception("Error executing grep \"query\" % for file %s " % (searchfor, filename))
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

def get_maximum_jobs(njobs):
    """Throttle the maximum number of jobs during typical office hours in order
not to overload the network"""
    ltime = time.localtime()
    if ltime.tm_hour > 9 and ltime.tm_hour < 20 and ltime.tm_wday < 5:
        return njobs
    else:
        return njobs

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
                if line.split()[5] != "H":
                    jobcount += 1
        # if the maximal job count is exceeded, sleep for a while
        if jobcount > njobs:
            print "Exceeded " + str(njobs) + " jobs, sleeping for 30 s..."
            time.sleep(30)
        else:
            wait = False


def create_config_file(template, replacements, destination):
    """Open the template file and write its content to destination, 
    replacing all keys in the dictionary 'replacements'"""
    opwd = os.getcwd()
    os.chdir("%s/config" % os.environ["ARASYS"])
    with open(destination, 'w') as f2:
        with open(template) as f:
            for line in f:
                for key in replacements:
                    line = line.replace(key, replacements[key])
                f2.write(line)
    os.chdir(opwd)

