#!/bin/env python
######################################################################
# dcache_transfer.py - transfer files from dcache to local directory
# only transfer files which are not there or different in size
#
# (C) RWTH Aachen University, III. Physikalisches Insitut A
# Author: M. Weber

import os
import sys
import optparse
import stat
import ara
import subprocess

def main():
    usage = """usage: %prog [options] dcache-dir

This program transfers files from the given directory on dcache to 
the current local directory. It checks the time stamp and sizes of
all files, and only transfers those files which are newer or not
completely transferred."""
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-f", "--force-deletion", action="store_true", 
                         dest="forcedeletion",
                         help="force deletion of local files",
                         default=False)
    (options, args) = optParser.parse_args()
    ara.config.read(options.cfgfile)

    # user must give path as argument
    if len(args) != 1:
        optParser.print_help()
        return 1

    # get file names from Grid storage element
    print "Reading files from grid storage, please be patient..."
    srcpath = ara.getUserStoragePath(args[0])
    srmfiles = ara.uberftpls(srcpath, True)
    # get local file names
    print "Reading local file names and stats..."
    filelist = os.listdir(".")
    localfiles = [ ]
    for localfile in filelist:
        file_stats = os.stat(localfile)
        if stat.S_ISREG(file_stats.st_mode):
            size = file_stats.st_size
            localfiles.append([size, localfile])

    # loop over all grid files and see if they are here
    print "Reconciling changes..."
    transferfiles = [ ] 
    for srmfile in srmfiles:
        if not srmfile[1].endswith(".root"):
            print "Warning: Found non-ROOT file on grid storage:", srmfile
            continue
        try:
            index = localfiles.index(srmfile)
        except ValueError:
            transferfiles.append(srmfile[1])
            print "Need to transfer file", srmfile[1]

    # determine if local ROOT files need to be deleted
    deletefiles = [ ]
    for localfile in localfiles:
        if not localfile[1].endswith(".root"):
            continue
        try:
            index = srmfiles.index(localfile)
        except ValueError:
            deletefiles.append(localfile[1])
            
    for localfile in deletefiles:
        if options.forcedeletion:
            print "Deleting local file", localfile
            os.unlink(localfile)
        else:
            print "Need to delete", localfile

    if len(deletefiles) > 0 and not options.forcedeletion:
        print "Use -f switch to delete files or remove files manually"
        return 1

    srcpath=ara.config.get('Grid', 'se')+srcpath
    destpath = 'file:///'+os.getcwd()+'/'
    for filename in transferfiles:
        print filename
        p = subprocess.Popen(["srmcp", srcpath+filename, destpath+filename],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if len(stderr) > 0:
            print "Errors found:", stderr
            return 1
        if len(stdout) > 0:
            print "Stdout", stdout

    print "Done."

    return 0

if __name__ == "__main__":
    sys.exit(main())
