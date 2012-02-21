#!/bin/env python

import os,sys,optparse
import ara

def merge(filelist, counter):
    outfilename="merged_" + str(counter) + ".root"
    print "Merging files", filelist, "in file", outfilename
    os.system("hadd " + outfilename + " " + " ".join(filelist))

def main():
    usage = "usage: %prog [options] number_of_files"
    optParser = optparse.OptionParser(usage)
    defaultdir="."
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-d", "--dir", dest="dir",
                         help="directory with ROOT files to merge",
                         default=defaultdir)
    (options, args) = optParser.parse_args()
    ara.config.read(options.cfgfile)

    # user must give path as argument
    if len(args) != 1:
        print "You must give the nunmber of files to merge as argument, use -h for help"
        return 1

    try:
        nmerge = int(args[0])
    except ValueError:
        print "You must give the nunmber of files to merge as argument, use -h for help"
        return 1
        
    if (nmerge <= 1):
        print "Argument must be a positive number > 1"
        return 2

    # get all ROOT files in directory
    filelist = os.listdir(os.path.abspath(options.dir))
    localfiles = [ ]
    for file in filelist:
        if file.endswith(".root"):
            localfiles.append(file);

    nfiles = len(localfiles)
    print "Found", nfiles , "ROOT files."
    counter = 0
    while counter*nmerge < nfiles:
        start = counter*nmerge
        end = min((counter+1)*nmerge, nfiles)
        merge(localfiles[start:end], counter)
        counter = counter + 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
