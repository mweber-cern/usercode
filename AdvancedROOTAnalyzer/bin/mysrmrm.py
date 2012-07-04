#!/bin/env python
######################################################################
# mysrmrm.py - remove files relative to user grid storage dir
#
# (C) RWTH Aachen University III. Physikalisches Insitut A
# Author: M. Weber

import os,sys
import optparse
import ara

def main():
    usage = "usage: %prog [options] dcache-dir-relative-to-user-dir"
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    (options, args) = optParser.parse_args()
    if len(args) != 1:
        optParser.print_help()
        return 1

    ara.config.read(options.cfgfile)

    # get file names from Grid storage element
    path = ara.getUserStoragePath(args[0])
    srmfiles = ara.srmls(path, True)
    filenames=[]
    for (size, filename) in srmfiles:
        filenames.append(filename)

    # remove individual files
    if len(filenames) > 0:
        print "Removing files"
        msg = ara.srmrm(path, filenames)
        if msg != '':
            print msg

    # remove directory
    print "Removing directory", path
    msg = ara.srmrmdir(path)
    if msg != '':
        print msg
    return 0

if __name__ == "__main__":
    sys.exit(main())
