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
    srmfiles = ara.srmls(args[0], True)
    filenames=[]
    for (size, filename) in srmfiles:
        filenames.append(filename)

    # remove individual files
    if len(filenames) > 0:
        print "Removing files"
        msg = ara.srmrm(args[0], filenames)
        if msg != '':
            print msg

    # remove directory
    print "Removing directory", args[0]
    msg = ara.srmrmdir(args[0])
    if msg != '':
        print msg
    return 0

if __name__ == "__main__":
    sys.exit(main())
