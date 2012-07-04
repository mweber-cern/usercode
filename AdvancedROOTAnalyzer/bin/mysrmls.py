#!/bin/env python
######################################################################
# mysrmls.py - list files relative to user grid storage dir
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
    srmfiles = ara.uberftpls(ara.getUserStoragePath(args[0]), True)
    # sum up size
    totalsize = 0
    for (size, filename) in srmfiles:
        totalsize = totalsize + size

    TB = int(totalsize/1000000000000)
    GB = int(totalsize/1000000000) % 1000
    MB = int(totalsize/1000000) % 1000
    KB = int(totalsize/1000) % 1000
    B  = int(totalsize % 1000)
    print "Found", len(srmfiles), "files:", TB, "TB", GB, "GB", MB, "MB", KB, "KB", B, "B"

    return 0

if __name__ == "__main__":
    sys.exit(main())
