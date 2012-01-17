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
    optParser = optparse.OptionParser()
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.configFileName)
    (options, args) = optParser.parse_args()
    if len(args) != 1:
        optParser.print_help()
        return 1

    ara.config.read(options.cfgfile)

    # get file names from Grid storage element
    srmfiles = ara.srmls(args[0], True)
    print "Found", len(srmfiles), "files"

    return 0

if __name__ == "__main__":
    sys.exit(main())
