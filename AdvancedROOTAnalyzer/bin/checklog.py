#!/bin/env python
######################################################################
# check_log.py - check log files for common errors and warnings
# and validate program finished successfully

import os,sys
import optparse
import ac3ana

def main():
    optParser = optparse.OptionParser()
    defaultconfig=os.path.expanduser("~/.ac3ana")
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=defaultconfig)
    (options, args) = optParser.parse_args()
    if len(args) < 1:
        optParser.print_help()
        return 1

    ac3ana.config.read(options.cfgfile)

    # check log files on command line
    for arg in args:
        if not ac3ana.check_log(arg):
            return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
