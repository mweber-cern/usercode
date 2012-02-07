#!/bin/env python
######################################################################
# check_log.py - check log files for common errors and warnings
# and validate program finished successfully

import os,sys
import optparse
import ara

def main():
    optParser = optparse.OptionParser()
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    (options, args) = optParser.parse_args()
    if len(args) < 1:
        optParser.print_help()
        return 1

    ara.config.read(options.cfgfile)

    errors = [ "exception", "segmentation", "error", "err:" ]
    warnings = [ "warning", "wrn" ]
    requirements = [ "End executing" ]

    # check log files on command line
    for arg in args:
        if not ara.check_log(arg, errors, warnings, requirements):
            return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
