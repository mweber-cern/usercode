#!/bin/env python
######################################################################
# gencutflow.py - generate cut flow histogram code
#
# (C) RWTH Aachen University III. Physikalisches Insitut A
# Author: M. Weber

import os,sys
import optparse
import ara

cutFlowTemplate = """
#include "TH1D.h"
#include <iostream>

const char * stages[] = { 
%s
};

TH1D * get_cutflow_histogram()
{
    UInt_t size = sizeof(stages)/sizeof(const TH1D *);
    TH1D * h = new TH1D("h1_cutflow", "cut flow", size, 0, size);
    // cout << " found " << size << " stages " << endl;
    for (UInt_t i = 0; i < size; i++) {
      // initialize bin
      h->Fill(stages[i], 0.);
    }
    // deflate -> remove unnecessary labels
    // h->LabelsDeflate();
    // Reset histogram statistics
    h->Reset();
    return h;
}
"""

def main():
    defaultOutputFileName = "CutFlow.C"
    usage = "usage: %prog [options] input-files"
    optParser = optparse.OptionParser(usage)
    optParser.add_option("-c", "--config", dest="cfgfile",
                         help="global configuration file",
                         default=ara.defaultConfigFileName)
    optParser.add_option("-f", "--file", dest="outputfile",
                         help="output file name",
                         default=defaultOutputFileName)
    (options, args) = optParser.parse_args()
    if len(args) < 1:
        optParser.print_help()
        return 1

    ara.config.read(options.cfgfile)

    # grep input files for cut flow
    files = " ".join(args)
    command = "grep 'Fill(\"cutflow\"' " + files  + " | cut -d'\"' -f4"
    stages = ara.getCommandOutput2(command)
    # insert C++ quotes 
    stages = '"' + '","'.join(stages.splitlines()) + '"'

    # open output file
    outputFile = open(options.outputfile, "w")
    outputFile.write(cutFlowTemplate % stages)
    outputFile.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
