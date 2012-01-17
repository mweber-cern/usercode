#!/bin/bash
HOST=`uname -n`
echo "job_starter running on host $HOST"
echo "Current working dir is $PWD"
echo "Current directory contents:"
ls -l
EXE=$1
shift
ARGS=$*
echo "Executable is $EXE"
echo "Arguments are $ARGS"
echo "Starting executable now..."
$EXE $ARGS
RC=$?
echo "Executable finished with return code $RC"
echo "Current directory contents:"
ls -l
echo "job_starter done."
exit $RC
