#!/bin/bash
######################################################################
# This script manages the extraction ("pick") of some few events     #
# from many datasets. Typically this is used for extracting events   #
# to be used for event scanning, e.g. with FireWorks                 #
#                                                                    #
# (C) 2011 Martin Weber                                              #
#                                                                    #
# v1.0                                                               #
######################################################################

function usage()
{
    NAME=`basename $0`
cat <<EOF
Pick some events from given datasets, e.g. for viewing with FireWorks

SYNOPSIS:  $NAME

OPTIONS:
  -h               This help message
  -e events.txt    Specify events in this file, one on each line.
                   use notation run:lumisection:eventnumber,
                   e.g. 140160:98:82878561
                   If not specified, defaults to "events.txt"
  -d datasets.txt  Specify the datasets in this file, one on each line,
                   e.g. /DoubleMu/Run2011A-PromptReco-v6/AOD
                   If not specified, defaults to "datasets.txt"
  -s               Submit the jobs
  -v               Retrieve job status
  -g               Get job output and merge all events
  -o merged.root   Output file name for fetch and merge operation
                   If not specified, defaults to "merged.root"
  -f               Force overwriting of existing files (use with -s or -g)

   Specifying both -s and -g does not make sense.

   You must have set your CMSSW environment and sourced your grid UI
   scripts before calling this program.

EXAMPLE:
   This script is run in three subsequent steps: 
   In a first step, submit the jobs: $NAME -s
   Check the status of the jobs: $NAME -v
   Once all jobs are DONE: $NAME -g

   The following is a cut-and-paste example for this script. Just copy
   the lines below and paste them in your shell. First, a file named
   datasets.txt is created that contains an example dataset file. Then an
   event file is created. Jobs are submitte, the status viewn. Once all
   jobs are completed (check manually), the files are retrieved and merged.

cat > datasets.txt <<EOD
/DoubleMu/Run2011A-PromptReco-v6/AOD
/DoubleMu/Run2011A-May10ReReco-v1/AOD
/DoubleMu/Run2011A-PromptReco-v4/AOD
/DoubleMu/Run2011A-05Aug2011-v1/AOD
/DoubleMu/Run2011B-PromptReco-v1/AOD
EOD
cat > events.txt <<EOD
172778:61:57905906
173389:28:34426455
173692:2162:2823842297
161008:65:40043881
EOD
# submit jobs
$NAME -d datasets.txt -e events.txt 
# look at job status
$NAME -v
# wait for jobs to finish
sleep 300
# look again at job status
$NAME -v
# once all jobs are finished
$NAME -g
EOF
}

function check_file()
{
    # file must have an absolut path
    FILE=$1
    DIRNAME=`dirname $FILE`
    if [[ "$DIRNAME" == "." ]]
    then
        FILE=$PWD/$FILE
    fi
    # test if file exists
    if [[ ! -e $FILE ]] 
    then
      echo "File $FILE not found"
      exit 2
    fi        
}

function check_environment()
{
    EXIT=
    # check CMSSW environment
    if [[ -z $CMSSW_VERSION ]]
    then
        echo "ERROR: You need to use \"cmsenv\" before invoking this script"
	EXIT=1
    fi

    # check CRAB environment
    if [[ -z $CRABLIBPYTHON ]] ; then
	echo "ERROR: You must setup your crab environment before invoking this script"
	EXIT=2
    fi
    if [[ -n $EXIT ]]; then
	exit $EXIT
    fi
}

# find association between datasets and events, which is unclear at this stage.
# an event could be found in multiple datasets or in only one, or even not at all.
# since the problem scales with number_of_datasets*number_of_events, submitting
# all events for all datasets is not a good idea. Therefore with this function
# we will find out which datasets need to be submitted. The resulting datasets
# will be stored in a cache file.
function find_assoc()
{
    if [[ -f $CACHE ]] ; then 
	echo "Regenerating datasets cache file..."
	rm $CACHE
    else
	echo "Generating datasets cache file..."
    fi
    touch $CACHE
    # check if there are some of the events in this dataset
    DSQUERY=
    for dataset in `cat $DATASETS`
    do
	if [[ -z $DSQUERY ]] ; then
	    DSQUERY="$dataset"
	else
	    DSQUERY="$DSQUERY,$dataset"
	fi
    done
    for event in `cat $EVENTS`
    do
	RUN=`echo $event | cut -d ':' -f 1`
	LS=`echo $event | cut -d ':' -f 2`
	EV=`echo $event | cut -d ':' -f 3`
	ID="Run: $RUN, Lumi section: $LS, event: $EV"
	RESFILE=`mktemp`
	dbsql "find dataset,file where dataset in ($DSQUERY) and run=$RUN and lumi=$LS" | grep "root" > $RESFILE
	NFOUND=`wc -l $RESFILE | cut -f 1 -d ' '`
	if [[ $NFOUND -eq 0 ]] ; then
	    echo "WARNING: Event $ID not found in any dataset"
	elif [[ $NFOUND -gt 1 ]] ; then
	    echo "WARNING: Event $ID found in multiple datasets ($NFOUND)"
	elif [[ $NFOUND -eq 1 ]] ; then
	    echo "Event $ID found"
	fi
	if [[ $NFOUND -gt 0 ]] ; then
	    cat $RESFILE | cut -f 1 -d ' ' >> $CACHE
	fi
	rm $RESFILE
    done
    # eliminate duplicates
    DSTMP=`mktemp`
    cat $CACHE | sort | uniq >> $DSTMP
    mv $DSTMP $CACHE
    echo "The following datasets will be processed:"
    cat $CACHE
}

function submit()
{
    echo "Submitting jobs (this will take some time)..."
    for dataset in `cat $CACHE`
    do
        DIR=`echo $dataset | sed -e s@/@_@g`
	if [[ -d $DIR ]] ; then
	    if [[ -n $FORCE ]] ; then
		rm -rf $DIR
	    else 
		cat <<EOF
ERROR: Directory already existing, $DIR
If you really want to submit, use the -f switch
or delete all previously created directories
EOF
		exit 3
	    fi
	fi
        mkdir -p $DIR
        cd $DIR
        edmPickEvents.py $dataset $EVENTS --crab
        crab -create -cfg pickevents_crab.config
        crab -submit
        cd ..
    done
}

function view()
{
    echo "Retrieving job status (this will take some time)..."
    for dataset in `cat $CACHE`
    do
        DIR=`echo $dataset | sed -e s@/@_@g`
        cd $DIR
	echo $DIR
	crab -status | grep "^[0-9].*"
	cd ..
    done
}

function fetch()
{
    echo "Fetching files from grid (This will take some time) ..."
    ROOTFILES=
    for dataset in `cat $CACHE`
    do
        DIR=`echo $dataset | sed -e s@/@_@g`
        cd $DIR
	# retrieve files from grid
        crab -getoutput
	# find root files and append to list
	FILES=`find . -name "*.root" | grep "/res/"`
	if [[ -n $FILES ]]
	then
	    for file in $FILES
	    do
		if [[ -z $ROOTFILES ]]
		then
		    ROOTFILES="file:$DIR/$file"
		else
 		    ROOTFILES="$ROOTFILES,file:$DIR/$file"
		fi
	    done
	fi
        cd ..
    done
    echo "Merging locally copied files..."
    # merge found files (if there are some...)
    if [[ -n $ROOTFILES ]]
    then
	edmCopyPickMerge inputFiles=$ROOTFILES outputFile=$OUTPUT maxSize=200000
    fi
}

function main()
{
    EVENTS="events.txt"
    DATASETS="datasets.txt"
    OUTPUT="merged.root"
    CACHE="pickEventsCache.txt"
    SUBMIT=
    VIEW=
    FETCH=
    while getopts "he:d:o:svgf" OPTION
    do
         case $OPTION in
             h)
                 usage
                 exit 1
                 ;;
             e)
                 EVENTS=$OPTARG
                 ;;
             d)
                 DATASETS=$OPTARG
                 ;;
             o)
                 OUTPUT=$OPTARG
                 ;;
             s)
                 SUBMIT=1
                 ;;
             v)
                 VIEW=1
                 ;;
             g)
                 FETCH=1
                 ;;
	     f)  FORCE=1
		 ;;
             ?)
                 usage
                 exit
                 ;;
         esac
    done
    
    if [[ -z $SUBMIT && -z $VIEW && -z $FETCH ]] || [[ -n $SUBMIT && -n $FETCH ]]
    then
         usage
         exit 1
    fi
    check_environment

    # check files
    if [[ -n $SUBMIT || -n $FETCH || -n $VIEW ]]
    then
        check_file $DATASETS
        DATASETS=$FILE
        echo "Using dataset file $DATASETS"
        check_file $EVENTS
        EVENTS=$FILE
        echo "Using event file $EVENTS"
	# check if files are newer than cache file
	if [[ $DATASETS -nt $CACHE ]] || [[ $EVENTS -nt $CACHE ]] ; then
	    find_assoc
	fi
    fi    
    if [[ -n $FETCH ]] 
    then
	if [[ -e $OUTPUT ]]
	then
	    if [[ -z $FORCE ]] ; then
		echo "Output file exists, please rename, move or remove it"
		echo "or use -f switch to overwrite"
		exit 4
	    else
		rm -f $OUTPUT
	    fi
	fi
	echo "Using output file $OUTPUT"
    fi

    # process commands
    if [[ -n $SUBMIT ]]
    then
	submit
    fi

    if [[ -n $VIEW ]]
    then
        view
    fi

    if [[ -n $FETCH ]]
    then
        fetch
    fi
}


# Call main 
main $@
