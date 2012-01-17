######################################################################
# You need to source this script in your .bashrc (or similar) in order
# to setup the ARA tools
#
# (C) 2011 Martin Weber

# prune_path takes two arguments. First is the actual part and second
# is the part of the path to get removed.
prune_path()
{
    # replace single occurence, at beginning, middle end,
    # and order such that no problem occurs
    echo "$1" | sed -e s@:$2:@:@g -e s@$2:@@g -e s@:$2@@g -e s@$2@@g
}


# prune old paths if they exist
if [[ -n $ARASYS ]] ; then
    export PATH=$(prune_path $PATH $ARASYS/bin)
    export PYTHONPATH=$(prune_path $PYTHONPATH $ARASYS/bin)
    export ARASYS=
fi

# set new path
DIR=$(dirname ${BASH_ARGV[0]})
export ARASYS=$(cd ${DIR}/..; pwd)

if [[ -z $PATH ]]; then
    export PATH=$ARASYS/bin
else
    export PATH=${PATH}:$ARASYS/bin
fi
if [[ -z $PYTHONPATH ]]; then
    export PYTHONPATH=$ARASYS/bin
else
    export PYTHONPATH=${PYTHONPATH}:$ARASYS/bin
fi
