#!/bin/bash
if [ ${#*} -le 2 ]; then
   echo "usage: addblock.sh pattern block file file1 [ file2 ... ]"
   echo
   echo "addblock.sh searches for the pattern in the given files and appends"
   echo "the file contents of the file block after the given line"
   exit 1
fi
SEDSCRIPT=`mktemp`
BLOCK=`sed -e '$ !s/$/\\\\/g' $2`
cat > $SEDSCRIPT <<EOF
# This sed script will append a block of lines after a given pattern
/$1/ a\\
$BLOCK
EOF
sed -i -f $SEDSCRIPT ${*:3}
rm $SEDSCRIPT
