#! /bin/bash

# before doing this, make sure you cleared up all the ME.root files..
# do this rm -f Trees_261012/*/*/*ME.root

export DIR=$1

if [ ! $# -eq 1 ]; then
	echo "USAGE: ./processall.sh DIR
	DIR - directory of the CJLST trees, e.g. Trees_261012/"
	exit 1
fi

rm -f $DIR/*/*/*ME.root

fullList=`ls $DIR/*/*/*.root`

for file in $fullList; do
    dirName=`dirname $file`
    fileName=`basename $file`
    echo "$dirName/$fileName"
    root -q -b runME_HZZ4l_FNAL.C\+\(\"$dirName/\",\"$fileName\",\"$dirName/\",-1,1\)
done

