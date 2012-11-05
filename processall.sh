#! /bin/bash

fullList=`ls Trees_261012/*/*/*.root`

for file in $fullList; do
    dirName=`dirname $file`
    fileName=`basename $file`
    echo "$dirName/$fileName"
    
    root -q -b runME_HZZ4l_FNAL.C\+\(\"$dirName/\",\"$fileName\",\"$dirName/\",-1,1\)
    
done

