#!/bin/bash

dlist=/data2/jmht/testset/testset_data/misc/files.list


i=0
sum=allresults.csv
for d in `cat $dlist`
do
    r=thresh/$d/AMPLE_0/benchmark/results.csv
    if [ $i = 0 ]; then
       cat $r | sed '/^\s*$/d' > $sum
    else
       cat $r | tail -n+2 | sed '/^\s*$/d' >> $sum
    fi
    (( i += 1 ))
done
