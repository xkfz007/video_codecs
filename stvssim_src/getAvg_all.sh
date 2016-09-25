#!/bin/bash
name=avgdata.txt

for qdir in *qcif
do
    test -d $qdir || continue
    cd $qdir
    test -f $name&&>$name||touch $name
    for i in avg_*.txt
    do
        dos2unix $i
        cat $i |sed -n  '9,17p'|awk '{print $4}'|tr -s "\n" " " >>$name
        bit=`cat $i|sed -n '19p'|awk '{print $8}'`
        echo ${bit%,}>>$name
        echo $i" -> "$name" done!"
    done
    cd ..
done
