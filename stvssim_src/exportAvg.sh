#!/bin/bash
data=avgdata_all.m
if [ -f "$data" ]
then
    >$data
else
    touch $data
fi
for i in *
do
    test -d $i || continue
    cd $i
    
    if [ ! -f "avgdata.txt" ]
    then
        echo "file avgdata.txt doesn't exist"
        exit
    fi
    test -f $data && \rm -v $data
    #cat avgdata.txt
    newname="${i}"
    echo "${newname}=[`cat avgdata.txt`];">>../$data
    cd -
done
