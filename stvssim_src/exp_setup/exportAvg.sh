#!/bin/bash
arg=$1
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
    if [ "$arg" = "ssim" ]
    then
        newname=${i}"_avg_ssimrdo"
    elif [ "$arg" = "stvssim" ]
    then
        newname=${i}"_avg_stvssimrdo"
    else
        newname=${i}_avg
    fi
    #echo $newname
    #cat avgdata.txt
    echo "${newname}=[`cat avgdata.txt`];">>../$data
    cd -
done
