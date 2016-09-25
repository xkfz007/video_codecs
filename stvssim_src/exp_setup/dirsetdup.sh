#!/bin/bash
#category="mserdo_lambda ssimrdo_lambda stvssimrdo_lambda"
category="mserdo ssimrdo stvssimrdo"
for i in $category
do
    dirname=${i}"_lambda"
    test -d $dirname || mkdir $dirname
    exe="lencod_"${i}".exe"
    cp -v $exe  $dirname/lencod.exe
    cp -v run.m encoder.cfg matlabRun.sh setupencoder.sh $dirname
	cd $dirname
	./setupencoder.sh
	cd ..
done

    
