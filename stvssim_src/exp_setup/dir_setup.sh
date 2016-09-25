#!/bin/bash
for i in lencod.exe*rdo*
do
    name=${i#lencod.exe_};
    echo $name
    test -d $name||mkdir $name
    cp -v *dll encoder.cfg run.m setupencoder.sh $name/
    cp -v $i $name/lencod.exe
    cd $name
    ./setupencoder.sh
    cd ..
done
