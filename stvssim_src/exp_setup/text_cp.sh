#!/bin/bash
for i in s*rdo*
do
   test -d /cygdrive/i/$i||mkdir /cygdrive/i/$i;
   cd $i;
   for j in *qcif;
   do
     test -d /cygdrive/i/$i/$j||mkdir /cygdrive/i/$i/$j;
	 cp -v $j/*txt /cygdrive/i/$i/$j/;
   done
   cd ..;
done
