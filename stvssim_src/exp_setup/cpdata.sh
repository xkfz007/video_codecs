#!/bin/bash
yuvname1="akiyo_qcif.yuv bridge_qcif.yuv carphone_qcif.yuv coastguard_qcif.yuv"
yuvname3="football_qcif.yuv foreman_qcif.yuv grandma_qcif.yuv highway_qcif.yuv"
yuvname2="news_qcif.yuv salesman_qcif.yuv silent_qcif.yuv soccer_qcif.yuv" 
#dirnames="mb_stvssim_dataM1 mb_stvssim_dataM2 mb_stvssim_dataM3"
dirnames="mb_stvssimw_dataM"

for i in `seq 1 3`
do
   dirn=${dirnames}$i
  eval yuvn="$""yuvname"$i
   test -d $dirn || exit
   cd $dirn
   for j in $yuvn
   do
     tname=${j%.yuv}
	 echo $tname
      test -d "temp" || mkdir temp
	  cp -rv $tname/*.txt temp/
	  cp -rv temp  ~/develop/sciExprmnts/$dirnames/$tname
	done
    cd ..
done
