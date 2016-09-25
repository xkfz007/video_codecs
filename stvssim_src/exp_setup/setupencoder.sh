#!/bin/bash
#yuvname="akiyo_qcif.yuv bridge_qcif.yuv carphone_qcif.yuv coastguard_qcif.yuv container_qcif.yuv football_qcif.yuv foreman_qcif.yuv grandma_qcif.yuv highway_qcif.yuv mobile_qcif.yuv news_qcif.yuv salesman_qcif.yuv silent_qcif.yuv soccer_qcif.yuv suzie_qcif.yuv" 
yuvname="football foreman news"
for i in $yuvname 
do
    dirname=${i}_cif
#    echo $dirname
    test -e $dirname || mkdir $dirname 
    cp -v encoder.cfg $dirname/
    cp -v lencod.exe $dirname/
done

