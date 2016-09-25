#!/bin/bash
name=avgdata.txt

for qdir in *cif
do
    test -d $qdir || continue
    cd $qdir
    test -f $name&&>$name||touch $name
    yuv_na=${qdir%_*}
    if [ "$yuv_na" = "bus" ]
    then
        for i in display_*.txt
        do
            dos2unix $i
            cat $i |sed -n  '114,116p'|awk '{print $11}'|tr -s ",\n" " " >>$name
            cat $i |sed -n  '117,119p'|awk '{print $4}'|tr -s "\n" " " >>$name
            bit=`cat $i|sed -n '121p'|awk '{print $8}'`
            echo ${bit%,}>>$name
            echo $i" -> "$name" done!"
        done
    else
        for i in display_*.txt
        do
            dos2unix $i
            cat $i |sed -n  '138,140p'|awk '{print $11}'|tr -s ",\n" " " >>$name
            cat $i |sed -n  '141,143p'|awk '{print $4}'|tr -s "\n" " " >>$name
            bit=`cat $i|sed -n '145p'|awk '{print $8}'`
            echo ${bit%,}>>$name
            echo $i" -> "$name" done!"
        done
    fi

    cd ..
done
