#!/bin/sh

for d in $PWD/uvt*[0-9]
do
    echo $d
    cd $d
    x=$(ls -d V_*)
    rename ${x: -2} ${PWD: -2} *
done
