#!/bin/bash
rm misfit.dat

if [ $# -lt 3 ]
then
    	printf "\n"
        printf "Give number of b-splines (bspl) and model smax\n"
        printf "USAGE: ./dorun-modemisfit2-nosome   #splines smax modes\n"
        printf "USAGE: ./dorun-modemisfit2-nosome  3 4 modes\n"
        exit
fi

spl=$1
maxdeg=$2
modes=("${@:3}")
maxcoeff=$(((${maxdeg}+1)*(${maxdeg}+2)/2))
params=$((${spl}*${maxcoeff}))

echo $maxcoeff $params

##loop through mode-folders and compute misfit
for dir in "${modes[@]}"; do
echo $dir
cd $dir

../../bin/mulAcst << %
$maxcoeff
$params
%

echo $dir >> ../misfit.dat
../../bin/misfit-permode2 << %
$maxcoeff
%

cd ../
done

../bin/avmisfit-permode2

