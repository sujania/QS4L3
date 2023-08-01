#!/bin/sh

spl=$1
maxdeg=$2
modes=("${@:3}")
maxcoeff=$(((${maxdeg}+1)*(${maxdeg}+2)/2))
params=$((${spl}*${maxcoeff}))

for dir in "${modes[@]}"
do
	echo $dir
	cd $dir
	smax=$(<smax.dat)
        if [ $(($smax%2)) -eq 0 ]; then
	rm A.dat
	../../bin/compu-deriv-dstmodel-bspl # derivatives
        ../../bin/buildATAcst << % 
$params
1.0
%
	fi
cd ..
done

