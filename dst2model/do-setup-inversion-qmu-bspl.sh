#!/bin/bash
spl=$1
maxdeg=$2
modes=("${@:3}")
maxcoeff=$(((${maxdeg}+1)*(${maxdeg}+2)/2))
params=$((${spl}*${maxcoeff}))

echo $params ../dta/mzero_${maxdeg}_l${spl}.in

cp ../dta/mzero_${maxdeg}_l${spl}.in qmu.sph

# works for any number of layers
../bin/sph2m-6l << % > mzero.dat
qmu.sph
%

cp mzero.dat mold.dat #mzero=mold

#addATAm input: number of model parameters, and number of mode
i=0
for mode in "${modes[@]}"; do
echo $i $mode
if (($i == 0)); then
cp ${mode}/ATA.dat . #initialize ATA.dat
cp ${mode}/ATd.dat . #initialize ATd.dat
cp ${mode}/dTd.dat . #initialize dTd.dat

echo $params > Ainput
echo $((${#modes[@]}-1)) >> Ainput # number of modes
else
echo ${mode}/ATA.dat >> Ainput
echo ${mode}/ATd.dat >> Ainput
echo ${mode}/dTd.dat >> Ainput
fi
i=$(($i+1))
done

../bin/addATAm < Ainput

cp ATA.dat ATAmatrix.dat
