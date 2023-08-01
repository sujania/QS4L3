#! /bin/bash
rm model-xyz.dat latlon.in

if [ $# -ne 2 ]
then
        printf "\n"
        printf "Parametrization number of layers/splines, model\n"
        printf "USAGE: ./do-modelmap-depths-qmu.sh 3 3bspl \n"
        exit
fi
layer=$1
spl=$2

ls qmu_${spl}*.sph

for mod in `ls qmu_${spl}*.sph`; do
for depth in 100 1500 2400; do

damp=${mod%.sph}

echo $mod $depth $layer $spl $damp

./do-sph2v-l $mod $depth $layer $spl #don't plot deg zero

mv model-xyz.dat model-xyz-$damp-$depth.dat

done
done


