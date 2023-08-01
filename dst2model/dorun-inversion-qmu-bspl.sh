#!/bin/bash
rm inversion.dat
rm inversion-simplemisfit.dat
spl=$1
maxdeg=$2
modes=("${@:4}")

for DAMP in 1000000000 500000000 100000000 70000000 50000000 \
            20000000 10000000 9000000 7000000 \
            5000000 2000000 1000000 900000 700000 \
            500000 400000 300000 200000 100000 70000  \
            50000 20000 10000 1000 100 10 1 0.1 0.01 
do
rm inversion.out
rm misfit.dat
echo "damping" $DAMP

../bin/invATA-lapack-layer << %
$DAMP
16
30
%

# works for any number of layers
../bin/m2sph-6l << %
qmu.sph
%

#average misfit from misfit per mode
bash dorun-modemisfit2.sh $spl $maxdeg "${modes[@]}"
cp misfit.dat misfit-permode-qmu-${spl}bspl-${DAMP}.dat

cp mnew.dat mnew-$DAMP.dat
cp qmu.sph qmu_${spl}bspl_${DAMP}.sph
done

cp inversion.dat  inversion-qmu-${spl}bspl.dat
