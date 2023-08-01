#!/bin/sh

# setting home directory
dir=`pwd`
echo 'HOME='$dir
sed -i "77s#.*#LIBDIR=${dir}#" src/libs/liblapack/make.inc

# setting dta pointers in derivatives coe
sed -i "14s#.*#      call chekcl('|-lu2:o:1:[${dir}/dta/foanis05.222]'#" src/compu-deriv-dstmodel-bspl.f
sed -i "16s#.*#     1          //'|-lu3:o:1:[${dir}/dta/m1084x2.htm] model on unit 3 (rdmdl)'#" src/compu-deriv-dstmodel-bspl.f
sed -i "17s#.*#     1          //'|-model:o:1:[${dir}/dta/mzero_4_l3.in] Model'#" src/compu-deriv-dstmodel-bspl.f
sed -i "144s#.*#      call openfl(1,'${dir}/dta/PREM222.BIN',1,0,0,istat,5364)#" src/compu-deriv-dstmodel-bspl.f

# compiling libraries
for lib in libandrew libgeo liblapack; do
cd $dir/src/libs/$lib
make
done

# compiling programs
cd $dir/src
for program in compu-deriv-dstmodel-bspl mulAcst \
               misfit-permode2 avmisfit-permode2 \
               invATA-lapack-layer m2sph-6l sph2m-6l bspl2v-3l \
               buildATAcst addATAm; do
echo $program
make $program
done
cd $dir
