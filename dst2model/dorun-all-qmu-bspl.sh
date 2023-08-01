#!/bin/sh

rm mnew* mold* *matrix* AT* dT* qmu.sph
rm ??s*/A* ??s*/dTd* ??s*/raw* ??s*/mdcpl.out
rm ??s*/adpk.dat ??s*/fort*


spl=3
maxdeg=4

modes=(00s05 00s06 00s07 01s04 01s05 01s06 01s07 01s08 01s09 01s10 02s04 02s05 02s06 02s12 02s13 03s09) # allmodes

bash dorun-syn-qmu-bspl.sh $maxdeg "${modes[@]}"
bash dorun-deriv-qmu-bspl.sh $spl $maxdeg "${modes[@]}"
bash do-setup-inversion-qmu-bspl.sh $spl $maxdeg "${modes[@]}"
bash dorun-inversion-qmu-bspl.sh $spl $maxdeg "${modes[@]}"

cat inversion.dat 

./do-modelmap-depths-qmu.sh ${spl} ${spl}bspl

mv inversion-* models/
mv qmu_*sph models/
mv misfit* models/
mv model-xyz-* models
