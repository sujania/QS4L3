#!/bin/sh

maxdeg=$1
modes=("${@:2}")
maxcoeff=$(((${maxdeg}+1)*(${maxdeg}+2)/2))

for dir in "${modes[@]}"
do
	cd $dir
	rm cstsyn.dat
	rm cstobs.dat
	rm cstweight.dat

	for (( c=1; c<=$maxcoeff; c++ )); do
		echo '0.000000' >> cstsyn.dat
	done

	smax=$(<smax.dat)
	smin=0
	echo "$dir" $maxdeg $maxcoeff $smax

	# copying uncert as weight
	cp cstobs_S20+CRUST.dat cstobs.dat
	awk '{print $2}' cstobs_S20+CRUST.dat > cstweight.dat #2nd kolom= max.error

	# filling the obs upto the maxdeg of the model
	if [ $smax -lt $maxdeg ]; then 
        for i in `seq $(($smax+2)) 2 $maxdeg`; do
        ind=$((2*$i+1))
        for j in `seq 1 1 $ind`; do
        echo "0.000000" >> cstobs.dat
        echo "1.000000" >> cstweight.dat
        done
	done
	fi
	sed -i '16,$d' cstobs.dat
	sed -i '16,$d' cstsyn.dat
	sed -i '16,$d' cstweight.dat
       cd ../
done

