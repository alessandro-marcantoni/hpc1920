#!/bin/bash

for f in ace.in
do
	inputName=$f
	hullName=${f:0:-3}.hull
	pngName=${f:0:-3}.png

	echo "processing $f _________________________________________________________________________________"
	#echo $hullName
	#echo $pngName

	mpirun -n 7 mpi-hull < $inputName > $hullName
	gnuplot -c plot-hull.gp $inputName $hullName $pngName
done
