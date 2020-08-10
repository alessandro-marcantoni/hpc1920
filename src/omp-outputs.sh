#!/bin/bash

#echo -e "ace.in:" > outputs/omp-ace.txt
#./omp-convex-hull < ace.in >> outputs/omp-ace.txt
#echo -e "box10k.in:" > outputs/omp-box10k.txt
#./omp-convex-hull < box10k.in >> outputs/omp-box10k.txt
#echo -e "box100k.in:" > outputs/omp-box100k.txt
#./omp-convex-hull < box100k.in >> outputs/omp-box100k.txt
#echo -e "box1M.in:" > outputs/omp-box1M.txt
#./omp-convex-hull < box1M.in >> outputs/omp-box1M.txt
#echo -e "gaus100k.in:" > outputs/omp-gaus100k.txt
#./omp-convex-hull < gaus100k.in >> outputs/omp-gaus100k.txt
#echo -e "gaus1M.in:" > outputs/omp-gaus1M.txt
#./omp-convex-hull < gaus1M.in >> outputs/omp-gaus1M.txt
#echo -e "circ1k.in:" > outputs/omp-circ1k.txt
#./omp-convex-hull < circ1k.in >> outputs/omp-circ1k.txt
#echo -e "circ10k.in:" > outputs/omp-circ10k.txt
#./omp-convex-hull < circ10k.in >> outputs/omp-circ10k.txt
#echo -e "circ100k.in:" > outputs/omp-circ100k.txt
#./omp-convex-hull < circ100k.in >> outputs/omp-circ100k.txt

for f in *.in
do
	inputName=$f
	hullName=outputs/${f:0:-3}.hull
	
	echo -e $f > outputs/omp-${f:0:-3}.txt
	./omp-inner-hull < $f >> outputs/omp-${f:0:-3}.txt
done
